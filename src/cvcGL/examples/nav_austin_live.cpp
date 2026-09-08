// nav_austin_live — a communications-denial-aware convoy on the real Austin bundle,
// running in real time on cvc::nav's sim_world (native drive: no Python, no libtorch).
//
// N vehicles drive a west->east route (global grid_nav::astar spine + learned CoefMLP
// local control, the same architecture as nav_finale). Over the run, time-gated RF
// JAM WINDOWS switch on and off: each is a jammer at a world position denying a disc
// of radius `radius_m` over [start_s, end_s) at some attenuation. While a window is
// active its disc is stamped into the A* PLANNING occupancy, so the convoy REROUTES
// around the denied zone in real time and straightens back when the window closes —
// the RF field bending the global plan, live. The overlay draws each active jammer as
// a translucent red disc + ring at its position, and the HUD reports the active
// channel and attenuation.
//
// The RF here is a deliberately simple flat-radius disc — the geometric FOOTPRINT of
// denial, not a path-loss field. A continuous native RF field solver (propagation +
// live solve) is a separate library on its own release cadence; this demo shows the
// nav RESPONSE to denial in real time and is self-contained.
//
// Jam windows come from a generic --rf-bundle JSON at runtime (nothing vendored):
//   { "jam_windows": [ { "label": "J1", "channel_id": "2",
//                        "start_s": 8, "end_s": 42, "attenuation_db": 70,
//                        "jammer_position_m": { "x": 120, "y": -80 },
//                        "radius_m": 420 }, ... ] }
// Flat "x"/"y"/"radius_m" keys are also accepted. With no --rf-bundle a synthetic
// jammer astride the convoy's route is generated so the demo stands alone.
//
//   nav_austin_live --bundle /path/to/scenes/austin_south --capture fly --offscreen \
//                   --frames 900 --out /tmp/austin_live [--rf-bundle jam.json]
//   ffmpeg -framerate 30 -i /tmp/austin_live/frame_%05d.png -c:v libx264 \
//          -pix_fmt yuv420p austin_live.mp4
// With no --bundle it falls back to a synthetic city so the scenario still runs.

#include "nav_common.h"

#include <algorithm>
#include <boost/program_options.hpp>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cvc/core/app.h>
#include <cvc/geometry/geometry.h>
#include <cvc/gl/CameraController.h>
#include <cvc/gl/GeometryNode.h>
#include <cvc/gl/ImGuiBinding.h>
#include <cvc/gl/ImGuiOverlay.h>
#include <cvc/gl/SceneGraph.h>
#include <cvc/gl/SceneRenderer.h>
#include <cvc/gl/ScreenTextHud.h>
#include <cvc/gl/TouchGestures.h>
#ifdef CVC_ENABLE_IMGUI
#include <imgui.h>
#endif
#include <cvc/image/image.h>
#include <cvc/nav/coef_mlp.h>
#include <cvc/nav/coef_train.h> // training_scene / city_scene (synthetic fallback)
#include <cvc/nav/grid_nav.h>   // inflate
#include <cvc/nav/sim_world.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <vector>
#include <vtkRenderer.h>

using cvc::gl::CameraController;
using cvc::gl::GeometryNode;
using cvc::gl::GeometryRenderMode;
using cvc::gl::SceneGraph;
using cvc::gl::SceneRenderer;

namespace {
const double PI = std::acos(-1.0);

// One RF jam window: a jammer denying a disc of radius `radius_m` centred at
// (x, y) world metres, active over [start_s, end_s) seconds of scenario time, at
// `atten_db` attenuation. `label`/`channel` are for the HUD only.
struct JamWindow {
  std::string label, channel;
  double start_s = 0, end_s = 0, atten_db = 0;
  double x = 0, y = 0, radius_m = 0;
  bool active(double t) const { return t >= start_s && t < end_s; }
};

// Pull the first JSON number that follows "key" (searched only within [from,to) of
// s). Returns false if the key or a number is not found in range.
bool json_num(const std::string &s, const char *key, double &v, std::size_t from = 0,
              std::size_t to = std::string::npos) {
  if (to == std::string::npos)
    to = s.size();
  const std::string k = std::string("\"") + key + "\"";
  const auto p = s.find(k, from);
  if (p == std::string::npos || p >= to)
    return false;
  const auto colon = s.find(':', p + k.size());
  if (colon == std::string::npos || colon >= to)
    return false;
  v = std::atof(s.c_str() + colon + 1);
  return true;
}

// Pull the first JSON string value that follows "key" within [from,to).
bool json_str(const std::string &s, const char *key, std::string &out, std::size_t from = 0,
              std::size_t to = std::string::npos) {
  if (to == std::string::npos)
    to = s.size();
  const std::string k = std::string("\"") + key + "\"";
  const auto p = s.find(k, from);
  if (p == std::string::npos || p >= to)
    return false;
  const auto colon = s.find(':', p + k.size());
  if (colon == std::string::npos || colon >= to)
    return false;
  const auto q1 = s.find('"', colon + 1);
  if (q1 == std::string::npos || q1 >= to)
    return false;
  const auto q2 = s.find('"', q1 + 1);
  if (q2 == std::string::npos || q2 > to)
    return false;
  out = s.substr(q1 + 1, q2 - q1 - 1);
  return true;
}

// Parse a generic jam-window bundle: the "jam_windows" (or "windows") array of
// objects. Each object supplies label/channel_id/start_s/end_s/attenuation_db and a
// position — either nested "jammer_position_m":{"x","y"} + "radius_m", or flat
// "x"/"y"/"radius_m". Objects missing a required numeric field are skipped. Returns
// the windows found (empty on any read/parse failure — the caller then synthesizes).
std::vector<JamWindow> read_jam_bundle(const std::string &path) {
  std::vector<JamWindow> out;
  std::ifstream f(path);
  if (!f)
    return out;
  const std::string s((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  auto arr = s.find("\"jam_windows\"");
  if (arr == std::string::npos)
    arr = s.find("\"windows\"");
  if (arr == std::string::npos)
    return out;
  const auto lb = s.find('[', arr);
  if (lb == std::string::npos)
    return out;
  // Walk top-level {...} objects inside the array (brace-depth 1 within the [ ]).
  int depth = 0;
  std::size_t objStart = std::string::npos;
  for (std::size_t i = lb; i < s.size(); ++i) {
    const char c = s[i];
    if (c == '[' && depth == 0 && i != lb)
      continue;
    if (c == ']' && depth == 0)
      break;
    if (c == '{') {
      if (depth == 0)
        objStart = i;
      ++depth;
    } else if (c == '}') {
      --depth;
      if (depth == 0 && objStart != std::string::npos) {
        const std::size_t a = objStart, b = i + 1;
        JamWindow w;
        double v = 0;
        const bool haveT = json_num(s, "start_s", w.start_s, a, b) &&
                           json_num(s, "end_s", w.end_s, a, b) &&
                           json_num(s, "attenuation_db", w.atten_db, a, b);
        // Position: nested jammer_position_m{x,y} first, else flat x/y.
        bool haveXY = false;
        const auto jp = s.find("\"jammer_position_m\"", a);
        if (jp != std::string::npos && jp < b) {
          const auto jb = s.find('}', jp);
          const std::size_t pe = (jb != std::string::npos && jb < b) ? jb + 1 : b;
          haveXY = json_num(s, "x", w.x, jp, pe) && json_num(s, "y", w.y, jp, pe);
        }
        if (!haveXY)
          haveXY = json_num(s, "x", w.x, a, b) && json_num(s, "y", w.y, a, b);
        if (!json_num(s, "radius_m", w.radius_m, a, b))
          w.radius_m = 300.0; // a sane default footprint if only the window is given
        json_str(s, "label", w.label, a, b);
        if (!json_str(s, "channel_id", w.channel, a, b)) {
          if (json_num(s, "channel_id", v, a, b))
            w.channel = std::to_string(static_cast<int>(v));
        }
        if (haveT && haveXY && w.end_s > w.start_s && w.radius_m > 0.0)
          out.push_back(w);
        objStart = std::string::npos;
      }
    }
  }
  return out;
}

// Stamp every ACTIVE jammer's disc into a COPY of the base planning occupancy as
// blocked cells, so A* routes around the denied zone. grid convention: r->y (row 0
// = min_y), c->x, index [r*cols + c] — matching occupancy_from_model / sim_world.
std::vector<std::uint8_t> stamp_jammers(const std::vector<std::uint8_t> &base, int rows, int cols,
                                        const navdemo::Bounds &b,
                                        const std::vector<JamWindow> &jams, double t) {
  std::vector<std::uint8_t> occ = base;
  const double sx = (cols - 1) / (b.max_x - b.min_x);
  const double sy = (rows - 1) / (b.max_y - b.min_y);
  for (const auto &j : jams) {
    if (!j.active(t))
      continue;
    // Denied radius in cells (use the larger axis scale for a conservative disc).
    const double rc = j.radius_m * std::max(sx, sy);
    const int cc = static_cast<int>(std::lround((j.x - b.min_x) * sx));
    const int cr = static_cast<int>(std::lround((j.y - b.min_y) * sy));
    const int ri = static_cast<int>(std::ceil(rc));
    for (int r = std::max(0, cr - ri); r <= std::min(rows - 1, cr + ri); ++r)
      for (int c = std::max(0, cc - ri); c <= std::min(cols - 1, cc + ri); ++c) {
        const double dr = r - cr, dc = c - cc;
        if (dr * dr + dc * dc <= rc * rc)
          occ[static_cast<std::size_t>(r) * cols + c] = 1;
      }
  }
  return occ;
}

// A signature of the active-jammer set at time t (which windows are on), so the
// route planner only replans when the denied region actually changes.
unsigned long long active_mask(const std::vector<JamWindow> &jams, double t) {
  unsigned long long m = 0;
  for (std::size_t i = 0; i < jams.size() && i < 64; ++i)
    if (jams[i].active(t))
      m |= (1ull << i);
  return m;
}
} // namespace

int main(int argc, char **argv) {
  namespace po = boost::program_options;
  std::string bundle, vehicle, rfBundle, capture = "fly", out = "frames", png, viewMode = "3d";
  int width = 1280, height = 720, nx = 384, N = 6;
  long frames = 0;
  double mouseSens = 0.25, moveSpeed = 0.0;
  double fps = 30.0, hz = 60.0;
  double standoffCells = 4.0;
  bool offscreen = false, no_shadows = false, no_minimap = false, no_ui = false;

  po::options_description desc("nav_austin_live — a comms-denial-aware convoy on Austin, in cvcGL");
  desc.add_options()("help,h", "show this help")(
      "bundle", po::value<std::string>(&bundle), "Austin bundle dir (terrain.json + buildings.glb)")(
      "rf-bundle", po::value<std::string>(&rfBundle),
      "jam-window JSON (generic; synthesized if omitted)")(
      "vehicles", po::value<int>(&N)->default_value(6), "convoy size")(
      "grid", po::value<int>(&nx)->default_value(384), "occupancy resolution")(
      "vehicle", po::value<std::string>(&vehicle),
      "vehicle model .glb (default: <bundle>/../../shared/Humvee.glb; else an arrow)")(
      "offscreen", po::bool_switch(&offscreen))("no-shadows", po::bool_switch(&no_shadows))(
      "no-minimap", po::bool_switch(&no_minimap), "hide the 2-D PiP minimap")(
      "frames", po::value<long>(&frames)->default_value(0))(
      "fps", po::value<double>(&fps)->default_value(30.0))(
      "hz", po::value<double>(&hz)->default_value(60.0))(
      "view", po::value<std::string>(&viewMode)->default_value("3d"),
      "3d (perspective) | map (top-down ortho)")(
      "capture", po::value<std::string>(&capture)->default_value("none"),
      "none (interactive) | fly | orbit (offscreen PNG capture)")(
      "mouse-sensitivity", po::value<double>(&mouseSens)->default_value(0.25))(
      "move-speed", po::value<double>(&moveSpeed)->default_value(0.0))(
      "width", po::value<int>(&width)->default_value(1280))(
      "height", po::value<int>(&height)->default_value(720))(
      "out", po::value<std::string>(&out)->default_value("frames"))(
      "png", po::value<std::string>(&png))("no-ui", po::bool_switch(&no_ui),
                                           "hide the ImGui overlay");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 0;
  }
  if (N < 1)
    N = 1;
  // Auto-probe for a scene bundle when none is given: $CVC_NAV_BUNDLE, then
  // ~/scenes/austin_south. The provenance goes on the HUD either way.
  if (bundle.empty()) {
    const char *envB = std::getenv("CVC_NAV_BUNDLE");
    if (envB && std::filesystem::exists(std::string(envB) + "/terrain.json"))
      bundle = envB;
    else if (const char *home = std::getenv("HOME")) {
      const std::string cand = std::string(home) + "/scenes/austin_south";
      if (std::filesystem::exists(cand + "/terrain.json"))
        bundle = cand;
    }
  }

  const bool capturing = (capture != "none");
  if (capturing) {
    offscreen = true;
    if (frames <= 0)
      frames = 900;
    std::filesystem::create_directories(out);
    std::printf("nav_austin_live: capturing %ld frames (%s, offscreen) -> %s/frame_*.png\n", frames,
                capture.c_str(), out.c_str());
  }
  bool minimap = !no_minimap;
  bool orthoNow = false;

  // 1. Occupancy + a mesh to render. Real bundle if given; else a synthetic city.
  navdemo::Bounds bounds;
  int ny = nx;
  std::vector<std::uint8_t> occ;
  cvc::geometry cityMesh;
  navdemo::Terrain terrain;
  bool haveMesh = false, usedBundle = false;

  if (!bundle.empty() &&
      navdemo::load_city_bundle(bundle, nx, ny, bounds, occ, &cityMesh, &terrain, /*prefer_flat=*/false)) {
    haveMesh = cityMesh.num_tris() > 0;
    usedBundle = true;
    std::printf("nav_austin_live: bundle %s  bounds [%.0f,%.0f]..[%.0f,%.0f]  %llu tris\n",
                bundle.c_str(), bounds.min_x, bounds.min_y, bounds.max_x, bounds.max_y,
                (unsigned long long)cityMesh.num_tris());
  }
  if (occ.empty()) { // synthetic fallback
    cvc::nav::training_scene ts = cvc::nav::city_scene(std::min(nx, 128));
    nx = ts.cols;
    ny = ts.rows;
    occ = ts.occ;
    bounds = {ts.min_x, ts.min_y, ts.max_x, ts.max_y};
    std::printf("nav_austin_live: no bundle; synthetic city %dx%d\n", nx, ny);
  }
  navdemo::add_border(occ.data(), ny, nx);
  const double span = std::max(bounds.max_x - bounds.min_x, bounds.max_y - bounds.min_y);

  // 2. The convoy: N vehicles entering from the west edge (spread in y), all bound
  //    for a shared rendezvous on the eastern side.
  const double sc = 0.05; // world->normalized (1 unit = 20 m); the tuned vehicle scale.
  auto norm = [&](double wx, double wy, float *o2) {
    o2[0] = static_cast<float>((wx - 0.5 * (bounds.min_x + bounds.max_x)) * sc);
    o2[1] = static_cast<float>((wy - 0.5 * (bounds.min_y + bounds.max_y)) * sc);
  };
  const double westX = bounds.min_x + 0.10 * (bounds.max_x - bounds.min_x);
  const double eastX = bounds.min_x + 0.88 * (bounds.max_x - bounds.min_x);
  const double midY = 0.5 * (bounds.min_y + bounds.max_y);
  std::vector<float> o(2 * N), goal(2 * N), color(3 * N);
  for (int i = 0; i < N; ++i) {
    const double frac = N > 1 ? static_cast<double>(i) / (N - 1) : 0.5;
    const double y = bounds.min_y + (0.40 + 0.20 * frac) * (bounds.max_y - bounds.min_y);
    norm(westX, y, &o[2 * i]);
    norm(eastX, midY, &goal[2 * i]);
    // One convoy, one hue family (a warm amber), lightening slightly down the column.
    const double lift = 0.12 * frac;
    color[3 * i] = static_cast<float>(std::min(1.0, 0.90 + lift));
    color[3 * i + 1] = static_cast<float>(std::min(1.0, 0.66 + lift));
    color[3 * i + 2] = static_cast<float>(std::min(1.0, 0.28 + lift));
  }

  cvc::nav::sim_world::config cfg;
  cfg.rows = ny;
  cfg.cols = nx;
  cfg.min_x = bounds.min_x;
  cfg.min_y = bounds.min_y;
  cfg.max_x = bounds.max_x;
  cfg.max_y = bounds.max_y;
  cfg.cx = 0.5 * (bounds.min_x + bounds.max_x);
  cfg.cy = 0.5 * (bounds.min_y + bounds.max_y);
  cfg.scale = sc;
  cfg.veh.rr = 0.15f;
  cfg.veh.d_hat = 0.35f;
  cfg.veh.dt = 0.06f;
  cfg.veh.vmax = 0.9f;
  cfg.freeze_sense = true; // the convoy knows the city map
  cfg.range_m = 200.0;
  cfg.n_rays = 200;
  cfg.reach_tol = 0.8f;

  std::vector<int> map_id(N);
  for (int i = 0; i < N; ++i)
    map_id[i] = i;
  cvc::nav::sim_world world(cfg, occ.data(), occ.data(), cvc::nav::coef_mlp::default_biased(),
                            o.data(), goal.data(), color.data(), N, map_id.data(), N);

  // 3. Jam windows. From --rf-bundle if given/parseable; else a synthetic pair astride
  //    the convoy's route so the demo stands alone. Synthetic timing spans the run.
  std::vector<JamWindow> jams;
  if (!rfBundle.empty()) {
    jams = read_jam_bundle(rfBundle);
    std::printf("nav_austin_live: rf-bundle %s -> %zu jam window(s)\n", rfBundle.c_str(),
                jams.size());
  }
  if (jams.empty()) {
    const double runS = (frames > 0 ? frames : 900) / fps;
    // A primary jammer straddling the route centre, on for the middle of the run, and
    // a second nearer the goal that flicks on late — two distinct reroutes.
    JamWindow a;
    a.label = "ALPHA";
    a.channel = "2";
    a.start_s = 0.12 * runS;
    a.end_s = 0.62 * runS;
    a.atten_db = 72;
    a.x = cfg.cx + 0.04 * span;
    a.y = midY;
    a.radius_m = 0.16 * span;
    JamWindow b;
    b.label = "BRAVO";
    b.channel = "5";
    b.start_s = 0.50 * runS;
    b.end_s = 0.92 * runS;
    b.atten_db = 64;
    b.x = bounds.min_x + 0.70 * (bounds.max_x - bounds.min_x);
    b.y = midY + 0.10 * span;
    b.radius_m = 0.12 * span;
    jams = {a, b};
    std::printf("nav_austin_live: synthetic jam windows: %zu\n", jams.size());
  }

  // 4. Scene.
  cvc::app app;
  app.properties("system.log_verbosity", "2");
  SceneGraph sg(app, "austin_live");

  if (haveMesh) {
    auto b = std::dynamic_pointer_cast<GeometryNode>(sg.addGraphics("buildings", cityMesh));
    if (b) {
      b->setUseSingleColor(true);
      b->setColor(0.62, 0.62, 0.66);
      b->setAmbient(0.4);
      b->setDiffuse(0.8);
    }
  } else {
    const double wall_rgb[3] = {0.58, 0.58, 0.64};
    sg.addGraphics("buildings", navdemo::occupancy_to_walls(occ.data(), ny, nx, bounds, 0.05 * span,
                                                            wall_rgb, /*vary=*/0.45));
  }
  const double ground_rgb[3] = {0.18, 0.21, 0.25};
  auto groundNode = std::dynamic_pointer_cast<GeometryNode>(
      sg.addGraphics("ground", navdemo::ground_quad(bounds, 0.0, ground_rgb)));
  cvc::image sat;
  bool haveSat = false;
  if (usedBundle) {
    const std::string satpath = bundle + "/satellite.png";
    if (std::filesystem::exists(satpath)) {
      sat = cvc::read_image(satpath);
      haveSat = sat.width() > 0 && sat.height() > 0;
    }
  }
  if (groundNode) {
    if (haveSat) {
      groundNode->setUseSingleColor(true);
      groundNode->setColor(1.0, 1.0, 1.0);
      groundNode->setAmbient(0.85);
      groundNode->setDiffuse(0.5);
      groundNode->setTexture(sat, /*zeroCopy=*/false);
    } else {
      groundNode->setAmbient(0.6);
      groundNode->setDiffuse(0.5);
    }
  }

  navdemo::AgentGlyphs glyphs;
  cvc::geometry agentGeom;
  std::vector<double> hv;
  std::vector<std::uint32_t> hvt;
  std::string vpath = vehicle;
  if (vpath.empty() && usedBundle)
    vpath = bundle + "/../../shared/Humvee.glb";
  const bool haveVehicle =
      !vpath.empty() && navdemo::load_vehicle_template(vpath, 15.0, hv, hvt);
  if (haveVehicle)
    agentGeom = glyphs.build_template(app, N, color.data(), hv, hvt, 0.002 * span);
  else
    agentGeom = glyphs.build(app, N, color.data(), 0.02 * span, 0.006 * span);
  auto agentNode = std::dynamic_pointer_cast<GeometryNode>(sg.addGraphics("agents", agentGeom));
  if (agentNode) {
    agentNode->setUseSingleColor(false);
    agentNode->setAmbient(0.7);
    agentNode->setDiffuse(0.8);
  }

  // RF jammer overlay: a translucent red disc per window, sitting just above the
  // ground. Visibility is toggled by the active window each frame; the disc alpha
  // rides the attenuation. disc_marker gives a flat triangle-fan disc.
  std::vector<std::shared_ptr<GeometryNode>> jamDisc(jams.size());
  for (std::size_t k = 0; k < jams.size(); ++k) {
    const double denied_rgb[3] = {0.85, 0.12, 0.10};
    auto d = std::dynamic_pointer_cast<GeometryNode>(sg.addGraphics(
        "jam_disc_" + std::to_string(k), navdemo::disc_marker(jams[k].radius_m, 2.0, denied_rgb)));
    if (d) {
      d->setUseSingleColor(true);
      d->setColor(0.85, 0.12, 0.10);
      d->setAmbient(1.0);
      d->setDiffuse(0.0);
      d->setOpacity(0.28);
      d->setDepthOffset(3.0);
      d->setPosition(jams[k].x, jams[k].y, 2.0);
      d->setVisible(false);
    }
    jamDisc[k] = d;
  }

  // A* route spines, one per vehicle, streamed as tube LINES in the convoy hue.
  auto make_lines_node = [&](const std::string &name, int segs, const double rgb[3], double lw,
                             double z) {
    cvc::geometry lg;
    for (int k = 0; k < segs; ++k) {
      for (int e = 0; e < 2; ++e) {
        lg.points().push_back({cfg.cx, cfg.cy, z});
        lg.colors().push_back({rgb[0], rgb[1], rgb[2]});
      }
      lg.lines().push_back({static_cast<cvc::geometry::index_t>(2 * k),
                            static_cast<cvc::geometry::index_t>(2 * k + 1)});
    }
    auto node = std::dynamic_pointer_cast<GeometryNode>(sg.addGraphics(name, lg));
    if (node) {
      node->setRenderMode(GeometryRenderMode::LINES);
      node->setUseSingleColor(true);
      node->setColor(rgb[0], rgb[1], rgb[2]);
      node->setLineWidth(lw);
      node->setRenderLinesAsTubes(true);
      node->setAmbient(1.0);
      node->setDiffuse(0.0);
      node->setDepthOffset(2.0);
    }
    return node;
  };
  const int SPINE_SEGS = 64;
  std::vector<std::shared_ptr<GeometryNode>> spineNodes(N);
  std::vector<std::vector<double>> spineXyz(N);
  for (int i = 0; i < N; ++i) {
    const double rgb[3] = {0.85 * color[3 * i], 0.85 * color[3 * i + 1], 0.85 * color[3 * i + 2]};
    spineNodes[i] = make_lines_node("spine_" + std::to_string(i), SPINE_SEGS, rgb, 3.0, 3.0);
    spineXyz[i].assign(static_cast<std::size_t>(3) * 2 * SPINE_SEGS, 0.0);
  }

  sg.addDirectionalLight(-40, 58, 1.0, 0.96, 0.88, 1.1);
  sg.addDirectionalLight(150, 32, 0.5, 0.58, 0.72, 0.45);

  SceneRenderer view(sg, width, height, offscreen, "main");
  const bool shadows = !no_shadows && sg.setShadowsEnabled(true);
  if (shadows) {
    sg.setShadowResolution(2048);
    sg.setShadowUpdateInterval(capturing ? 2 : 4);
  }
  view.renderer()->GradientBackgroundOn();
  view.renderer()->SetBackground(0.20, 0.26, 0.36);
  view.renderer()->SetBackground2(0.55, 0.66, 0.82);

  cvc::gl::ScreenTextHud banner(view, "banner");
  banner.setPosition(0.5, 0.92);
  banner.setFontSize(20);
  banner.setColor(1.0, 0.55, 0.45);
  cvc::gl::ScreenTextHud status(view, "status");
  status.setCentered(false);
  status.setPosition(0.015, 0.03);
  status.setFontSize(13);
  status.setColor(0.85, 0.88, 0.92);
  {
    char st[240];
    std::snprintf(st, sizeof st,
                  "%d-vehicle convoy · %s · global A* spine + learned CoefMLP drive · "
                  "reroutes around RF denial live",
                  N, usedBundle ? "real Austin (runtime bundle)" : "synthetic city — pass --bundle");
    status.setText(st);
  }

  CameraController cam(view);
  cam.frameBounds(bounds.min_x, bounds.min_y, 0.0, bounds.max_x, bounds.max_y, 0.05 * span);
  cam.setMouseSensitivity(mouseSens);
  if (moveSpeed > 0.0)
    cam.setMoveSpeed(moveSpeed);
  if (!capturing && viewMode == "map")
    navdemo::set_ortho_topdown(view, bounds, 10.0, &cam);
  sg.setDiagnosticChromeVisible(false);
  sg.processEvents();

  // Minimap base: dimmed satellite (or occupancy grey), north-up.
  const int MM = 180;
  cvc::image mmbase(MM, MM, cvc::image::pixel_format::RGB, cvc::image::data_type::u8);
  {
    unsigned char *p = mmbase.data();
    const unsigned char *sd = haveSat ? sat.data() : nullptr;
    const int sw = haveSat ? sat.width() : 0, sh = haveSat ? sat.height() : 0,
              schn = haveSat ? sat.channels() : 0;
    for (int r = 0; r < MM; ++r)
      for (int c = 0; c < MM; ++c) {
        const long i = (static_cast<long>(r) * MM + c) * 3;
        if (haveSat) {
          const long si =
              (static_cast<long>(r * (sh - 1) / (MM - 1)) * sw + c * (sw - 1) / (MM - 1)) * schn;
          p[i] = static_cast<unsigned char>(0.4 * sd[si]);
          p[i + 1] = static_cast<unsigned char>(0.4 * sd[si + 1]);
          p[i + 2] = static_cast<unsigned char>(0.4 * sd[si + 2]);
        } else {
          const int oc = c * (nx - 1) / (MM - 1), orr = (MM - 1 - r) * (ny - 1) / (MM - 1);
          const unsigned char v = occ[static_cast<std::size_t>(orr) * nx + oc] ? 90 : 32;
          p[i] = v;
          p[i + 1] = v + 6;
          p[i + 2] = v + 12;
        }
      }
  }

  // 5. Run. Sim on the render thread (a small convoy). The base planning occupancy is
  //    a lightly-inflated copy of the truth; each time the ACTIVE-jammer set changes we
  //    stamp the active discs into it and replan every route (the live RF->nav coupling).
  const std::vector<std::uint8_t> planBase = cvc::nav::inflate(occ.data(), ny, nx, 1);
  std::vector<navdemo::Route> routes(N);
  std::vector<std::size_t> curWp(N, static_cast<std::size_t>(-1));
  std::vector<float> pos(2 * N), head(N), spd(N);
  std::vector<int> md(N);
  std::vector<std::uint8_t> rch(N);
  world.snapshot(pos.data(), head.data(), spd.data(), md.data(), rch.data());
  unsigned long long lastMask = ~0ull; // force an initial plan on frame 0

  const long total = frames > 0 ? frames : 900;
  const int sub = std::max(1, static_cast<int>(hz / std::max(1.0, fps)));
  double eye[3], focal[3];
  long frame = 0;
  const auto tw0 = std::chrono::steady_clock::now();
  double lastT = 0.0;

  bool uiPaused = false, ui2D = false, uiSpines = true, uiMinimap = minimap;
#ifdef CVC_ENABLE_IMGUI
  cvc::gl::ImGuiOverlay ui(view);
  ui.attachCamera(cam);
  cvc::gl::TouchGestures touch(view, cam);
  ui.setVisible(!no_ui && !capturing);
  int uiActive = 0;
  ImGui::SetCurrentContext(ui.imguiContext());
  ui.setDrawCallback([&] {
    if (ImGui::BeginMainMenuBar()) {
      if (ImGui::BeginMenu("Sim")) {
        ImGui::MenuItem("Paused", nullptr, &uiPaused);
        ImGui::EndMenu();
      }
      if (ImGui::BeginMenu("View")) {
        ImGui::MenuItem("2-D map", nullptr, &ui2D);
        ImGui::MenuItem("Route spines", nullptr, &uiSpines);
        ImGui::MenuItem("Minimap", nullptr, &uiMinimap);
        ImGui::EndMenu();
      }
      ImGui::EndMainMenuBar();
    }
    ImGui::SetNextWindowPos(ImVec2(10, 30), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(320, 0), ImGuiCond_FirstUseEver);
    ImGui::Begin("Austin live");
    ImGui::Text("%d vehicles | %s", N, usedBundle ? "real Austin" : "synthetic city");
    ImGui::Text("frame %ld | jammers active: %d", frame, uiActive);
    if (!usedBundle)
      ImGui::TextDisabled("pass --bundle for the real Austin map");
    ImGui::End();
  });
#else
  cvc::gl::TouchGestures touch(view, cam);
#endif

  while (!view.windowClosed()) {
    double t, dt;
    if (capturing) {
      t = frame / fps;
      dt = 1.0 / fps;
    } else {
      t = std::chrono::duration<double>(std::chrono::steady_clock::now() - tw0).count();
      dt = t - lastT;
      lastT = t;
    }
#ifdef CVC_ENABLE_IMGUI
    for (int i = 0; i < N; ++i)
      if (spineNodes[i])
        spineNodes[i]->setVisible(uiSpines);
    if (ui2D != orthoNow) {
      orthoNow = ui2D;
      if (orthoNow)
        navdemo::set_ortho_topdown(view, bounds, 10.0, &cam);
      else {
        cam.setMode(CameraController::Mode::Orbit);
        cam.frameBounds(bounds.min_x, bounds.min_y, 0.0, bounds.max_x, bounds.max_y, 0.05 * span);
      }
    }
    minimap = uiMinimap;
    if (uiPaused) {
      view.processUIEvents();
      touch.update();
      view.render();
      continue;
    }
#endif

    // RF: update the active-jammer set. Replan only when it changes (a window opens or
    // closes) — replanning every frame would reset each carrot FSM's escape state.
    const unsigned long long mask = active_mask(jams, t);
    int nActive = 0;
    for (std::size_t k = 0; k < jams.size(); ++k) {
      const bool on = jams[k].active(t);
      nActive += on;
      // Pulse the disc alpha a little while active so it reads as "live".
      const double pulse = 0.24 + 0.08 * std::sin(6.0 * t);
      if (jamDisc[k]) {
        jamDisc[k]->setVisible(on);
        if (on)
          jamDisc[k]->setOpacity(pulse * std::min(1.0, jams[k].atten_db / 70.0));
      }
    }
#ifdef CVC_ENABLE_IMGUI
    uiActive = nActive;
#endif
    if (mask != lastMask) {
      lastMask = mask;
      const std::vector<std::uint8_t> planOcc = stamp_jammers(planBase, ny, nx, bounds, jams, t);
      for (int i = 0; i < N; ++i) {
        const double gx = goal[2 * i] / sc + cfg.cx, gy = goal[2 * i + 1] / sc + cfg.cy;
        routes[i] = navdemo::plan_route(planOcc.data(), ny, nx, bounds, pos[2 * i], pos[2 * i + 1],
                                        gx, gy, /*inflate_cells=*/0, standoffCells);
        curWp[i] = static_cast<std::size_t>(-1);
      }
    }

    // Follow the route: advance waypoint on arrival; retarget only when it changes.
    for (int i = 0; i < N; ++i) {
      navdemo::Route &rt = routes[i];
      const double arr = 0.03 * span;
      while (rt.idx + 1 < rt.wp.size()) {
        const double dx = rt.wp[rt.idx][0] - pos[2 * i], dy = rt.wp[rt.idx][1] - pos[2 * i + 1];
        if (dx * dx + dy * dy < arr * arr)
          ++rt.idx;
        else
          break;
      }
      if (!rt.wp.empty() && rt.idx != curWp[i]) {
        float w2[2];
        norm(rt.wp[rt.idx][0], rt.wp[rt.idx][1], w2);
        world.retarget(i, w2[0], w2[1]);
        curWp[i] = rt.idx;
      }
    }

    for (int s = 0; s < sub; ++s)
      world.step();
    world.snapshot(pos.data(), head.data(), spd.data(), md.data(), rch.data());
    const auto &xyz = glyphs.pack(pos.data(), head.data());
    if (agentNode)
      agentNode->updateVertices(xyz);

    if (frame % 5 == 0)
      for (int i = 0; i < N; ++i) {
        if (!spineNodes[i])
          continue;
        std::vector<double> &X = spineXyz[i];
        const navdemo::Route &rt = routes[i];
        double ax = pos[2 * i], ay = pos[2 * i + 1];
        std::size_t sseg = 0;
        for (std::size_t w = rt.idx; w < rt.wp.size() && sseg < SPINE_SEGS; ++w, ++sseg) {
          const std::size_t b6 = 6 * sseg;
          X[b6] = ax;
          X[b6 + 1] = ay;
          X[b6 + 2] = 3.0;
          X[b6 + 3] = rt.wp[w][0];
          X[b6 + 4] = rt.wp[w][1];
          X[b6 + 5] = 3.0;
          ax = rt.wp[w][0];
          ay = rt.wp[w][1];
        }
        for (; sseg < SPINE_SEGS; ++sseg) {
          const std::size_t b6 = 6 * sseg;
          X[b6] = X[b6 + 3] = ax;
          X[b6 + 1] = X[b6 + 4] = ay;
          X[b6 + 2] = X[b6 + 5] = 3.0;
        }
        spineNodes[i]->updateVertices(X);
      }

    // Banner: name the active denial. Quiet when clear.
    if (nActive > 0) {
      std::string act;
      double maxdb = 0;
      for (const auto &j : jams)
        if (j.active(t)) {
          if (!act.empty())
            act += ", ";
          act += (j.label.empty() ? std::string("?") : j.label);
          if (!j.channel.empty())
            act += " (ch " + j.channel + ")";
          maxdb = std::max(maxdb, j.atten_db);
        }
      char bn[200];
      std::snprintf(bn, sizeof bn, "RF DENIAL ACTIVE  \xe2\x80\x94  %s  \xe2\x88\x92%.0f dB  \xe2\x80\x94  convoy rerouting",
                    act.c_str(), maxdb);
      banner.setText(bn);
    } else {
      banner.setText("");
    }

    if (capturing) {
      const double az = 0.4 + 0.25 * std::sin(0.12 * t);
      navdemo::orbit_camera(bounds, 0.05 * span, az, 42.0 * PI / 180.0, 1.9, eye, focal);
      view.setCamera(eye[0], eye[1], eye[2], focal[0], focal[1], focal[2], 0, 0, 1, 30);
    } else {
      cam.update(dt);
    }

    if (capturing) {
      view.render();
      std::vector<unsigned char> rgb = view.frameRGB();
      const int W = view.frameWidth(), H = view.frameHeight();
      cvc::image frameImg(W, H, cvc::image::pixel_format::RGB, cvc::image::data_type::u8);
      unsigned char *fp = frameImg.data();
      for (int r = 0; r < H; ++r)
        std::copy(&rgb[(static_cast<long>(H - 1 - r) * W) * 3],
                  &rgb[(static_cast<long>(H - r) * W) * 3], &fp[static_cast<long>(r) * W * 3]);
      if (minimap) {
        const int pad = 12, x0 = W - MM - pad, y0 = pad;
        navdemo::blit_clamped(fp, W, H, mmbase.data(), MM, MM, 3, x0, y0);
        auto toPx = [&](double wx, double wy, int &px, int &py) {
          px = x0 + static_cast<int>((wx - bounds.min_x) / (bounds.max_x - bounds.min_x) * (MM - 1));
          py = y0 + static_cast<int>((bounds.max_y - wy) / (bounds.max_y - bounds.min_y) * (MM - 1));
        };
        // Active jammers: a red disc on the minimap.
        for (const auto &j : jams)
          if (j.active(t)) {
            int jx, jy;
            toPx(j.x, j.y, jx, jy);
            const int rr = static_cast<int>(j.radius_m / (bounds.max_x - bounds.min_x) * (MM - 1));
            navdemo::plot_disc(fp, W, H, jx, jy, std::max(2, rr), 200, 40, 36);
          }
        for (int i = 0; i < N; ++i) {
          const navdemo::Route &rt = routes[i];
          const unsigned char lr = static_cast<unsigned char>(90 + color[3 * i] * 165);
          const unsigned char lg = static_cast<unsigned char>(90 + color[3 * i + 1] * 165);
          const unsigned char lb = static_cast<unsigned char>(90 + color[3 * i + 2] * 165);
          int ax_, ay_;
          toPx(pos[2 * i], pos[2 * i + 1], ax_, ay_);
          for (std::size_t w = rt.idx; w < rt.wp.size(); ++w) {
            int bx_, by_;
            toPx(rt.wp[w][0], rt.wp[w][1], bx_, by_);
            navdemo::plot_line(fp, W, H, ax_, ay_, bx_, by_, lr, lg, lb);
            ax_ = bx_;
            ay_ = by_;
          }
          int px, py;
          toPx(pos[2 * i], pos[2 * i + 1], px, py);
          navdemo::plot_disc(fp, W, H, px, py, 4, 12, 12, 12);
          navdemo::plot_disc(fp, W, H, px, py, 3, static_cast<unsigned char>(color[3 * i] * 255),
                             static_cast<unsigned char>(color[3 * i + 1] * 255),
                             static_cast<unsigned char>(color[3 * i + 2] * 255));
        }
      }
      char path[1024];
      std::snprintf(path, sizeof path, "%s/frame_%05ld.png", out.c_str(), frame);
      cvc::write_image(frameImg, path);
    } else {
      view.processUIEvents();
      view.render();
    }

    ++frame;
    if (frames > 0 && frame >= frames)
      break;
    if (frames == 0 && !capturing && view.windowClosed())
      break;
    if (frames == 0 && capturing && frame >= total)
      break;
  }

  cam.detach();
  if (!png.empty())
    view.writePNG(png);
  std::printf("nav_austin_live: done (%ld frames)\n", frame);
  return 0;
}
