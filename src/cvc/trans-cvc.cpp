#include <cvc/app.h>
#include <cvc/volume_file_info.h>
#include <cvc/geometry.h>
#include <cvc/algorithm.h>

#include <boost/regex.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/current_function.hpp>

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cstdlib>

namespace
{
  /*
   * Utils
   */

  bool is_geometry(const boost::any& data)
  {
    return data.type() == typeid(CVC_NAMESPACE::geometry);
  }

  bool is_volume(const boost::any& data)
  {
    return data.type() == typeid(CVC_NAMESPACE::volume);
  }

  bool is_volume_file_info(const boost::any& data)
  {
    return data.type() == typeid(CVC_NAMESPACE::volume_file_info);
  }

  bool is_geometry_filename(const std::string& filename)
  {
    std::string errors;
    boost::regex file_extension("^(.*)(\\.\\S*)$");
    boost::smatch what;

    if(boost::regex_match(filename, what, file_extension))
      if(what[2].compare(".raw") == 0 ||
         what[2].compare(".rawn") == 0 ||
         what[2].compare(".rawnc") == 0 ||
         what[2].compare(".rawc") == 0 ||
         what[2].compare(".off") == 0)
        return true;
    return false;
  }

  bool is_volume_filename(const std::string& filename)
  {
    std::string errors;
    boost::regex file_extension("^(.*)(\\.\\S*)$");
    boost::smatch what;

    if(boost::regex_match(filename, what, file_extension))
      {
	std::vector<std::string> exts = 
	  CVC_NAMESPACE::volume_file_io::getExtensions();
	BOOST_FOREACH(std::string& ext, exts)
	  if(what[2] == ext)
	    return true;
      }
    return false;
  }

  boost::any load(const std::string& filename)
  {
    std::string errors;
    boost::regex file_extension("^(.*)(\\.\\S*)$");
    boost::smatch what;
    
    if(is_geometry_filename(filename))
      {
        CVC_NAMESPACE::geometry geo(filename);
        return geo;
      }
    else if(is_volume_filename(filename))
      {
        CVC_NAMESPACE::volume vol(filename);
        return vol;
      }
    else
      throw CVC_NAMESPACE::read_error(BOOST_CURRENT_FUNCTION);

    return boost::any();
  }

  void save(const boost::any& data, const std::string& filename)
  {
    if(is_geometry(data))
      {
        CVC_NAMESPACE::geometry geo = boost::any_cast<CVC_NAMESPACE::geometry>(data);
        geo.write(filename);
      }
    else if(is_volume(data))
      {
        CVC_NAMESPACE::volume vol = boost::any_cast<CVC_NAMESPACE::volume>(data);
        vol.write(filename);
      }
    else
      throw CVC_NAMESPACE::write_error(BOOST_CURRENT_FUNCTION);
  }

  /*
   * Commands
   */
  
  typedef boost::function<void (const std::vector<std::string>&)> command_func;
  typedef boost::tuple<command_func, std::string>                 command;
  typedef std::map<std::string, command>                          command_map;
  command_map commands;

  void help(const std::vector<std::string>& args)
  {
    using namespace std;
    cout << "Usage: trans-cvc <command> <command args>" << endl << endl;
    for(command_map::iterator i = commands.begin();
        i != commands.end();
        ++i)
      {
        cout << " - " << i->first << endl;
        cout << i->second.get<1>() << endl << endl;
      }
  }

  void copy(const std::vector<std::string>& args)
  {
    using namespace std;
    using namespace boost;
    using namespace CVC_NAMESPACE;
    thread_info ti(BOOST_CURRENT_FUNCTION);

    if(args.empty()) throw command_line_error("Missing input filename.");
    if(args.size() < 2)  throw command_line_error("Missing output volume filename");
    save(load(args[0]),args[1]);
  }

  void info(const std::vector<std::string>& args)
  {
    using namespace std;
    using namespace boost;
    using namespace CVC_NAMESPACE;
    thread_info ti(BOOST_CURRENT_FUNCTION);

    if(args.empty()) throw command_line_error("Missing input filename.");

    if(is_geometry_filename(args[0]))
      {
        geometry geo(args[0]);
        if(geo.num_points() > 0) cout << "Num vertices: " << geo.num_points() << endl;
        if(geo.num_lines() > 0) cout << "Num lines: " << geo.num_lines() << endl;
        if(geo.num_tris() > 0) cout << "Num tris: " << geo.num_tris() << endl;
        if(geo.num_quads() > 0) cout << "Num quads: " << geo.num_quads() << endl;

        point_t min_pt = geo.min_point();
        point_t max_pt = geo.max_point();
        cout << "min_point: " << 
          str(format("%1%,%2%,%3%") % min_pt[0] % min_pt[1] % min_pt[2]) << endl;
        cout << "max_point: " << 
          str(format("%1%,%2%,%3%") % max_pt[0] % max_pt[1] % max_pt[2]) << endl;
      }
    else if(is_volume_filename(args[0]))
      {
        volume_file_info volinfo;
        volinfo.read(args[0]);
        cout << volinfo.filename() << ":" <<endl;
        cout << "Num Variables: " << volinfo.numVariables() << endl;
        cout << "Num Timesteps: " << volinfo.numTimesteps() << endl;
        cout << "Dimension: " << volinfo.XDim() << "x" << volinfo.YDim() << "x" << volinfo.ZDim() << endl;
        cout << "Bounding box: ";
        cout << "(" << volinfo.boundingBox().minx << "," << volinfo.boundingBox().miny << "," << volinfo.boundingBox().minz << ") ";
        cout << "(" << volinfo.boundingBox().maxx << "," << volinfo.boundingBox().maxy << "," << volinfo.boundingBox().maxz << ") ";
        cout << endl;
        cout << "Span: " << "(" << volinfo.XSpan() << "," << volinfo.YSpan() << "," << volinfo.ZSpan() << ") " << endl;
        double volmin = volinfo.min(), volmax = volinfo.max();
        for(unsigned int i = 0; i<volinfo.numVariables(); i++)
          {
            cout << "Name of var " << i << ": " << volinfo.name(i) << endl;
            cout << "Voxel type of var " << i << ": " << volinfo.voxelTypeStr(i) << endl;
            for(unsigned int j = 0; j<volinfo.numTimesteps(); j++)
              {
                if(volmin > volinfo.min(i,j)) volmin = volinfo.min(i,j);
                if(volmax < volinfo.max(i,j)) volmax = volinfo.max(i,j);
                cout << "Min voxel value of var " << i << ", timestep " << j << ": " << volinfo.min(i,j) << endl;
                cout << "Max voxel value of var " << i << ", timestep " << j << ": " << volinfo.max(i,j) << endl;
              }
          }
        cout << "Min voxel value (of whole dataset): " << volmin << endl;
        cout << "Max voxel value (of whole dataset): " << volmax << endl;       
      }
    else
      throw read_error("Unknown data");
  }

  // 01/11/2014 - Joe R. - moved bounding box to last argument and made optional.
  void sdf(const std::vector<std::string>& args)
  {
    using namespace std;
    using namespace CVC_NAMESPACE;
    thread_info ti(BOOST_CURRENT_FUNCTION);
    if(args.empty())     throw command_line_error("Missing geometry filename");
    if(args.size() < 2)  throw command_line_error("Missing dimension");
    if(args.size() < 3)  throw command_line_error("Missing output volume filename");

    geometry geom(args[0]);
    bounding_box bbox = geom.extents();
    if(args.size() >= 4) bbox = bounding_box(args[3]);

    sdf(geom,
        dimension(args[1]),
        bbox)
      .write(args[2]);
  }

  void iso(const std::vector<std::string>& args)
  {
    using namespace std;
    using namespace boost;
    using namespace CVC_NAMESPACE;
    thread_info ti(BOOST_CURRENT_FUNCTION);
    if(args.empty())     throw command_line_error("Missing volume filename");
    if(args.size() < 2)  throw command_line_error("Missing isovalue");
    if(args.size() < 3)  throw command_line_error("Missing output geometry filename");
    iso(volume(args[0]),
        lexical_cast<double>(args[1]))
      .write(args[2]);
  }

  class init_commands
  {
  public:
    init_commands()
    {
      using namespace std;
      using namespace boost;
      
      commands["copy"] =
        make_tuple(command_func(copy),
                   string("copy <input filename> <output filename>\n"
                          "Copies input to output file, possibly converting file formats."));
      commands["help"] = 
        make_tuple(command_func(help),
		   string("Prints command list."));
      commands["iso"] = 
        make_tuple(command_func(iso),
                   string("iso <volume filename> <isovalue> <output geometry filename>\n"
                          "Computes isosurface geometry for a given volume and isovalue.\n"
                          "<isovalue> - value describing the surface to be computed."));                          
      commands["sdf"] = make_tuple(command_func(sdf),
                                   string("sdf <geometry filename> <dimension> <output volume filename> [bounding_box]\n"
                                          "Computes a signed distance function volume of the geometry specified."
                                          "<dimension> - comma separated list of 3 integers specifying output volume dimension\n"
                                          "[bounding_box] - comma separated list of 6 floats specifying output volume bounding box"
					  "                 If unspecified, defaults to the extents of the input geometry."));
      commands["info"] = make_tuple(command_func(info),
                                    string("info <filename>\n"
                                           "Prints info about the specified file."));
    }
  } static_init;
}

int main(int argc, char **argv)
{
  using namespace std;
  using namespace CVC_NAMESPACE;
  
  vector<string> args;
  for(int i = 0; i < argc; i++)
    args.push_back(argv[i]);

  try
    {
      args.erase(args.begin()); //no need for the first arg
      if(args.empty()) throw command_line_error("Missing command string");
      string cmd = args[0];
      args.erase(args.begin()); //erase the command string

      if(commands.find(cmd)==commands.end())
        throw command_line_error("Invalid command.");

      commands[cmd].get<0>()(args);
    }
  catch(command_line_error& e)
    {
      if(!e.what_str().empty()) cout << "Error: " << e.what_str() << endl;
      cout << "Usage: " << argv[0] << " <command> <command args>" << endl;
      return EXIT_FAILURE;
    }
  catch(CVC_NAMESPACE::exception& e)
    {
      cerr << "Exception: " << e.what() << endl;
    }

  return EXIT_SUCCESS;
}
