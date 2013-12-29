#include <stdexcept>
#include <multi_sdf/sdf.h>

namespace multi_sdf
{

float
sdf( const Point& q, const Mesh &mesh, const vector<double>& weights,
     KdTree& kd_tree, const Triangulation& triang )
{
  VECTOR3 query (CGAL::to_double(q.x()),
                 CGAL::to_double(q.y()),
                 CGAL::to_double(q.z()));
  kd_tree.queryPosition(query);
  // Initialize the search structure, and search all N points
  int n_vid = kd_tree.getNeighbourPositionIndex(0);
  if(n_vid == -1) throw std::runtime_error("No nearest neighbor. MDS empty?");
  CGAL_assertion( ! mesh.vert_list[n_vid].iso());
  MVertex nv = mesh.vert_list[n_vid];
  double min_sq_d = HUGE;
  int n_fid = -1;
  for(int i = 0; i < nv.num_inc_face; i ++)
  {
     MFace f = mesh.face_list[nv.inc_face(i)];
     Point p[3] = {mesh.vert_list[f.get_corner(0)].point(),
                   mesh.vert_list[f.get_corner(1)].point(),
                   mesh.vert_list[f.get_corner(2)].point()};
     Triangle_3 t (p[0], p[1], p[2]);
     Plane_3 H (p[0], p[1], p[2]);
     Point _q = H.projection(q);
     // check if _q is inside t.
     if( t.has_on(_q) )
     {
        double sq_d = CGAL::to_double((q-_q)*(q-_q));
        if( sq_d < min_sq_d ) { min_sq_d = sq_d; n_fid = nv.inc_face(i); }
     }
     else
     {
        for(int j = 0; j < 3; j ++)
        {
           double _d = CGAL::to_double(CGAL::squared_distance(_q,Segment(p[j], p[(j+1)%3])));
           double sq_d = CGAL::to_double((q-_q)*(q-_q)) + _d;
           if( sq_d < min_sq_d ) { min_sq_d = sq_d; n_fid = nv.inc_face(i); }
        }
     }
  }

  // locate the query point in the triang which is already tagged
  // with in-out flag by the reconstruction.
  bool is_q_outside = false;
  Triangulation::Locate_type lt;
  int u = -1, v = -1;
  Cell_handle c = triang.locate(q, lt, u, v);
  if( lt == Triangulation::OUTSIDE_AFFINE_HULL )
  {
     is_q_outside = true;
     cerr << "Point " << q << " is outside the affine hull." << endl;
  }
  else if( lt == Triangulation::OUTSIDE_CONVEX_HULL )
     is_q_outside = true;
  else
  {
     if( lt == Triangulation::CELL )
     {
        if( c->outside )
           is_q_outside = true;
        else
           is_q_outside = false;
     }
     else if( lt == Triangulation::FACET )
     {
        Cell_handle _c = c->neighbor(u);
        if( c->outside && _c->outside )
           is_q_outside = true;
        else
           is_q_outside = false;
     }
     else if( lt == Triangulation::EDGE )
     {
        if( is_outside_VF(triang, Edge(c, u, v)) )
           is_q_outside = true;
        else
           is_q_outside = false;
     }
     else
     {
        CGAL_assertion( lt == Triangulation::VERTEX );
        is_q_outside = false;
     }
  }

  double w;
  if(weights.size() && mesh.face_list[n_fid].label != -1) w = weights[mesh.face_list[n_fid].label];
  else w = 1.0;
  // double w = mesh.face_list[n_fid].w;

  double gen_sdf = 0;
  if( is_q_outside ) gen_sdf = w*sqrt(min_sq_d);
  else gen_sdf = -w*min_sq_d;

#if 0
  double MAX = 10, MIN = -10;

  if( gen_sdf > MAX ) gen_sdf = MAX;
  if( gen_sdf < MIN ) gen_sdf = MIN;
#endif
  return gen_sdf;
}

}
