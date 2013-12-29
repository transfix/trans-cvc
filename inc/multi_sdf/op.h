#ifndef OP_H
#define OP_H

#include <multi_sdf/datastruct.h>
#include <multi_sdf/util.h>

namespace multi_sdf
{

void
draw_ray(const Ray& ray, 
         const double& r, 
         const double& g, 
         const double& b, 
         const double& a,
         ofstream& fout);

void
draw_segment(const Segment& segment,
             const double& r, 
             const double& g, 
             const double& b, 
             const double& a,
             ofstream& fout);

void
draw_poly(const vector<Point>& poly,
          const double& r, 
          const double& g, 
          const double& b, 
          const double& a,
          ofstream& fout);

void
draw_VF(const Triangulation& triang,
        const Edge& dual_e, 
        const double& r, 
        const double& g, 
        const double& b, 
        const double& a,
        ofstream& fout);

void
draw_tetra(const Cell_handle& cell,
           const double& r, 
           const double& g, 
           const double& b, 
           const double& a,
           ofstream& fout);

void 
write_wt( const Triangulation &triang,
	  const char* file_prefix);

void 
write_iobdy( const Triangulation &triang,
	     const char* file_prefix);

}

#endif // OP_H

