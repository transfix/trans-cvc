/* $Id: dt.cpp 1527 2010-03-12 22:10:16Z transfix $ */

#include <multi_sdf/dt.h>

namespace multi_sdf
{

// -----------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------

// -- flatness marking ---
const double DEFAULT_ANGLE = M_PI / 8.0;      // Half of the co-cone angle.
const double DEFAULT_SHARP = 2 * M_PI / 3.0;  // Angle of sharp edges.
const double DEFAULT_RATIO = 1.2 * 1.2;       // Squared thinness factor.
const double DEFAULT_FLAT  = M_PI / 3.0;      // Angle for flatness Test

// -- robust cocone ---
const double DEFAULT_BIGBALL_RATIO  = (1./4.)*(1./4.);  // parameter to choose big balls.
const double DEFAULT_THETA_IF_d  = 5.0;       // parameter for infinite-finite deep intersection.
const double DEFAULT_THETA_FF_d  = 10.0;      // parameter for finite-finite deep intersection.

void
recon( const vector<Point>& pts, Triangulation& triang)
{
 // for flatness marking (in cocone)
 double flatness_ratio = DEFAULT_RATIO;
 double cocone_phi = DEFAULT_ANGLE;
 double flat_phi = DEFAULT_FLAT;

 for(int i = 0; i < pts.size(); i ++)
    triang.insert(pts[i]);

 // --- Init ----
 cerr << "Init ";
 initialize(triang);
 cerr << ".";
 // compute voronoi vertex
 compute_voronoi_vertex_and_cell_radius(triang);
 cerr << ". done." << endl;

 // ---- Reconstruction -----
 cerr << "TC ";
 tcocone(cocone_phi, DEFAULT_SHARP, flat_phi, flatness_ratio, triang);
 cerr << " done." << endl;
 write_wt(triang, "test");

 return;
}

}
