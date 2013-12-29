#ifndef SDF_H
#define SDF_H

#include <multi_sdf/mds.h>
#include <multi_sdf/matrix.h>
#include <multi_sdf/kdtree.h>
#include <multi_sdf/util.h>

namespace multi_sdf
{

float
sdf( const Point& q, const Mesh &mesh, const vector<double>& weights,
     KdTree& kd_tree, const Triangulation& triang );

}

#endif // SDF_H

