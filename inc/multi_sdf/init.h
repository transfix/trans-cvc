#ifndef INIT_H
#define INIT_H

#include <multi_sdf/datastruct.h>

namespace multi_sdf
{

void
initialize(Triangulation &triang);
void
compute_voronoi_vertex_and_cell_radius(Triangulation &triang);

}

#endif // INIT_H

