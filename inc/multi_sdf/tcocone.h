#ifndef TCOCONE_H
#define TCOCONE_H

#include <multi_sdf/datastruct.h>
#include <multi_sdf/util.h>
#include <multi_sdf/robust_cc.h>

namespace multi_sdf
{

void 
compute_poles( Triangulation &triang);

void 
mark_flat_vertices( Triangulation &triang,
		    double ratio, double cocone_phi, double flat_phi);

void
tcocone(const double DEFAULT_ANGLE,
        const double DEFAULT_SHARP,
	const double DEFAULT_FLAT,
	const double DEFAULT_RATIO,
	Triangulation &triang);

}

#endif // TCOCONE_H

