#ifndef RCOCONE_H
#define RCOCONE_H

#include <multi_sdf/datastruct.h>
#include <multi_sdf/util.h>
#include <multi_sdf/robust_cc.h>

namespace multi_sdf
{

void
robust_cocone(const double bb_ratio,
	      const double theta_ff,
	      const double theta_if,
	      Triangulation &triang,
	      const char *outfile_prefix);

}

#endif // RCOCONE_H

