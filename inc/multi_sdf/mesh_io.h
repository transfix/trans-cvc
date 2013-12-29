#ifndef MESH_IO_H
#define MESH_IO_H

#include <multi_sdf/mds.h>

namespace multi_sdf
{

enum FILE_TYPE
{
   OFF,
   COFF,
   RAW,
   RAWN,
   RAWC,
   RAWNC,
   STL,
   SMF
};

void
read_labeled_mesh(Mesh &mesh, const string& ip_filename, 
                  const FILE_TYPE& ftype, 
                  const bool& read_color_opacity, 
                  const bool& is_uniform);

void
read_labeled_mesh(Mesh &mesh,
		  int nverts, float* verts, float* colors, int ntris, int* tris);

void
write_mesh(const Mesh& mesh, const char* ofname, FILE_TYPE ftype, 
           bool write_color_opacity, bool use_input_mesh_color, 
           float r, float g, float b, float a);

}

#endif
