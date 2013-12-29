/* $Id: mesh_io.cpp 1527 2010-03-12 22:10:16Z transfix $ */

#include <multi_sdf/mesh_io.h>

namespace multi_sdf
{

// ------------------------
// read_mesh
// ------------------------
void
read_labeled_mesh(Mesh &mesh, const string& ip_filename, 
                  const FILE_TYPE& ftype, 
                  const bool& read_color_opacity, 
                  const bool& is_uniform)
{
  FILE* fp = fopen(ip_filename.c_str(), "r");
  char header[10];
  switch(ftype)
  {
     case OFF: 
        fscanf(fp, "%s", header);
        if( strcmp(header, "OFF") != 0 )
        {
           cerr << "Error : Missing header OFF" << endl;
           exit(1);
        }
        break;
     case RAW: 
     case RAWC: 
     case RAWNC: 
        break;
     default:
        break;
  };

  int nv = -1, nf = -1;
  fscanf(fp, "%d %d", &nv, &nf);
  mesh.set_nv(nv);
  mesh.set_nf(nf);
  if( ftype == OFF )
  {
     int skip_int = -1;
     fscanf(fp, "%d ", &skip_int);
  }
  for(int i = 0; i < nv; i ++)
  {
     float x,y,z;
     fscanf(fp, "%f %f %f", &x,&y,&z);
     float r=1,g=1,b=1;
     if(ftype == RAWC) fscanf(fp,"%f %f %f", &r,&g,&b);
     float nx,ny,nz;
     if(ftype == RAWNC) fscanf(fp,"%f %f %f %f %f %f", &nx,&ny,&nz,&r,&g,&b);

     mesh.vert_list.push_back(MVertex(Point(x,y,z)));
     mesh.vert_list[mesh.vert_list.size() - 1].id = i;
     mesh.vert_list[mesh.vert_list.size() - 1].set_iso(true);
     if( ftype == RAWC || ftype == RAWNC )
     {
        mesh.vert_list[mesh.vert_list.size() - 1].vert_color[0] = r;
        mesh.vert_list[mesh.vert_list.size() - 1].vert_color[1] = g;
        mesh.vert_list[mesh.vert_list.size() - 1].vert_color[2] = b;
     }
  }

  int ne = 0;
  for(int i = 0; i < nf; i ++)
  {
     int np = -1, v1 = -1, v2 = -1, v3 = -1, label = -1;
     float r = 1, g = 1, b = 1, a = 1;
     if( ftype == OFF && read_color_opacity && ! is_uniform)
        fscanf(fp,"%d %d %d %d %f %f %f %f %d", &np, &v1, &v2, &v3, &r, &g, &b, &a, &label);
     else if( ftype == OFF && read_color_opacity )
        fscanf(fp,"%d %d %d %d %f %f %f %f", &np, &v1, &v2, &v3, &r, &g, &b, &a);
     else if( ftype == OFF && ! is_uniform )
        fscanf(fp,"%d %d %d %d %d", &np, &v1, &v2, &v3, &label);
     else if( ftype == OFF )
        fscanf(fp,"%d %d %d %d", &np, &v1, &v2, &v3);
     else if( ! is_uniform )
        fscanf(fp,"%d %d %d %d", &v1, &v2, &v3, &label);
     else
        fscanf(fp,"%d %d %d", &v1, &v2, &v3);

     if( ftype == OFF && np != 3 ) 
     {
        cerr << endl << "Please check the file format. " << endl;
        cerr << "This program accepts OFF format with triangles plus weights." << endl;
        exit(1);
     }

     mesh.face_list.push_back(MFace(v1,v2,v3));
     mesh.face_list[mesh.face_list.size() - 1].id = i;
     mesh.face_list[mesh.face_list.size() - 1].label = label;
     // mesh.face_list[mesh.face_list.size() - 1].label = a;
     // mesh.face_list[mesh.face_list.size() - 1].w = a;
     // mesh.face_list[mesh.face_list.size() - 1].w = 1;

     // add the incident faces in all three vertices.
     mesh.vert_list[v1].add_inc_face(i);
     mesh.vert_list[v1].set_iso(false);
     mesh.vert_list[v2].add_inc_face(i);
     mesh.vert_list[v2].set_iso(false);
     mesh.vert_list[v3].add_inc_face(i);
     mesh.vert_list[v3].set_iso(false);


     if( read_color_opacity )
     {
        mesh.face_list[mesh.face_list.size() - 1].set_color(0, r);
        mesh.face_list[mesh.face_list.size() - 1].set_color(1, g);
        mesh.face_list[mesh.face_list.size() - 1].set_color(2, b);
        mesh.face_list[mesh.face_list.size() - 1].set_color(3, a);
     }
     // book keeping for incident vertices for a vertex.
     if(!mesh.vert_list[v1].is_inc_vert(v2))
     {
        mesh.vert_list[v1].add_inc_vert(v2);
	assert(!mesh.vert_list[v2].is_inc_vert(v1));
	mesh.vert_list[v2].add_inc_vert(v1);
	// here create an edge between v1 and v2.
	// add the face to the edge
	MEdge me = MEdge(v1,v2);
	me.add_inc_face(i);
	mesh.edge_list.push_back(me);
	ne ++;
	// insert the id of the incident edge in the vertices
	mesh.vert_list[v1].add_inc_edge(mesh.edge_list.size() - 1);
	mesh.vert_list[v2].add_inc_edge(mesh.edge_list.size() - 1);
     }
     else
     {
        // the edge is already there.
	// find the edge id using the vertex indices.
	MVertex mv = mesh.vert_list[v1];
	int eid = -1;
	assert(mv.get_eid(v2, eid));
	assert(eid != -1);
	// add the face to the edge.
	if(mesh.edge_list[eid].num_inc_face <= 1)
	{
	   assert(mesh.edge_list[eid].num_inc_face == 1);
	   mesh.edge_list[eid].add_inc_face(i);
	}
	else
	   mesh.edge_list[eid].num_inc_face ++; // nonmanifold.
     }

     if(!mesh.vert_list[v2].is_inc_vert(v3))
     {
        mesh.vert_list[v2].add_inc_vert(v3);
	assert(!mesh.vert_list[v3].is_inc_vert(v2));
	mesh.vert_list[v3].add_inc_vert(v2);
	// here create an edge between v2 and v3.
	// add the face to the edge
	MEdge me = MEdge(v2,v3);
	me.add_inc_face(i);
	mesh.edge_list.push_back(me);
	ne ++;

	// insert the id of the incident edge in the vertices
	mesh.vert_list[v2].add_inc_edge(mesh.edge_list.size() - 1);
	mesh.vert_list[v3].add_inc_edge(mesh.edge_list.size() - 1);
     }
     else
     {
	// the edge is already there.
	// find the edge id using the vertex indices.
	MVertex mv = mesh.vert_list[v2];
	int eid = -1;
	assert(mv.get_eid(v3, eid));
	assert(eid != -1);
	// add the face to the edge.
	if(mesh.edge_list[eid].num_inc_face <= 1)
	{
	   assert(mesh.edge_list[eid].num_inc_face == 1);
	   mesh.edge_list[eid].add_inc_face(i);
	}
	else
	   mesh.edge_list[eid].num_inc_face ++; // nonmanifold.
     }

     if(!mesh.vert_list[v3].is_inc_vert(v1))
     {
        mesh.vert_list[v3].add_inc_vert(v1);
	assert(!mesh.vert_list[v1].is_inc_vert(v3));
	mesh.vert_list[v1].add_inc_vert(v3);
	// here create an edge between v3 and v1.
	// add the face to the edge
	MEdge me = MEdge(v3,v1);
	me.add_inc_face(i);
	mesh.edge_list.push_back(me);
	ne ++;
	// insert the id of the incident edge in the vertices
	mesh.vert_list[v3].add_inc_edge(mesh.edge_list.size() - 1);
	mesh.vert_list[v1].add_inc_edge(mesh.edge_list.size() - 1);
     }
     else
     {
        // the edge is already there.
	// find the edge id using the vertex indices.
	MVertex mv = mesh.vert_list[v3];
	int eid = -1;
	assert(mv.get_eid(v1, eid));
	assert(eid != -1);
	// add the face to the edge.
	if(mesh.edge_list[eid].num_inc_face <= 1)
	{
	   assert(mesh.edge_list[eid].num_inc_face == 1);
	   mesh.edge_list[eid].add_inc_face(i);
	}
	else
	   mesh.edge_list[eid].num_inc_face ++; // nonmanifold.
     }
  }
  fclose(fp);

  mesh.set_ne(ne);
  // add the edges in the face so that
  // if v1 is the ith corner of a face 
  // then v2<->v3 is the ith edge.
  for(int i = 0; i < nf; i ++)
  {
     for(int j = 0; j < 3; j ++)
     {
        int u = mesh.face_list[i].get_corner((j+1)%3);
	int w = mesh.face_list[i].get_corner((j+2)%3);
	// find the edge id connecting u amd w.
	int eid = -1;
	assert(mesh.vert_list[u].get_eid(w, eid));
	assert(eid != -1);
	// this edge should be the jth edge of the face.
	mesh.face_list[i].set_edge(j, eid);
     }
  }
}

// ------------------------
// read_mesh
// ------------------------
void
read_labeled_mesh(Mesh &mesh,
		  int nverts, float* verts, float* colors, int ntris, int* tris)
{
  mesh.set_nv(nverts);
  mesh.set_nf(ntris);
  for(int i = 0; i < nverts; i++)
  {
    float x,y,z;
    x = verts[i*3+0];
    y = verts[i*3+1];
    z = verts[i*3+2];
    float r=1.0,g=1.0,b=1.0;
    if(colors)
      {
	r = colors[i*3+0];
	g = colors[i*3+1];
	b = colors[i*3+2];
      }

     mesh.vert_list.push_back(MVertex(Point(x,y,z)));
     mesh.vert_list[mesh.vert_list.size() - 1].id = i;
     mesh.vert_list[mesh.vert_list.size() - 1].set_iso(true);
     mesh.vert_list[mesh.vert_list.size() - 1].vert_color[0] = r;
     mesh.vert_list[mesh.vert_list.size() - 1].vert_color[1] = g;
     mesh.vert_list[mesh.vert_list.size() - 1].vert_color[2] = b;
  }

  int ne = 0;
  for(int i = 0; i < ntris; i ++)
  {
     int v1 = -1, v2 = -1, v3 = -1, label = -1;
     
     v1 = tris[i*3+0];
     v2 = tris[i*3+1];
     v3 = tris[i*3+2];

     mesh.face_list.push_back(MFace(v1,v2,v3));
     mesh.face_list[mesh.face_list.size() - 1].id = i;
     mesh.face_list[mesh.face_list.size() - 1].label = label;
     
     // add the incident faces in all three vertices.
     mesh.vert_list[v1].add_inc_face(i);
     mesh.vert_list[v1].set_iso(false);
     mesh.vert_list[v2].add_inc_face(i);
     mesh.vert_list[v2].set_iso(false);
     mesh.vert_list[v3].add_inc_face(i);
     mesh.vert_list[v3].set_iso(false);

     // book keeping for incident vertices for a vertex.
     if(!mesh.vert_list[v1].is_inc_vert(v2))
     {
        mesh.vert_list[v1].add_inc_vert(v2);
	assert(!mesh.vert_list[v2].is_inc_vert(v1));
	mesh.vert_list[v2].add_inc_vert(v1);
	// here create an edge between v1 and v2.
	// add the face to the edge
	MEdge me = MEdge(v1,v2);
	me.add_inc_face(i);
	mesh.edge_list.push_back(me);
	ne ++;
	// insert the id of the incident edge in the vertices
	mesh.vert_list[v1].add_inc_edge(mesh.edge_list.size() - 1);
	mesh.vert_list[v2].add_inc_edge(mesh.edge_list.size() - 1);
     }
     else
     {
        // the edge is already there.
	// find the edge id using the vertex indices.
	MVertex mv = mesh.vert_list[v1];
	int eid = -1;
	assert(mv.get_eid(v2, eid));
	assert(eid != -1);
	// add the face to the edge.
	if(mesh.edge_list[eid].num_inc_face <= 1)
	{
	   assert(mesh.edge_list[eid].num_inc_face == 1);
	   mesh.edge_list[eid].add_inc_face(i);
	}
	else
	   mesh.edge_list[eid].num_inc_face ++; // nonmanifold.
     }

     if(!mesh.vert_list[v2].is_inc_vert(v3))
     {
        mesh.vert_list[v2].add_inc_vert(v3);
	assert(!mesh.vert_list[v3].is_inc_vert(v2));
	mesh.vert_list[v3].add_inc_vert(v2);
	// here create an edge between v2 and v3.
	// add the face to the edge
	MEdge me = MEdge(v2,v3);
	me.add_inc_face(i);
	mesh.edge_list.push_back(me);
	ne ++;

	// insert the id of the incident edge in the vertices
	mesh.vert_list[v2].add_inc_edge(mesh.edge_list.size() - 1);
	mesh.vert_list[v3].add_inc_edge(mesh.edge_list.size() - 1);
     }
     else
     {
	// the edge is already there.
	// find the edge id using the vertex indices.
	MVertex mv = mesh.vert_list[v2];
	int eid = -1;
	assert(mv.get_eid(v3, eid));
	assert(eid != -1);
	// add the face to the edge.
	if(mesh.edge_list[eid].num_inc_face <= 1)
	{
	   assert(mesh.edge_list[eid].num_inc_face == 1);
	   mesh.edge_list[eid].add_inc_face(i);
	}
	else
	   mesh.edge_list[eid].num_inc_face ++; // nonmanifold.
     }

     if(!mesh.vert_list[v3].is_inc_vert(v1))
     {
        mesh.vert_list[v3].add_inc_vert(v1);
	assert(!mesh.vert_list[v1].is_inc_vert(v3));
	mesh.vert_list[v1].add_inc_vert(v3);
	// here create an edge between v3 and v1.
	// add the face to the edge
	MEdge me = MEdge(v3,v1);
	me.add_inc_face(i);
	mesh.edge_list.push_back(me);
	ne ++;
	// insert the id of the incident edge in the vertices
	mesh.vert_list[v3].add_inc_edge(mesh.edge_list.size() - 1);
	mesh.vert_list[v1].add_inc_edge(mesh.edge_list.size() - 1);
     }
     else
     {
        // the edge is already there.
	// find the edge id using the vertex indices.
	MVertex mv = mesh.vert_list[v3];
	int eid = -1;
	assert(mv.get_eid(v1, eid));
	assert(eid != -1);
	// add the face to the edge.
	if(mesh.edge_list[eid].num_inc_face <= 1)
	{
	   assert(mesh.edge_list[eid].num_inc_face == 1);
	   mesh.edge_list[eid].add_inc_face(i);
	}
	else
	   mesh.edge_list[eid].num_inc_face ++; // nonmanifold.
     }
  }

  mesh.set_ne(ne);
  // add the edges in the face so that
  // if v1 is the ith corner of a face 
  // then v2<->v3 is the ith edge.
  for(int i = 0; i < ntris; i ++)
  {
     for(int j = 0; j < 3; j ++)
     {
        int u = mesh.face_list[i].get_corner((j+1)%3);
	int w = mesh.face_list[i].get_corner((j+2)%3);
	// find the edge id connecting u amd w.
	int eid = -1;
	assert(mesh.vert_list[u].get_eid(w, eid));
	assert(eid != -1);
	// this edge should be the jth edge of the face.
	mesh.face_list[i].set_edge(j, eid);
     }
  }
}

void
write_mesh(const Mesh& mesh, const char* ofname, FILE_TYPE ftype, 
           bool write_color_opacity, bool use_input_mesh_color, 
           float r, float g, float b, float a)
{
  ofstream fout;
  fout.open(ofname);
  if( ! fout )
  {
     cerr << "Can not open " << ofname << " for writing." << endl;
     exit(0);
  }

  switch(ftype)
  {
     case OFF: 
        fout << "OFF" << endl;
        fout << mesh.get_nv() << " " << mesh.get_nf() << " 0" << endl;
        break;
     case COFF: 
        fout << "COFF" << endl;
        fout << mesh.get_nv() << " " << mesh.get_nf() << " 0" << endl;
        break;
     case RAW: 
     case RAWN: 
     case RAWC: 
     case RAWNC: 
        fout << mesh.get_nv() << " " << mesh.get_nf() << endl;
        break;
     case STL: 
        break;
     default:
        break;
  };

  for(int i = 0; i < mesh.get_nv(); i ++)
  {
     fout << mesh.vert_list[i].point() << " ";
     switch(ftype)
     {
        case OFF: 
           fout << endl;
           break;
        case COFF:
           if( use_input_mesh_color )
              fout << mesh.vert_list[i].vert_color[0] << " "
                   << mesh.vert_list[i].vert_color[1] << " "
                   << mesh.vert_list[i].vert_color[2] << " "  
                   << mesh.vert_list[i].vert_color[3] << endl;
           else
              fout << r << " " << g << " " << b << " " << a << endl;
           break;
        case RAW: 
           fout << endl;
           break;
        case RAWN: 
           fout << mesh.vert_list[i].mesh_normal() << endl;
           break;
        case RAWC: 
           if( use_input_mesh_color )
              fout << endl;
           else
              fout << r << " " << g << " " << b << endl;
           break;
        case RAWNC: 
           fout << endl;
           break;
        default:
           break;
     };
  }

  for(int i = 0; i < mesh.get_nf(); i ++)
  {
     if( ftype == OFF || ftype == COFF )
        fout << "3\t";
     fout << mesh.face_list[i].get_corner(0) << " "
          << mesh.face_list[i].get_corner(1) << " "
          << mesh.face_list[i].get_corner(2) << " ";
     if( ftype == OFF )
     {
        if( write_color_opacity )
        {
           if( use_input_mesh_color )
              fout << mesh.face_list[i].get_color(0) << " " 
                   << mesh.face_list[i].get_color(1) << " " 
                   << mesh.face_list[i].get_color(2) << " " 
                   << mesh.face_list[i].get_color(3);
           else
              fout << r << " " << g << " " << b << " " << a;
        }
     }
     fout << endl;
  }
  fout.close();
}

}
