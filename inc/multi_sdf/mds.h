#ifndef MDS_H
#define MDS_H

#include <multi_sdf/datastruct.h>

#warning REMOVE using namespace std declaration at some point!
using namespace std;

namespace multi_sdf
{

class MVertex
{
public:
	MVertex()
	{ init(); }
	MVertex(const Point& p) : coord(p)
	{ init(); }

	void set_point(const Point& p)			{ coord = p; }
	Point point() const				{ return coord; }
 
        void set_mesh_normal( const Vector& v)          { mv_normal = v; }
        Vector mesh_normal() const                      { return mv_normal; }

        void set_iso(const bool& b)                     { mv_iso = b; }
        bool iso() const                                { return mv_iso; }

	void add_inc_vert(const int i)			{ inc_vert_list.push_back(i);
							  num_inc_vert ++; }
	int inc_vert(int i) const			{ return inc_vert_list[i];}
	bool is_inc_vert(const int v)			
	{
	   for(int i = 0; i < num_inc_vert; i ++)
	      if(inc_vert_list[i] == v)
	         return true;
	   return false;
	}
	void add_inc_edge(const int i)			{ inc_edge_list.push_back(i);
                                                          num_inc_edge ++; }
	int inc_edge(int i) const			{ return inc_edge_list[i];}
	bool get_eid(const int v, int &eid)			
	{
	   eid = -1;
	   assert(num_inc_vert == num_inc_edge);
	   for(int i = 0; i < num_inc_vert; i ++)
	      if(inc_vert_list[i] == v)
	      {
	         eid = inc_edge_list[i];
		 return true;
	      }
	      return false;
	}

	void add_inc_face(const int i)			{ inc_face_list.push_back(i);
							  num_inc_face ++; }
	int inc_face(int i) const			{ return inc_face_list[i];}

	int id;
	bool visited;
	int num_inc_edge;
	int num_inc_vert;
	int num_inc_face;

        vector<double> vert_color;

	inline void init()
	{
	   id = -1;
           mv_normal = CGAL::NULL_VECTOR;
           mv_iso = true;
	   num_inc_vert = 0;
	   num_inc_edge = 0;
	   num_inc_face = 0;
	   inc_vert_list.clear();
	   inc_edge_list.clear();
           inc_face_list.clear();
           vert_color.resize(4,1);
	}

private:
	Point coord;
        Vector mv_normal;
        bool mv_iso;
	vector<int> inc_vert_list;
	vector<int> inc_edge_list;
	vector<int> inc_face_list;
};

class MEdge
{
public:
	MEdge()
	{ init(); }

	MEdge(const int v1, const int v2)
	{ init(); endpoint[0] = v1; endpoint[1] = v2; }

	void set_endpoint(const int i, const int val)	{ endpoint[i] = val; }
	int get_endpoint(int i) const			{ return endpoint[i]; }

	void add_inc_face(const int fid)		{ inc_face[num_inc_face] = fid;
       							  num_inc_face ++; }
	void get_inc_face(int &f1, int& f2)		{ f1 = inc_face[0]; f2 = inc_face[1]; }

	int num_inc_face;
	inline void init()
	{
	   endpoint[0] = endpoint[1] = -1;
	   inc_face[0] = inc_face[1] = -1;
	   num_inc_face = 0;
	}

private:
	int endpoint[2];
	int inc_face[2];
};


class MFace
{
public:
	MFace()
	{ init(); }
	MFace(int v1, int v2, int v3)
	{ init(); corner[0] = v1; corner[1] = v2; corner[2] = v3; }

	void set_corner(const int i, const int val)	{ corner[i] = val; }
	int get_corner(int i) const			{ return corner[i]; }

	void set_edge(const int i, const int val)	{ edge_array[i] = val; }
	int get_edge(int i) const			{ return edge_array[i]; }

	void set_color(const int i, const double c)	{ rgba[i] = c;}
	double get_color(int i) const			{ return rgba[i];}
	
        int id;
        int label;
        double w;

	inline void init()
	{
	   corner[0] = corner[1] = corner[2] = -1;
	   edge_array[0] = edge_array[1] = edge_array[2] = -1;
	   rgba[0] = rgba[1] = rgba[2] = rgba[3] = 1;
	   id = -1;
           label = -1;
           w = 1;
	}

private:
	int corner[3];
	int edge_array[3];
	double rgba[4];
};

class Mesh
{
public:
	Mesh()
	{ init(); }
	Mesh(int v, int f)
	{ init(); nv = v; nf = f; }

	void set_nv(const int n)  	{nv = n;}
	int get_nv() const		{return nv;}
	void set_nf(const int n)	{nf = n;}
	int get_nf() const 		{return nf;}
	void set_ne(const int n)	{ne = n;}
	int get_ne() const 		{return ne;}

	void add_vertex(MVertex v)      {vert_list.push_back(v); }
	MVertex vertex(int i) const	{return vert_list[i];}
	
	void add_face(const MFace f)	{face_list.push_back(f); }
	MFace face(int i) const		{return face_list[i];}

	vector<MVertex> vert_list;
	vector<MEdge> edge_list;
	vector<MFace> face_list;

	inline void init()
	{
	   nv = 0;
	   nf = 0;
	   ne = 0;
	   vert_list.clear();
	   face_list.clear();
	}

private:
	int nv;
	int nf;
	int ne;
};

}

#endif //MDS_H
