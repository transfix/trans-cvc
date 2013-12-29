/* $Id: util.h 1527 2010-03-12 22:10:16Z transfix $ */

#ifndef UTIL_H
#define UTIL_H

#include <multi_sdf/datastruct.h>

namespace multi_sdf
{
  double 
    cosine( const Vector& v, const Vector& w);
  
  void
    normalize(Vector& v);

  double
    length_of_seg(const Segment& s);

  bool
    is_obtuse(const Point& p0, const Point& p1, const Point& p2);

  int 
    find_third_vertex_index( const Facet& f, Vertex_handle v, Vertex_handle w);

  int 
    edge_index( const int facet_index, const int first_vertex_index, 
		const int second_vertex_index);

  void 
    vertex_indices( const int facet_index, const int edge_index,
		    int& first_vertex, int& second_vertex);

  bool
    is_same_side_of_ray(const Point& p0, const Point& p1,
			const Point& a, const Point& b);

  void
    is_contained_in_inf_tr(const Point& p0, const Point& p1, const Point& p2,
			   const Point& a, const Point& b,
			   const vector<int>& coincidence_vector,
			   bool* contained);
  bool
    is_outside_bounding_box(const Point& p, 
			    const vector<double>& bounding_box);

  bool
    is_outside_bounding_box(const vector<Point>& points, 
			    const vector<double>& bounding_box);

  bool
    is_VF_outside_bounding_box(const Triangulation& triang,
			       const Edge& e,
			       const vector<double>& bounding_box);
  bool
    is_outside_VF(const Triangulation& triang, 
		  const Edge& e);

  bool
    is_inside_VF(const Triangulation& triang, 
		 const Edge& e);

  bool
    is_surf_VF(const Triangulation& triang, 
	       const Edge& e);

  bool
    is_cospherical_pair(const Triangulation& triang, const Facet& f);

  bool
    identify_cospherical_neighbor(Triangulation &triang);

  // void
  // mark_VF_visited(Triangulation& triang,
  //                 Cell_handle& c, int uid, int vid);

  // void
  // mark_VF_medax_flag(Triangulation& triang,
  //                    Cell_handle& c, int uid, int vid, const bool& b);

  bool          
    is_inf_VF(const Triangulation& triang,
	      const Cell_handle& c, const int uid, const int vid);

  bool
    is_there_any_common_element(const vector<int>& vec1, const vector<int>& vec2);

}

#endif // UTIL_H

