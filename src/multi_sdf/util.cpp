/* $Id: util.cpp 1527 2010-03-12 22:10:16Z transfix $ */

#include <multi_sdf/util.h>

namespace multi_sdf
{

// ---------------------------------------
// cosine
// -------------
// Given two vectors v and w find the cosine
// of the angle between v and w
// ---------------------------------------
double 
cosine( const Vector& v, const Vector& w) 
{
  return CGAL::to_double( v * w) / 
    sqrt( CGAL::to_double( v * v) * CGAL::to_double( w * w));
}

// -----------
// normalize
// -----------
// Given a vector v, normalize it
// -----------
void
normalize(Vector& v)
{
   v = (1./(sqrt(CGAL::to_double(v*v))))*v;
   return;
}

// check if the angle <p0,p2,p1 > 90 degree
bool
is_obtuse(const Point& p0, const Point& p1, const Point& p2)
{
   Vector v0 = p0 - p2;
   Vector v1 = p1 - p2;
   return (CGAL::to_double(v0 * v1) < 0);
}

double
length_of_seg(const Segment& s)
{
   return sqrt(CGAL::to_double((s.point(0) - s.point(1))*(s.point(0) - s.point(1))));
}


// Find the index of the third vertex for a facet where this vertex 
// is neither *v nor *w.
int 
find_third_vertex_index( const Facet& f, Vertex_handle v, Vertex_handle w) {
  int id = f.second;
  for ( int i =1; i <= 3; ++i) {
    if ( f.first->vertex((id+i)%4) != v && 
	 f.first->vertex((id+i)%4) != w) 
      return (id+i)%4;
  }
  return -1;
}

// Compute the index of an edge. The index of an edge in a facet
// is defined as the index of a facet in cell, i.e. it is the index
// of the opposite vertex. The arguments are the facet index of the
// facet witth respect to the cell and the indices of the vertices
// incident to the edge also with respect to the cell.
int 
edge_index( const int facet_index, const int first_vertex_index, 
	    const int second_vertex_index) {
  return 6 - facet_index - first_vertex_index - second_vertex_index;
}


// ---------------------------------------
// vertex_indices
// --------------
// Compute the indices of the vertices 
// incident to an edge. Given, facet_id
// and edge_id.
// ---------------------------------------
void 
vertex_indices( const int facet_index, const int edge_index,
		int& first_vertex, int& second_vertex) {
  if ( (facet_index == 0 && edge_index == 1) ||
       (facet_index == 1 && edge_index == 0)) {
    first_vertex = 2; second_vertex = 3;
  } else if ( (facet_index == 0 && edge_index == 2) ||
	      (facet_index == 2 && edge_index == 0)) {
    first_vertex = 1; second_vertex = 3;
  } else if ( (facet_index == 0 && edge_index == 3) ||
	      (facet_index == 3 && edge_index == 0)) {
    first_vertex = 1; second_vertex = 2;
  } else if ( (facet_index == 1 && edge_index == 2) ||
	      (facet_index == 2 && edge_index == 1)) {
    first_vertex = 0; second_vertex = 3;
  } else if ( (facet_index == 1 && edge_index == 3) ||
	      (facet_index == 3 && edge_index == 1)) {
    first_vertex = 0; second_vertex = 2;
  } else if ( (facet_index == 2 && edge_index == 3) ||
	      (facet_index == 3 && edge_index == 2)) {
    first_vertex = 0; second_vertex = 1;
  }
}

// ------------------------
// is_same_side_of_ray
// ------------------------
// Given a ray and two points a and b find out if a and 
// b is in the same side of the ray.
// 
// Assumption: The ray, a and b are in the same plane.
// 
// [*] Instead of the Ray we will work on a pair of points
// (p0, p1) where p1 is dummy and can be shifted. [*]
// -------------------------
bool
is_same_side_of_ray(const Point& p0, const Point& p1,
                    const Point& a, const Point& b)
{
   return (CGAL::to_double(CGAL::cross_product(p1 - p0, a - p1)*
                           CGAL::cross_product(p1 - p0, b - p1)) >= 0);
}

// --------------------------
// is_contained_in_inf_tr
// --------------------------
// Given two rays and a segment find out if the segment
// is fully/partially contained in the infinite triangle.
// 
// Parameters: ray 1 = (p0, p1), ray 2 = (p0, p2). 
// Both p1 and p2 are just two points on the rays.
//             segment = (a, b)
// [*] Note: a or b can coincide with p1 or p2.
// coincidence_vector records that information.
// coincidence_vector[0] = {0,1,2} - matches with none, p1, p2.
// coincidence_vector[1] = {0,1,2} - matches with none, p1, p2.
//             two booleans co0 and co1 indicate that.
// Result will be stored in contained[].
// contained[0] says if a is within the fan.
// contained[1] says if b is within the fan.
// If both entries true the segment is fully contained.
// If both entries false the segment is fully outside.
// Otherwise it's partially contained.
// ---------------------------
void
is_contained_in_inf_tr(const Point& p0, const Point& p1, const Point& p2,
                       const Point& a, const Point& b,
                       const vector<int>& coincidence_vector,
                       bool* contained)
{
   if( coincidence_vector[0] == 0 && coincidence_vector[1] == 0)
   {
      // point a is within the opening (or visible by p0) if 
      // it is in the same side of ray1 (p0, p1) as p2 is and
      // it is in the same side of ray2 (p0, p2) as p1 is.
      contained[0] = is_same_side_of_ray(p0, p1, p2, a) &&
                     is_same_side_of_ray(p0, p2, p1, a);
      // same for b.
      contained[1] = is_same_side_of_ray(p0, p1, p2, b) &&
                     is_same_side_of_ray(p0, p2, p1, b);
      return;
   }
   else
   {
      // point a = p1 and b doesn't match.
      if(coincidence_vector[0] == 1 && coincidence_vector[1] == 0)
      {
         if( ! is_same_side_of_ray(p0, p1, p2, b) )
            contained[0] = contained[1] = false;
         else
         {
            contained[0] = true;
            if( is_same_side_of_ray(p0, p2, p1, b) )
               contained[1] = true;
            else
               contained[1] = false;
         }
         return;
      }
      // point a = p2 and b doesn't match.
      if(coincidence_vector[0] == 2 && coincidence_vector[1] == 0)
      {
         if( ! is_same_side_of_ray(p0, p2, p1, b) )
            contained[0] = contained[1] = false;
         else
         {
            contained[0] = true;
            if( is_same_side_of_ray(p0, p1, p2, b) )
               contained[1] = true;
            else
               contained[1] = false;
         }
         return;
      }
      // point b = p1 and a doesn't match.
      if(coincidence_vector[0] == 0 && coincidence_vector[1] == 1)
      {
         if( ! is_same_side_of_ray(p0, p1, p2, a) )
            contained[0] = contained[1] = false;
         else
         {
            contained[1] = true;
            if( is_same_side_of_ray(p0, p2, p1, a) )
               contained[0] = true;
            else
               contained[0] = false;
         }
         return;
      }
      // point b = p2 and a doesn't match.
      if(coincidence_vector[0] == 0 && coincidence_vector[1] == 2)
      {
         if( ! is_same_side_of_ray(p0, p2, p1, a) )
            contained[0] = contained[1] = false;
         else
         {
            contained[1] = true;
            if( is_same_side_of_ray(p0, p1, p2, a) )
               contained[0] = true;
            else
               contained[0] = false;
         }
         return;
      }
      // otherwise the segment(a,b) coincides with the segment(p1,p2).
      // the segment must be contained in the inf_tr (p0, p1, p2).
      contained[0] = true; contained[1] = true;
      return;
   }
   return;
}

// ------------------------
// is_outside_bounding_box
// ------------------------
// Given a bounding_box and a set of points
// find if any point is outside the bounding_box
// ------------------------
bool
is_outside_bounding_box(const Point& p, 
                       const vector<double>& bounding_box)
{
   if(CGAL::to_double(p.x()) < bounding_box[0] ||
      CGAL::to_double(p.x()) > bounding_box[1] ||
      CGAL::to_double(p.y()) < bounding_box[2] ||
      CGAL::to_double(p.y()) > bounding_box[3] ||
      CGAL::to_double(p.z()) < bounding_box[4] ||
      CGAL::to_double(p.z()) > bounding_box[5] )
      return true;
   return false;
}



// ------------------------
// is_outside_bounding_box
// ------------------------
// Given a bounding_box and a set of points
// find if any point is outside the bounding_box
// ------------------------
bool
is_outside_bounding_box(const vector<Point>& points, 
                        const vector<double>& bounding_box)
{
   for(int i = 0; i < (int)points.size(); i ++)
   {
      Point p = points[i];
      if( is_outside_bounding_box(p, bounding_box) ) return true;
   }
   return false;
}

// ------------------------
// is_VF_outside_bounding_box
// ------------------------
// Given a bounding_box and a VF
// find if any VV of VF is outside the bounding_box
// ------------------------
bool
is_VF_outside_bounding_box(const Triangulation& triang, 
                           const Edge& e,
                           const vector<double>& bounding_box)
{
   Facet_circulator fcirc = triang.incident_facets(e);
   Facet_circulator begin = fcirc;
   do{
      if(is_outside_bounding_box((*fcirc).first->voronoi(), bounding_box)) return true;
      fcirc ++;
   } while(fcirc != begin);
   return false;
}

bool
is_outside_VF(const Triangulation& triang, 
              const Edge& e)
{
   Facet_circulator fcirc = triang.incident_facets(e);
   Facet_circulator begin = fcirc;
   do{
      if( ! (*fcirc).first->outside) return false;
      fcirc ++;
   } while(fcirc != begin);
   return true;
}

bool
is_inside_VF(const Triangulation& triang, 
             const Edge& e)
{
   Facet_circulator fcirc = triang.incident_facets(e);
   Facet_circulator begin = fcirc;
   do{
      if((*fcirc).first->outside) return false;
      fcirc ++;
   } while(fcirc != begin);
   return true;
}

bool
is_surf_VF(const Triangulation& triang, 
           const Edge& e)
{
   Facet_circulator fcirc = triang.incident_facets(Edge(e));
   Facet_circulator begin = fcirc;
   bool in = false, out = false;
   do{
      Cell_handle cur_c = (*fcirc).first;
      if( cur_c->outside ) out = true;
      else in = true;
      // int cur_fid = (*fcirc).second;
      // if(cur_c->cocone_flag(cur_fid) ) return true;
      fcirc ++;
   } while(fcirc != begin);
   return (in && out);
}

bool          
is_inf_VF(const Triangulation& triang,
          const Cell_handle& c, const int uid, const int vid)
{
   Facet_circulator fcirc = triang.incident_facets(Edge(c,uid,vid));
   Facet_circulator begin = fcirc;
   do{
      Cell_handle cur_c = (*fcirc).first;
      if( triang.is_infinite( cur_c ) ) return true;
      fcirc ++;
   } while(fcirc != begin);
   return false;
}


bool
is_there_any_common_element(const vector<int>& vec1, const vector<int>& vec2)
{
   for(int i = 0; i < (int) vec1.size(); i ++)
      for(int j = 0; j < (int) vec2.size(); j ++)
         if( vec1[i] == vec2[j] ) return true;
   return false;
}

}
