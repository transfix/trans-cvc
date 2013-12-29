
#include <multi_sdf/op.h>

namespace multi_sdf
{

extern vector<double> bounding_box;

void
draw_ray(const Ray& ray,
         const double& r, 
         const double& g, 
         const double& b, 
         const double& a, 
         ofstream& fout)
{
   fout << "{OFF" << endl;
   fout << "2 1 0" << endl;
   fout << ray.source() << endl;
   fout << (ray.source() - CGAL::ORIGIN) + ray.to_vector() << endl;
   fout << "2\t0 1 " << r << " " << g << " " << b << " " << a << endl;
   fout << "}" << endl;
}

void
draw_segment(const Segment& segment, 
             const double& r, 
             const double& g, 
             const double& b, 
             const double& a,
             ofstream& fout)
{
   fout << "{OFF" << endl;
   fout << "2 1 0" << endl;
   fout << segment.point(0) << endl;
   fout << segment.point(1) << endl;
   fout << "2\t0 1 " << r << " " << g << " " << b << " " << a << endl;
   fout << "}" << endl;
}


void
draw_poly(const vector<Point>& poly,
          const double& r, 
          const double& g, 
          const double& b, 
          const double& a,
          ofstream& fout)
{
   fout << "{OFF" << endl;
   fout << (int)poly.size() << " 1 0" << endl;
   for(int i = 0; i < (int)poly.size(); i ++)
      fout << poly[i] << endl;
   fout << (int)poly.size() << "\t";
   for(int i = 0; i < (int)poly.size(); i ++)
      fout << i << " ";
   fout << r << " " << g << " " << b << " " << a << endl;
   fout << "}" << endl;
}


void
draw_VF(const Triangulation& triang,
        const Edge& dual_e, 
        const double& r, 
        const double& g, 
        const double& b, 
        const double& a,
        ofstream& fout)
{
   Facet_circulator fcirc = triang.incident_facets(dual_e);
   Facet_circulator begin = fcirc;
   vector<Point> vvset;
   do
   {
      Cell_handle cc = (*fcirc).first;
      vvset.push_back(cc->voronoi());
      fcirc ++;
   } while(fcirc != begin);

   fout << "{OFF" << endl;
   fout << (int)vvset.size() << " 1 0" << endl;
   for(int i = 0; i < (int)vvset.size(); i ++)
      fout << vvset[i] << endl;
   fout << (int)vvset.size() << "\t";
   for(int i = 0; i < (int)vvset.size(); i ++)
      fout << i << " ";
   fout << r << " " << g << " " << b << " " << a << endl;
   fout << "}" << endl;
}

void
draw_tetra(const Cell_handle& cell,
           const double& r, 
           const double& g, 
           const double& b, 
           const double& a,
           ofstream& fout)
{
   fout << "{OFF" << endl;
   fout << "4 4 0" << endl;
   for(int i = 0; i < 4; i ++)
      fout << cell->vertex(i)->point() << endl;
   for(int i = 0; i < 4; i ++)
      fout << "3\t" << (i+1)%4 << " " 
                    << (i+2)%4 << " "
                    << (i+3)%4 << " "
                    << r << " " << g << " " << b << " " << a << endl;
   fout << "}" << endl;
}

void 
write_wt( const Triangulation &triang,
	  const char* file_prefix) 
{
  char filename[100];
  strcat(strcpy(filename, file_prefix), ".surf");
  ofstream fout;
  fout.open(filename);
  if(! fout)
  {
     cerr << "Can not open " << filename << " for writing. " << endl;
     exit(1);
  }

  // Count number of facets on the surface.
  int num_facets = 0;
  for ( FFI fit = triang.finite_facets_begin();
	fit != triang.finite_facets_end(); ++fit) 
     if((*fit).first->cocone_flag((*fit).second) )
        num_facets ++;

  fout <<"OFF" << endl;
  fout << triang.number_of_vertices()
       << " " << num_facets << " 0" << endl;

  // Write the vertices.
  for ( FVI vit = triang.finite_vertices_begin();
	vit != triang.finite_vertices_end(); ++vit) 
     fout << vit->point() << endl;
  // Write the facets.
  for ( FFI fit = triang.finite_facets_begin();
	fit != triang.finite_facets_end(); ++fit) 
  {
     Cell_handle c[2] = {(*fit).first, (*fit).first->neighbor((*fit).second)};
     int id[2] = {c[0]->index(c[1]), c[1]->index(c[0])};
     if(! c[0]->cocone_flag(id[0]) )
        continue;
     CGAL_assertion( c[0]->bdy[id[0]] && c[1]->bdy[id[1]] );
     CGAL_assertion( c[0]->outside != c[1]->outside );

     Vertex_handle vh[3] = { c[0]->vertex((id[0]+1)%4),
                             c[0]->vertex((id[0]+2)%4),
                             c[0]->vertex((id[0]+3)%4) };
     if( ! c[0]->outside )
        if( CGAL::is_negative( Tetrahedron( vh[0]->point(), vh[1]->point(), vh[2]->point(), 
                                            c[0]->vertex(id[0])->point() ).volume() ) )
           fout << "3\t" << vh[0]->id << " " << vh[1]->id << " " << vh[2]->id << " ";
        else
           fout << "3\t" << vh[1]->id << " " << vh[0]->id << " " << vh[2]->id << " ";
     else
        if( CGAL::is_negative( Tetrahedron( vh[0]->point(), vh[1]->point(), vh[2]->point(), 
                                            c[1]->vertex(id[1])->point() ).volume() ) )
           fout << "3\t" << vh[0]->id << " " << vh[1]->id << " " << vh[2]->id << " ";
        else
           fout << "3\t" << vh[1]->id << " " << vh[0]->id << " " << vh[2]->id << " ";
     fout << endl;
  }
  fout.close();
  return;
}


void 
write_iobdy( const Triangulation &triang,
	     const char* file_prefix) 
{

  char filename[100];
  strcat(strcpy(filename, file_prefix), ".io");

  ofstream fout;
  fout.open(filename);

  if(! fout)
  {
	  cerr << "Can not open " << filename << " for writing. " << endl;
	  exit(1);
  }

  // Count
  int num_facets = 0;
  for ( FFI fit = triang.finite_facets_begin();
	fit != triang.finite_facets_end(); ++fit) {
    Cell_handle ch = (*fit).first;
    int id = (*fit).second;
    if(ch->outside != ch->neighbor(id)->outside)
	    num_facets ++;
  }
  
  // The header of the output file

  fout <<"OFF";
  fout <<"  " << triang.number_of_vertices(); //The number of points
  fout <<" " << num_facets; //The number of facets
  fout <<" 0" << endl;

  // Write the vertices.
  for ( FVI vit = triang.finite_vertices_begin();
	vit != triang.finite_vertices_end(); ++vit) 
    fout << vit->point() << endl;
  
  for ( FFI fit = triang.finite_facets_begin();
	fit != triang.finite_facets_end(); ++fit) 
  {
    Cell_handle ch = (*fit).first;
    int id = (*fit).second;

    if(ch->outside == ch->neighbor(id)->outside)
	    continue;
    
    fout << " 3\t";
    for (int i = 1; i <= 3; i++)
    	fout << " " << ch->vertex((id+i)%4)->id;
    fout << "\t " <<  "1 1 1 1 \n"; ;
  }

  fout.close();
}

}

