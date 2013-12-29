/* $Id: init.cpp 1527 2010-03-12 22:10:16Z transfix $ */

#include <multi_sdf/init.h>
#include <multi_sdf/robust_cc.h>
#include <multi_sdf/op.h>

namespace multi_sdf
{

// ---------------------------------------------------------
// initialize
// ----------
// Initialize some of the attributes of the triangulation.
// ---------------------------------------------------------
void
initialize(Triangulation &triang)
{
   // set vertex id.
   int id = 0;
   for(FVI vit = triang.finite_vertices_begin();
	vit != triang.finite_vertices_end(); vit ++)
   {
	   vit->id = id++;
	   vit->visited = false;
    	   vit->bad = false;
    	   vit->bad_neighbor = false;
   }

   // set cell id.
   id = 0;
   for(ACI cit = triang.all_cells_begin();
	cit != triang.all_cells_end(); cit ++)
   {
	   cit->id = id++;
    	   cit->visited = false;
    	   cit->outside = false;
    	   cit->transp = false;

	   for(int id = 0 ; id < 4; id++)
	   {
      		cit->set_cocone_flag(id,false);
      		cit->neighbor(id)->set_cocone_flag(cit->neighbor(id)->index(cit),false);
      		cit->bdy[id] = false;
      		cit->opaque[id] = false;
      		for(int k = 0; k < 4; k ++)
	    		cit->umbrella_member[id][k] = -1;
    	   }
    	   
	   // set the convex hull points.
    	   if(! triang.is_infinite(cit)) continue;

    	   for(int i = 0; i < 4; i ++)
    	   {
	    	   if(! triang.is_infinite(cit->vertex(i))) continue;
	    	   cit->vertex((i+1)%4)->set_convex_hull(true);
	    	   cit->vertex((i+2)%4)->set_convex_hull(true);
	    	   cit->vertex((i+3)%4)->set_convex_hull(true);
    	   }
   }
}



// ---------------------------------------------------------
// compute_voronoi_vertex_and_cell_radius
// --------------------------------------
// Compute all the finite voronoi vertices and the circumradius
// of all the finite cells.
// ---------------------------------------------------------
void
compute_voronoi_vertex_and_cell_radius(Triangulation &triang)
{
  
   bool is_there_any_problem_in_VV_computation = false;
       
   for(FCI cit = triang.finite_cells_begin();
        cit != triang.finite_cells_end(); cit ++)
   {
           // For each cell call the circumcenter computation.
           // The function will set a boolean variable passed
           // as a parameter showing if the computation is correct.
           bool is_correct_computation = true;
           cit->set_voronoi(nondg_voronoi_point( cit->vertex(0)->point(),
                                                 cit->vertex(1)->point(),
                                                 cit->vertex(2)->point(),
                                                 cit->vertex(3)->point(),
                                                 is_correct_computation
                                               )
                           );
           is_there_any_problem_in_VV_computation |= !is_correct_computation;

           // The problem arises when the computation is not correct
           // that is, when is_correct_computation = false.

           // If is_correct_computation = true then we need to check
           // the radius. Sometimes, even if the circumcenter computation
           // is correct, radius of the cell is too big to fit in double
           // and that creates a problem. Sometimes, (possibly due to 
           // overflow) it is <= 0. In the next loop we will check for
           // the validity of radius and if that also works well we will
           // return. Otherwise, we will go the next part where the 
           // degeneracies are taken care of.
           
           if( is_correct_computation )
           {
                bool is_correct_radius = true;
                
                // check if the radius fits good in double.
          
                double r = CGAL::to_double((cit->voronoi() - cit->vertex(0)->point()) *
                                           (cit->voronoi() - cit->vertex(0)->point())
                                          );

                if( isnan(r) || isinf(r) )
                      is_correct_radius = false;
                if( r <= 0 )
                      is_correct_radius = false;

                // if it does, go back and collect the next cell
                // and do the same computation.
                if( is_correct_radius )
                      continue;
                is_there_any_problem_in_VV_computation |= !is_correct_radius;
           }

           cerr << " < bad > ";
           cit->set_dirty(true);

           // The flow comes here means either the voronoi computation 
           // is incorrect or the cell radius is junk.
           // Either case, we need to take measures.
           // Our assumption is 
           // This happens when the four points of the cell are coplanar
           // atleast that is reflected by the choice of arithmatic.
           // Our approach is to approximate the circumcenter of the cell
           // by the circumcenter of one of the triangular facets.

           cit->set_voronoi(dg_voronoi_point( cit->vertex(0)->point(),
                                              cit->vertex(1)->point(),
                                              cit->vertex(2)->point(),
                                              cit->vertex(3)->point(),
                                              is_correct_computation
                                            )
                           );

           // debug
           
           // Certain checks to make sure we are not setting anything bad
           // about the voronoi information in any cell.
           
           double cx = CGAL::to_double(cit->voronoi().x());
           double cy = CGAL::to_double(cit->voronoi().y());
           double cz = CGAL::to_double(cit->voronoi().z());

           // check 1 : the coordinates of the voronoi point fits well in double.
           CGAL_assertion( ( ( ! isnan(cx)) && ( ! isinf(cx)) ) &&
                           ( ( ! isnan(cy)) && ( ! isinf(cy)) ) &&
                           ( ( ! isnan(cz)) && ( ! isinf(cz)) ) 
                         );
           // check 2 : the radius of the cell fits well in double and is non-negative.

           double r = CGAL::to_double((cit->voronoi() - cit->vertex(0)->point()) *
                                      (cit->voronoi() - cit->vertex(0)->point())
                                     );
           CGAL_assertion( r > 0 &&
                           ! isnan(r) &&
                           ! isinf(r)
                         );

           // end debug
   }

   if(is_there_any_problem_in_VV_computation)
   {
      ofstream fout;
      fout.open("degen_tetra");
      for(FCI cit = triang.finite_cells_begin();
         cit != triang.finite_cells_end(); cit ++)
      {
         if(cit->dirty())
         {
            draw_tetra(cit, 1, 1, 0, 1, fout); continue;
         }
         if(cit->voronoi() == CGAL::ORIGIN)
         {
            draw_tetra(cit, 0, 1, 0, 1, fout); continue;
         }
      }
      fout.close();
   }
   return;
}

}

