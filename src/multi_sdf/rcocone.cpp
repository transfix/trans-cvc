/* $Id: rcocone.cpp 1527 2010-03-12 22:10:16Z transfix $ */

#include <multi_sdf/rcocone.h>

namespace multi_sdf
{

// -------------------------
// initialize
// -------------------------
// Initializes the attributes associated with the voronoi and
// delaunay simplices
// -------------------------
void
compute_cell_radius(Triangulation &triang)
{
    for(FCI cit = triang.finite_cells_begin();
	cit != triang.finite_cells_end(); cit ++)
    {
	    double r = CGAL::to_double( (cit->voronoi() - cit->vertex(0)->point())*
			                (cit->voronoi() - cit->vertex(0)->point()) );
	    if(isnan(r) || isinf(r))
	    {
		    cerr << "cell radius nan or inf " << endl;
		    cit->set_cell_radius(HUGE);
	    }
	    else if(r <= 0)
		    cit->set_cell_radius(0);
	    else
		    cit->set_cell_radius(r);

    }
}



// ------------------------------
// mark_big_tetrahedral_balls
// ------------------------------
// Mark the tetrahedra which have the circumspherical
// balls quite big compared to the distance to the 
// nearest neighbor for every corner of the tetrahedron
// ------------------------------
void
mark_big_tetrahedral_balls(Triangulation &triang,
		           const double &bb_ratio)
{

    // find the nearest neighbor distance(nnd) for each vertex.
    // nnd(v) = average of distance to first three nearest neighbors.
    for(FEI eit = triang.finite_edges_begin();
	eit != triang.finite_edges_end(); eit ++)
    {
	    Cell_handle ch = (*eit).first;
	    Vertex_handle u = ch->vertex((*eit).second);
	    Vertex_handle v = ch->vertex((*eit).third);

	    // find the length of the edge.
	    double le = CGAL::to_double((u->point() - v->point())*
			    		(u->point() - v->point()) );
	    CGAL_assertion(le > 0);
	    
	    for(int i = 0; i < (int)u->nnd_vector.size(); i ++)
	    {
		    if(u->nnd_vector[i] < le) continue;
		    // insert le into the ith position.
		    double push_val = le;
		    for(int j = i; j < (int)u->nnd_vector.size(); j ++)
		    {
			    double temp = u->nnd_vector[j];
			    u->nnd_vector[j] = push_val;
			    push_val = temp;
		    }
		    break;
	    }
	    
	    for(int i = 0; i < (int)v->nnd_vector.size(); i ++)
	    {
		    if(v->nnd_vector[i] < le) continue;
		    // insert le into the ith position.
		    double push_val = le;
		    for(int j = i; j < (int)v->nnd_vector.size(); j ++)
		    {
			    double temp = v->nnd_vector[j];
			    v->nnd_vector[j] = push_val;
			    push_val = temp;
		    }
		    break;
	    }
    }



    // for each cell check if the circumradius is big compared to
    // the nnd of atleast one of its corners.
    for(FCI cit = triang.finite_cells_begin();
	cit != triang.finite_cells_end(); cit ++)
    {
	    double r = cit->cell_radius();
	    for(int i = 0; i < 4; i ++)
	    {
		    double nnd = (1./3.)*(cit->vertex(i)->nnd_vector[0] +
				          cit->vertex(i)->nnd_vector[1] +
					  cit->vertex(i)->nnd_vector[2] );
		    // if the radius is sufficiently big i.e.,
		    // radius is bigger than a certain fraction of 
		    // the average nearest neighbor distance
		    if(r > bb_ratio * nnd)
			    cit->set_bb(true);
	    }
    }
}







// -------------------------------
// identify_deep_intersection
// -------------------------------
// Identify which two big balls really 
// intersect deeply.
// -------------------------------
void
identify_deep_intersection(Triangulation &triang,
			   const double theta_ff,
			   const double theta_if)
{


   for(AFI fit = triang.all_facets_begin();
	fit != triang.all_facets_end(); fit ++)
   {
	   Cell_handle c[2]; int id[2];
	   c[0] = (*fit).first; id[0] = (*fit).second;
	   c[1] = c[0]->neighbor(id[0]); id[1] = c[1]->index(c[0]);

	   // there are three cases depending on whether the tetrahedra
	   // are finite or infinite.

	   // case 1 : both infinite .. not needed .. we assume they 
	   //          all are outside, big and intersect deeply.
	   if(triang.is_infinite(c[0]) && triang.is_infinite(c[1]))
		   continue;

	   // case 2 : one infinite and one finite.
	   if(triang.is_infinite(c[0]) || triang.is_infinite(c[1]))
	   {
		   // assign c[0] the finite tetrahedron.
		   if(triang.is_infinite(c[0]))
		   {
			   Cell_handle temp = c[0];
			   c[0] = c[1];
			   c[1] = temp;
			   id[0] = c[0]->index(c[1]);
			   id[1] = c[1]->index(c[0]);
		   }
		   CGAL_assertion(! triang.is_infinite(c[0]) &&
				    triang.is_infinite(c[1]) );

		   if(! c[0]->bb()) continue;

		   Point p[4];
		   p[0] = c[0]->vertex(id[0])->point();
		   p[1] = c[0]->vertex((id[0]+1)%4)->point();
		   p[2] = c[0]->vertex((id[0]+2)%4)->point();
		   p[3] = c[0]->vertex((id[0]+3)%4)->point();
		   // if the cell is not flowing into the infinite cell
		   // then the intersection is not deep.
		   if(CGAL::orientation(p[1],p[2],p[3],p[0]) ==
		      CGAL::orientation(p[1],p[2],p[3],c[0]->voronoi()))
			   continue;

		   // check the deepness of the intersection.
		   // let r_f = circumradius of the facet
		   //     r_t = circumradius of the tetrahedron
		   // then sin(theta_current) = r_f/r_t.
		   // so deep intersection <=> r_f/r_t < sin(theta).    
		   double r_f = sq_cr_tr_3(p[1],p[2],p[3]);
		   double r_t = c[0]->cell_radius();

		   if(r_f/r_t < sin(theta_if)*sin(theta_if))
		   {
		   	c[0]->set_deep_int(id[0], true);
		   	c[1]->set_deep_int(id[1], true);
		   }

		   continue;
	   }

	   // case 3 : both finite.
	   
	   if(!c[0]->bb() || !c[1]->bb())
		   continue;

	   // check if the circumscribing balls of the two cells 
	   // intersect deeply.
	   // Let (p,r_p) and (q,r_q) be the two balls. 
	   // Check the expression : ||p-q||*||p-q|| - r_p*r_p - r_q*r_q.
	   Point p = c[0]->voronoi();
	   Point q = c[1]->voronoi();
	   double r_p = c[0]->cell_radius();        
	   double r_q = c[1]->cell_radius();        
	   double d = CGAL::to_double((p-q)*(p-q)); // = ||p-q||*||p-q||

	   // we call an intersection deep when the angle made
	   // is less than 60 degree.

	   double r1 = (r_p > r_q) ? r_p : r_q;
	   double r2 = r_p + r_q - r1;

	   CGAL_assertion(r1 >= r2);

	   // d <= r1 + r2 - 2*sqrt(r1*r2)*cos(theta).
	   // when theta = 60 degree.
	   // d <= r1 + r2 - sqrt(r1*r2).
	   if(d <= r1 + r2 - 2*sqrt(r1*r2)*cos(theta_ff))
	   {
		   c[0]->set_deep_int(id[0], true);
		   c[1]->set_deep_int(id[1], true);
	   }
   }
}

// -------------------------------
// walk
// ------
// Walk among the cells to collect 
// the connected components.
// -------------------------------
void
walk(Triangulation &triang)
{

   for(ACI cit = triang.all_cells_begin();
	cit != triang.all_cells_end(); cit ++)
   {
	   cit->visited = false;
	   if(triang.is_infinite(cit))
		   cit->outside = true;
	   else
		   cit->outside = false;
   }

   // start from a cell which has big circumscribing ball
   // and which is flowing into infinity.
   for(FCI cit = triang.finite_cells_begin();
	cit != triang.finite_cells_end(); cit ++)
   {
	   if(cit->visited) continue;
	   if(!cit->bb()) continue;

	   bool start_flag = false;
	   Cell_handle start_cell = cit;

	   for(int i = 0; i < 4; i ++)
	   {
		   if(!triang.is_infinite(cit->neighbor(i)))
			   continue;
		   // so this cell has a face in the convex hull.
		   // check if the circumscribing ball of the cell
		   // intersects the infinite ball deeply.
		   if(cit->deep_int(i))
			   start_flag = true;
	   }
	   
	   if(! start_flag) continue;

	   vector<Cell_handle> cell_stack;
	   cell_stack.push_back(start_cell);
	   start_cell->visited = true;

	   while(!cell_stack.empty())
	   {
		   Cell_handle c = cell_stack.back();
		   cell_stack.pop_back();

		   c->outside = true;

		   // get the other neighboring cells into the stack.
		   // restrictions :
		   // 1. cell has to have big circumscribing balls.
		   // 2. cell must not be visited already.
		   // 3. deep_intersection flag has to be true.
		   // 4. cell has to be finite.
		   for(int i = 0; i < 4; i ++)
		   {
			   if(c->neighbor(i)->visited) continue;
			   if(triang.is_infinite(c->neighbor(i))) continue;
			   if(!c->neighbor(i)->bb()) continue;
			   CGAL_assertion(c->deep_int(i) ==
					  c->neighbor(i)->deep_int(c->neighbor(i)->index(c)));
			   if(!c->deep_int(i)) continue;

			   cell_stack.push_back(c->neighbor(i));
			   c->neighbor(i)->visited = true;
		   }
	   }
   }

   
   for(FCI cit = triang.finite_cells_begin();
	cit != triang.finite_cells_end(); cit ++)
   {
	   
	   if(! cit->outside) 
	   {
		   for(int i = 0; i < 4; i ++)
		   {
			   if(triang.is_infinite(cit->neighbor(i)))
			   {
				   cit->vertex((i+1)%4)->set_on_smooth_surface(true);
				   cit->vertex((i+2)%4)->set_on_smooth_surface(true);
				   cit->vertex((i+3)%4)->set_on_smooth_surface(true);
			   }
		   }
		   continue;
	   }
	   for(int i = 0; i < 4; i ++)
		   cit->vertex(i)->set_on_smooth_surface(true);
   }

}

// -------------------------------
// write_tcip
// ---------
// Write the points in the outside tetrahedron. 
// These will be input to tight cocone.
// -------------------------------
void
write_tcip(const Triangulation &triang,
           const char* file_prefix)
{
    ofstream fout;
    char filename[100];
    
    strcat(strcpy(filename, file_prefix), ".tcip");
    fout.open(filename);

    if(!fout)
    {
	    cerr << "can not open output file " << filename << endl;
	    exit(1);
    
    }
    
    for(FVI vit = triang.finite_vertices_begin();
	vit != triang.finite_vertices_end(); vit ++)
	    if(vit->on_smooth_surface())
	    	fout << vit->point() << endl;
}

void
robust_cocone(const double bb_ratio,
	      const double theta_ff,
	      const double theta_if,
	      Triangulation &triang,
	      const char *outfile_prefix)
{
  // computation of cell radius.
  compute_cell_radius(triang);
  cerr << ".";

  // mark big balls.
  mark_big_tetrahedral_balls(triang, bb_ratio);
  cerr << ".";

  // identify deep intersection between the balls.
  identify_deep_intersection(triang, theta_ff, theta_if);
  cerr << ".";

  // walk among the cells to collect the components.
  walk(triang);
  cerr << ".";

  // write the points as input to the tight cocone.
  write_tcip(triang, outfile_prefix);
  cerr << ".";

}

}
