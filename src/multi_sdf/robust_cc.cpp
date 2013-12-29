
#include <multi_sdf/robust_cc.h>

namespace multi_sdf
{

// -----------------------------------------------------------------------
// Functions to handle the degeneracy and potential problem in voronoi 
// computation
// -----------------------------------------------------------------------

///////////////////////////////////////////////////////////////
// these next functions are used to get the voronoi point if
// points are coplanar. Coplanarity makes CGAL very sad. :-(

//a flag that will be set if the call to circumsphere fails
static bool cgal_failed;

//the handler that will flag the operation as having failed
//oh, what the heck, let's make it static while we're at it
static void failure_func(
        const char *type,
        const char *exp,
        const char *file,
        int line,
        const char *expl)
{
        //bad CGAL!
        cgal_failed = true;

        return;
}



// -------------------------------------------------------------
// nondg_voronoi_point
// -------------------------------------------------------------
// This function will not assume any degeneracy and compute the
// voronoi point. If the computation fails it will set a boolean
// variable passed as a parameter by which the calling function 
// will know that the computation was not correct.
// -------------------------------------------------------------
Point
nondg_voronoi_point(const Point &a,
                    const Point &b,
                    const Point &c,
                    const Point &d,
                    bool &is_correct_computation)
{
 
        
        //we tell CGAL to call our function if there is a problem
	//we also tell it not to die if things go haywire
        //CGAL::Failure_function old_ff = CGAL::set_error_handler(failure_func);
        //CGAL::Failure_behaviour old_fb = CGAL::set_error_behaviour(CGAL::CONTINUE);

	// be optimistic :-)
	//this is a global
        cgal_failed = false;

        Point cc = CGAL::circumcenter( a,b,c,d );

	is_correct_computation = !cgal_failed;

        if(cgal_failed)
        {
            // set cc a junk value. It's not the duty of
            // this function to find out what the correct 
            // circumcenter should be.

            cc = CGAL::ORIGIN;
        }
	//put everything back the way we found it,
	//CGAL::set_error_handler(old_ff);
        //CGAL::set_error_behaviour(old_fb);

	//see ya! I'm outta here!
	
        if(cc == CGAL::ORIGIN) cerr << "<nondg_vp : returning cc as ORIGIN>" << endl;
        return cc;
}

// -------------------------------------------------------------
// dg_voronoi_point
// -------------------------------------------------------------
// This function computes a voronoi point when some degeneracies
// have been discovered about the four points of a tetrahedron.
// 
// In this function we assume that the degeneracies occured 
// because of the coplanarity. We approximate the circumcenter
// of the tetrahedron by the circumcenter of one of the triangular
// facets. 
// -------------------------------------------------------------
Point
dg_voronoi_point(const Point &a,
                 const Point &b,
                 const Point &c,
                 const Point &d,
                 bool &is_correct_computation)
{

        // First, we check if our assumption is correct.
        // This is more of a debugging purpose than of 
        // actual computation.
        Tetrahedron t = Tetrahedron(a,b,c,d);
        if( ! t.is_degenerate() )
        {
              // cerr << "error in the assumption of coplanarity." << endl;
        /*
        // debug
        cout << "{OFF " << endl;
        cout << "4 4 0" << endl;
        cout << a << "\n" << b << "\n" << c << "\n" << d << endl;
        cout << "3\t0 1 2" << endl;
        cout << "3\t0 2 3" << endl;
        cout << "3\t0 3 1" << endl;
        cout << "3\t1 2 3" << endl;
        cout << "}" << endl;
        // end debug
        */
        }

        // Approximate the circumcenter of the tetrahedron with that of
        // one of its triangular facets. The facet should not be collinear.
        // The following boolean variable will keep track if we have found
        // a valid circumcenter (of a triangle) to replace that of a
        // tetrahedron.
        Point cc = CGAL::ORIGIN;
        is_correct_computation = false;
        Point p[4] = {a,b,c,d};
        for(int i = 0; i < 4; i ++)
        {
             // first check if the facet is degenerate or not.
             Triangle_3 t = Triangle_3(p[(i+1)%4], p[(i+2)%4], p[(i+3)%4]);
             if( t.is_degenerate() ) continue;

             // since we found a non-degenerate triangle we can now compute
             // its circumcenter and we will be done.
             cc = nondg_cc_tr_3(p[(i+1)%4], p[(i+2)%4], p[(i+3)%4], is_correct_computation);

             if( is_correct_computation )
                   break;
        }

        if( is_correct_computation )
        {
             if(cc == CGAL::ORIGIN) cerr << "<dg_vp : returning cc as ORIGIN>" << endl;
             return cc;
        }

        cerr << "four points are colinear. " << endl;

        // for the time being, I just average the four points. What should be
        // the circumcenter of a tetrahedron whose four points are collinear ?
        cc = CGAL::ORIGIN +
             ( 0.25*Vector
                       (
                          a[0]+b[0]+c[0]+d[0],
                          a[1]+b[1]+c[1]+d[1],
                          a[2]+b[2]+c[2]+d[2]
                       )
             );

        if(cc == CGAL::ORIGIN) cerr << "<dg_vp : returning cc as ORIGIN>" << endl;
        return cc;
}




// -------------------------------------------------------------
// nondg_cc_tr_3
// -------------------------------------------------------------
// This function computes the circumcenter of a triangle in 3D.
// It doesn't assume any degeneracy but in case it happens it 
// will set a flag which has been passed as parameter.
// -------------------------------------------------------------
Point
nondg_cc_tr_3(const Point &a,
              const Point &b,
              const Point &c,
              bool &is_correct_computation)
{
        
        //we tell CGAL to call our function if there is a problem
	//we also tell it not to die if things go haywire
        //CGAL::Failure_function old_ff = CGAL::set_error_handler(failure_func);
        //CGAL::Failure_behaviour old_fb = CGAL::set_error_behaviour(CGAL::CONTINUE);

	// be optimistic :-)
	//this is a global
        cgal_failed = false;

        Point cc = CGAL::circumcenter( a,b,c );

	is_correct_computation = !cgal_failed;

        if(cgal_failed)
        {
            // set cc a junk value. It's not the duty of
            // this function to find out what the correct 
            // circumcenter should be.

            cc = CGAL::ORIGIN;
        }
	//put everything back the way we found it,
	//CGAL::set_error_handler(old_ff);
        //CGAL::set_error_behaviour(old_fb);

        return cc;
}


// -------------------------------------------------------------
// cc_tr_3
// -------------------------------------------------------------
// This function computes the circumcenter of a triangle in 3D.
// It doesn't assume any degeneracy and calls nondg function to
// do the job. The nondg function checks if there is any deg and
// sets a flag. Then this function will approximate the cc as the
// average of the three points passed as the parameters.
// -------------------------------------------------------------
Point
cc_tr_3(const Point &a,
        const Point &b,
        const Point &c)
{
        bool is_correct_computation = true;
        
        // call nondg function which will set is_correct_computation
        // depending on the degeneracy.
        Point cc = nondg_cc_tr_3(a,b,c,is_correct_computation);

        if( is_correct_computation ) return cc;
        cerr << "< Bad cc of triangle > ";

        // if the three points were degenerate (collinear) then 
        // we approximate the circumcenter by the average of the
        // three points.
        cc = CGAL::ORIGIN +
             ( (1./3.)*Vector
                       (
                          a[0]+b[0]+c[0],
                          a[1]+b[1]+c[1],
                          a[2]+b[2]+c[2]
                       )
             );
        
        return cc;
}


 
// -------------------------------------------------------------
// sq_cr_tr_3
// -------------------------------------------------------------
// This function computes the squared circumradius of a triangle in 3D.
// It calls cc_tr_3 to get the circumcenter and then computes the
// radius.
// -------------------------------------------------------------
double
sq_cr_tr_3(const Point &a,
           const Point &b,
           const Point &c)
{
        Point cc = cc_tr_3(a,b,c);
        double r = CGAL::to_double((cc - a)*(cc - a));

        CGAL_assertion( (r > 0) &&
                        (! isnan(r)) &&
                        (! isinf(r)) );

        return r;
}

Point
circumcenter(const Facet& f)
{
   Cell_handle c = f.first;
   int id = f.second;
   Point p[3];
   for(int i = 0; i < 3; i ++)
      p[i] = c->vertex((id + (i+1))%4)->point();
   return cc_tr_3(p[0], p[1], p[2]);
}

double
circumradius(const Facet& f)
{
   Point c = circumcenter(f);
   Point p = f.first->vertex((f.second+1)%4)->point();
   return sqrt(CGAL::to_double((p-c)*(p-c)));
}

}
