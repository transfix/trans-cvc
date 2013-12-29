/* $Id: main.cpp 2539 2010-08-06 21:40:47Z transfix $ */

// Samrat
// Started Dec 4th, 2007.

#include <cstdio>
#include <multi_sdf/mesh_io.h>
#include <multi_sdf/sdf.h>
#include <multi_sdf/kdtree.h>
#include <multi_sdf/matrix.h>
#include <multi_sdf/dt.h>
#include <multi_sdf/multi_sdf.h>

using namespace std;

namespace multi_sdf
{

const int DIMX = 100;
const int DIMY = 100;
const int DIMZ = 100;
const double BB_SCALE = 0.2;

static inline void
construct_bbox(const Mesh& mesh, vector<double>& bbox)
{
   double x_min = HUGE, x_max = -HUGE,
          y_min = HUGE, y_max = -HUGE,
          z_min = HUGE, z_max = -HUGE;
   for(int i = 0; i < mesh.get_nv(); i ++)
   {
      if( mesh.vert_list[i].iso() ) continue;
      Point p = mesh.vert_list[i].point();
      
      // check x-span
      if(CGAL::to_double(p.x()) < x_min) 
         x_min = CGAL::to_double(p.x());
      if(CGAL::to_double(p.x()) > x_max) 
         x_max = CGAL::to_double(p.x());
      // check y-span
      if(CGAL::to_double(p.y()) < y_min) 
         y_min = CGAL::to_double(p.y());
      if(CGAL::to_double(p.y()) > y_max) 
         y_max = CGAL::to_double(p.y());
      // check z-span
      if(CGAL::to_double(p.z()) < z_min) 
         z_min = CGAL::to_double(p.z());
      if(CGAL::to_double(p.z()) > z_max) 
         z_max = CGAL::to_double(p.z());
   }
   bbox.push_back(x_min - BB_SCALE*(x_max-x_min));
   bbox.push_back(y_min - BB_SCALE*(y_max-y_min));
   bbox.push_back(z_min - BB_SCALE*(z_max-z_min));

   bbox.push_back(x_max + BB_SCALE*(x_max-x_min));
   bbox.push_back(y_max + BB_SCALE*(y_max-y_min));
   bbox.push_back(z_max + BB_SCALE*(z_max-z_min));
}

static inline void
assign_sdf_weight(Mesh& mesh, vector<double>& weights)
{
   // map the color of each facet to a weight (scalar).
   // for the time being the information is in file called "weights".
   ifstream fin;
   fin.open("weights");
   istream_iterator<double> input(fin);
   istream_iterator<double> beyond;
   double tw = 0;
   while(input != beyond) { tw += *input; weights.push_back(*input); input++; }  
   for(int i = 0; i < (int)weights.size(); i ++) weights[i] /= tw;
}

}


int main(int argc, char** argv)
{
  using namespace multi_sdf;

  std::string ifname;
  std::string ofname;
  FILE_TYPE in_ftype = OFF;
  bool read_color_opacity = false;
  bool is_uniform = false;

  int dimx = DIMX, dimy = DIMY, dimz = DIMZ;

  // Check commandline options.
  bool help = false;
  if(argc == 1)
     help = true;
  for (int i = 1; i < argc; i++) {
    if ( argc < 3 ) {
      help = true;
      break;
    }
    if ( (strcmp( "-h", argv[i]) == 0) || 
	 (strcmp( "-help", argv[i]) == 0)) {
      help = true;
      break;
    }
    else if(strcmp("-u", argv[i]) == 0 )
       is_uniform = true;
    else if(strcmp("-ca", argv[i]) == 0 )
       read_color_opacity = true;
    else if(strcmp("-off", argv[i]) == 0 ||
          strcmp("-OFF", argv[i]) == 0)
       in_ftype = OFF;
    else if(strcmp("-raw", argv[i]) == 0 ||
            strcmp("-RAW", argv[i]) == 0)
       in_ftype = RAW;
    else if(strcmp("-rawn", argv[i]) == 0 ||
            strcmp("-RAWN", argv[i]) == 0)
       in_ftype = RAWN;
    else if(strcmp("-rawc", argv[i]) == 0 ||
            strcmp("-RAWC", argv[i]) == 0)
       in_ftype = RAWC;
    else if(strcmp("-rawnc", argv[i]) == 0 ||
            strcmp("-RAWNC", argv[i]) == 0)
       in_ftype = RAWNC;
    else if ( i+1 >= argc) 
      help = true;
    else 
    {
       ifname = argv[i];
       ofname = argv[i+1];
       i++;
    }
  }

  if ( help) {
    cerr << "Usage: " << argv[0] 
	 << " [FILE_TYPE <default -off>] "
         << " infile outfile " << endl;
    exit( 1);
  }

  // read the annotated input file and
  Mesh mesh;
  cerr << "Reading input mesh ";
  read_labeled_mesh(mesh, ifname, in_ftype, read_color_opacity, is_uniform);
  cerr << "done." << endl;

  // build a bounding box around the input and store the
  // origin, span etc.
  vector<double> bbox; bbox.clear();
  construct_bbox(mesh, bbox);

  // construct a kd-tree of all the non-isolated mesh_vertices.
  vector<VECTOR3> points;
  vector<Point> pts;
  for(int i = 0; i < mesh.get_nv(); i ++)
  {
     if( mesh.vert_list[i].iso() ) continue;
     Point p = mesh.vert_list[i].point();
     pts.push_back(p);
     points.push_back(VECTOR3(CGAL::to_double(p.x()),
                              CGAL::to_double(p.y()),
                              CGAL::to_double(p.z())));
  }
  KdTree kd_tree(points, 20);
  kd_tree.setNOfNeighbours(1);

  // Now perform a reconstruction to build a tetrahedralized solid
  // with in-out marked.
  Triangulation triang;
  recon(pts, triang);

  // assign weight to each triangle.
  vector<double> weights; weights.clear();
  // assign_sdf_weight(mesh, weights); // comment out for uniform weight.

  cerr << "SDF ";
  try
    {
      //VolMagick::VolumeCache volcache;
      // VolMagickOpStatus status;

      // VolMagick::setDefaultMessenger(&status);

      VolMagick::Volume vol;

      vol.dimension(VolMagick::Dimension(dimx,dimy,dimz));
      vol.voxelType(VolMagick::Float);
      vol.boundingBox(VolMagick::BoundingBox(bbox[0],bbox[1],bbox[2],   // xmin,ymin,zmin,
                                             bbox[3],bbox[4],bbox[5])); // xmax,ymax,zmax

      for(unsigned int k=0; k<vol.ZDim(); k++)
      {
	for(unsigned int j=0; j<vol.YDim(); j++)
        {
	  for(unsigned int i=0; i<vol.XDim(); i++)
	    {
              double x = vol.XMin() + i*vol.XSpan();
              double y = vol.YMin() + j*vol.YSpan();
              double z = vol.ZMin() + k*vol.ZSpan();
              double fn_val = sdf(Point(x,y,z), mesh, weights, kd_tree, triang);
	      vol(i,j,k, fn_val);
	    }
        }
        cerr << ".";
      }

      VolMagick::createVolumeFile(ofname,
				  vol.boundingBox(),
				  vol.dimension(),
				  std::vector<VolMagick::VoxelType>(1, vol.voxelType()));
      VolMagick::writeVolumeFile(vol, ofname);
    }
  catch(VolMagick::Exception &e)
    {
      cerr << e.what() << endl;
    }
  catch(std::exception &e)
    {
      cerr << e.what() << endl;
    }

  cerr << "done." << endl;

  return 0;
}
