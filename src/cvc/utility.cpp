/*
  Copyright 2007-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolMagick.

  VolMagick is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolMagick is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <cvc/utility.h>
#include <cvc/app.h>

namespace CVC_NAMESPACE
{
  void calcGradient(std::vector<volume>& grad, const volume& vol, data_type vt)
  {
    thread_info ti(BOOST_CURRENT_FUNCTION);

    double dx,dy,dz,length;
    int i, j, k;

    volume gradx(vol.voxel_dimensions(),vt,vol.boundingBox());
    volume grady(vol.voxel_dimensions(),vt,vol.boundingBox());
    volume gradz(vol.voxel_dimensions(),vt,vol.boundingBox());

    //central differences algorithm
    for(k=0; k<int(vol.ZDim()); k++)
      {
	for(j=0; j<int(vol.YDim()); j++)
	  for(i=0; i<int(vol.XDim()); i++)
	    {
	      dx = (vol(MIN(i+1,int(vol.XDim())-1),j,k) - vol(MAX(i-1,0),j,k))/2.0;
	      dy = (vol(i,MIN(j+1,int(vol.YDim())-1),k) - vol(i,MAX(j-1,0),k))/2.0;
	      dz = (vol(i,j,MIN(k+1,int(vol.ZDim())-1)) - vol(i,j,MAX(k-1,0)))/2.0;
	      length = sqrt(dx*dx+dy*dy+dz*dz);
	      if(length>0.0)
		{
		  dx /= length;
		  dy /= length;
		  dz /= length;
		}
	      
	      switch(vt)
		{
		case UChar:
		  dx = dx*double((~char(0))>>1)+double((~char(0))>>1);
		  dy = dy*double((~char(0))>>1)+double((~char(0))>>1);
		  dz = dz*double((~char(0))>>1)+double((~char(0))>>1);
		  break;
		case UShort:
		  dx = dx*double((~short(0))>>1)+double((~short(0))>>1);
		  dy = dy*double((~short(0))>>1)+double((~short(0))>>1);
		  dz = dz*double((~short(0))>>1)+double((~short(0))>>1);
		  break;
		case UInt:
		  dx = dx*double((~int(0))>>1)+double((~int(0))>>1);
		  dy = dy*double((~int(0))>>1)+double((~int(0))>>1);
		  dz = dz*double((~int(0))>>1)+double((~int(0))>>1);
		  break;
		default: break;
		}
	      
	      gradx(i,j,k, dx);
	      grady(i,j,k, dy);
	      gradz(i,j,k, dz);
	    }

	cvcapp.threadProgress(float(k)/float(vol.ZDim()));
      }

    grad.clear();
    grad.push_back(gradx);
    grad.push_back(grady);
    grad.push_back(gradz);

    cvcapp.threadProgress(1.0f);
  }

  void sub(volume& dest, const volume& vol, 
	   uint64 off_x, uint64 off_y, uint64 off_z,
	   const dimension& subvoldim)
  {
    thread_info ti(BOOST_CURRENT_FUNCTION);

    if(!(dimension(off_x+subvoldim[0],off_y+subvoldim[1],off_z+subvoldim[2]) <= vol.voxel_dimensions()))
      throw index_out_of_bounds("Subvolume offset and dimension exceeds the boundary of input volume.");
    
    dest.unsetMinMax();
    dest.voxel_dimensions(subvoldim);
    dest.voxelType(vol.voxelType());
    dest.boundingBox(bounding_box(vol.XMin()+off_x*vol.XSpan(),
				  vol.YMin()+off_y*vol.YSpan(),
				  vol.ZMin()+off_z*vol.ZSpan(),
				  vol.XMin()+(off_x+subvoldim[0]-1)*vol.XSpan(),
				  vol.YMin()+(off_y+subvoldim[1]-1)*vol.YSpan(),
				  vol.ZMin()+(off_z+subvoldim[2]-1)*vol.ZSpan()));

    for(uint64 k=0; k<subvoldim[2]; k++)
      for(uint64 j=0; j<subvoldim[1]; j++)
	for(uint64 i=0; i<subvoldim[0]; i++)
	  dest(i,j,k, vol(i+off_x,j+off_y,k+off_z));
  }

  void volconvert(const std::string& input_volume_file,
                  const std::string& output_volume_file)
  {
    using namespace boost;

    thread_info ti(BOOST_CURRENT_FUNCTION);

    cvcapp.log(2,str(format("%s :: out-of-core convert\n")
                     % BOOST_CURRENT_FUNCTION));
      
    volume_file_info volinfo;
    volinfo.read(input_volume_file);

    //createVolumeFile in Utlity.h 
    createVolumeFile(output_volume_file,volinfo); 

    //read in slice by slice
    for(unsigned int k = 0; k < volinfo.ZDim(); k++)
      {
        for(unsigned int var=0; var<volinfo.numVariables(); var++)
          for(unsigned int time=0; time<volinfo.numTimesteps(); time++)
            {
              volume vol;
              readVolumeFile(vol,input_volume_file,
                             var,time,
                             0,0,k,
                             dimension(volinfo.XDim(),volinfo.YDim(),1));
                  
              vol.desc(volinfo.name(var));
              writeVolumeFile(vol,output_volume_file,
                              var,time,
                              0,0,k);
            }
        cvcapp.threadProgress(((float)k)/((float)((int)(volinfo.ZDim()-1))));
      }
  }
}
