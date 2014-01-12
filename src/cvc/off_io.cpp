#include <cvc/geometry_file_io.h>
#include <cvc/app.h>
#include <cvc/exception.h>
#include <cvc/utility.h>

#include <fstream>

namespace CVC_NAMESPACE
{
  CVC_DEF_EXCEPTION(invalid_off_file);

  struct off_io : public geometry_file_io
  {
    // --------------
    // off_io::off_io
    // --------------
    // Purpose:
    //   Initializes the extension list and id.
    // ---- Change History ----
    // 01/11/2014 -- Joe R. -- Creation.
    off_io()
      : _id("off_io : v1.0")
    {
      _extensions.push_back(".off");
    }

    // ----------
    // off_io::id
    // ----------
    // Purpose:
    //   Returns a string that identifies this geometry_file_io object.  This should
    //   be unique, but is freeform.
    // ---- Change History ----
    // 01/11/2014 -- Joe R. -- Creation.
    virtual const std::string& id() const
    {
      return _id;
    }

    // ------------------
    // off_io::extensions
    // ------------------
    // Purpose:
    //   Returns a list of extensions that this geometry_file_io object supports.
    // ---- Change History ----
    // 01/11/2014 -- Joe R. -- Creation.
    virtual const extension_list& extensions() const
    {
      return _extensions;
    }

    // ---------------
    // off_io::read
    // ---------------
    // Purpose:
    //   Reads a geomview OFF file.  This is not fully functional but it does handle
    //   most common variants of off files...
    // ---- Change History ----
    // 11/16/2011 -- arand -- Creation.
    // 01/11/2014 -- Joe R. -- Adapted for geometry_file_io.
    virtual geometry read(const std::string& filename) const
    {
      using namespace std;
      using namespace boost;

      geometry geom;

      std::ifstream inf(filename.c_str());
      if(!inf)
	throw read_error(string("Could not open ") + filename);

      unsigned int num_verts, num_elems;
      string line;
      vector<std::string> split_line;
      int line_num = 0;

      getline(inf, line); line_num++; // junk the first line
      if(!inf)
	throw read_error(str(format("Error reading file %1%, line %2%")
			     % filename
			     % line_num));

      //string headword;
      //inf >> headword;    
      //if (headword.compare("OFF") != 0 && 
      //  headword.compare("COFF") != 0) {
      //  throw read_error(str(format("Error reading header for file %1%")
      //                % filename));
      //}

      getline(inf, line); line_num++;
      if(!inf)
	throw read_error(str(format("Error reading file %1%, line %2%")
			     % filename
			     % line_num));
      trim(line);
      split(split_line,
	    line,
	    is_any_of(" "),
	    token_compress_on);
      if(split_line.size() != 3)
	throw read_error(str(format("Not an OFF file (wrong number of tokens in line 2: [%1%])")
			     % split_line.size()));

      try
	{
	  num_verts = lexical_cast<unsigned int>(split_line[0]);
	  num_elems = split_line.size() > 1 ?
	    lexical_cast<unsigned int>(split_line[1]) : 0;
	  
	  for(unsigned int vt = 0; vt < num_verts; vt++) 
	    {
	      getline(inf, line); line_num++;
	      if(!inf)
		throw read_error(str(format("Error reading file %1%, line %2%")
				     % filename
				     % line_num));
	      trim(line);
	      split(split_line,
		    line,
		    is_any_of(" "),
		    token_compress_on);       
	      
	      switch(split_line.size())
		{
		case 6: // colors
		case 7: // colors with alpha
		  // arand: ignoring transparency...
		  color_t color;
		  for(int i = 3; i < 6; i++)
		    color[i-3] = lexical_cast<double>(split_line[i]);
		  geom.colors().push_back(color);
		case 3: // no colors
		  point_t point;
		  for(int i = 0; i < 3; i++)
		    point[i] = lexical_cast<double>(split_line[i]);
		  geom.points().push_back(point);
		  break;
		default:
		  throw read_error(str(format("Error reading file %1%, line %2%")
				       % filename
				       % line_num));
		}
	    }
	}
      catch(std::exception& e)
	{
	  throw read_error(str(format("Error reading file %1%, line %2%, contents: '%3%', reason: %4%")
			       % filename
			       % line_num
			       % line
			       % string(e.what())));
	}
        
      for(unsigned int tri = 0; tri < num_elems; tri++)
	{
	  getline(inf, line); line_num++;
	  if(!inf)
	    throw read_error(str(format("Error reading file %1%, line %2%")
				 % filename
				 % line_num));
	  trim(line);
	  split(split_line,
		line,
		is_any_of(" "),
		token_compress_on);
	  
	  int size = lexical_cast<int>(split_line[0]);
	  
	  // don't handle segments yet... just planar facets
	  if (size < 3) 
	    {
	      throw read_error(str(format("Error reading file %1%, line %2%")
				   % filename
				   % line_num));
	    }

	  // arand: just triangulate every surface so that
	  //        we don't have special cases...
	  for (int i=2; i<split_line.size()-1; i++) 
	    {
	      tri_t tri;
	      tri[0] = lexical_cast<unsigned int>(split_line[1]);
	      tri[1] = lexical_cast<unsigned int>(split_line[i]);
	      tri[2] = lexical_cast<unsigned int>(split_line[i+1]);
	      geom.tris().push_back(tri);
	    }
	}

      return geom;
    }

    // -------------
    // off_io::write
    // -------------
    // Purpose:
    //   Supposed to write geometry to OFF file format, but is unimplemented so throws an error for now.
    // ---- Change History ----
    // 01/11/2014 -- Joe R. -- Creation.
    virtual void write(const geometry& geom, const std::string& filename) const
    {
      throw write_error(BOOST_CURRENT_FUNCTION);
    }

  protected:
    std::string _id;
    extension_list _extensions;
  };
}

namespace
{
  class off_io_init
  {
  public:
    off_io_init()
    {
      CVC_NAMESPACE::geometry_file_io::insert_handler(
        CVC_NAMESPACE::geometry_file_io::ptr(new CVC_NAMESPACE::off_io)
      );
    }
  } static_init;
}
