/*
  Copyright 2008-2011 The University of Texas at Austin

  Authors: Joe Rivera <transfix@ices.utexas.edu>
  Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolumeRover.

  VolumeRover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeRover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <cvc/geometry.h>
#include <cvc/utility.h>

#ifdef CVC_GEOMETRY_ENABLE_PROJECT
//Requires CGAL
#include <cvc/project_verts.h>
#endif

#ifdef CVC_GEOMETRY_ENABLE_BUNNY
#include <cvc/bunny.h>
#endif

#include <map>
#include <list>
#include <limits>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/regex.hpp>

//#define CVC_GEOMETRY_CORRECT_INDEX_START
//#define CVC_GEOMETRY_ENABLE_BUNNY

namespace CVC_NAMESPACE
{
  geometry::geometry()
    : _extents_set(false)
  {
    for(uint64_t i = 0; i < 3; i++)
      _min[i] = _max[i] = 0.0;

    init_ptrs();
  }

  geometry::geometry(const geometry& geom)
    : _points(geom._points), _boundary(geom._boundary),
      _normals(geom._normals), _colors(geom._colors),
      _lines(geom._lines), _tris(geom._tris),
      _quads(geom._quads), _extents_set(geom._extents_set),
      _min(geom._min), _max(geom._max) 
  {
    //make sure all our pointers are valid
    init_ptrs();
  }

  geometry::~geometry() {}

  void geometry::copy(const geometry &geom)
  {
    _points = geom._points;
    _boundary = geom._boundary;
    _normals = geom._normals;
    _colors = geom._colors;
    _lines = geom._lines;
    _tris = geom._tris;
    _quads = geom._quads;
    _extents_set = geom._extents_set;
    _min = geom._min;
    _max = geom._max;

    //make sure all our pointers are valid
    init_ptrs();
  }

  geometry& geometry::operator=(const geometry& geom)
  {
    copy(geom); return *this;
  }

  point_t geometry::min_point() const
  {
    if(!_extents_set) calc_extents();
    return _min;
  }

  point_t geometry::max_point() const
  {
    if(!_extents_set) calc_extents();
    return _max;
  }

  uint64_t geometry::num_points() const
  {
    //TODO: throw an exception if _points.size() != _boundary.size() 
    //!= _normals.size() != _colors.size()
    return _points ? _points->size() : 0;
  }

  uint64_t geometry::num_lines() const
  {
    return _lines ? _lines->size() : 0;
  }

  uint64_t geometry::num_tris() const
  {
    return _tris ? _tris->size() : 0;
  }

  uint64_t geometry::num_quads() const
  {
    return _quads ? _quads->size() : 0;
  }

  bool geometry::empty() const
  {
    if(num_points() == 0) // ||   (num_lines() == 0 && num_tris() == 0 && num_quads() == 0))
      return true;
    return false;
  }

  geometry& geometry::merge(const geometry& geom)
  {
    using namespace std;

    //append vertex info
    points().insert(points().end(),
                    geom.points().begin(),
                    geom.points().end());
    normals().insert(normals().end(),
                     geom.normals().begin(),
                     geom.normals().end());
    colors().insert(colors().end(),
                    geom.colors().begin(),
                    geom.colors().end());
        
    //we apparently don't have normal iterators for boundary_t :/
    unsigned int old_size = boundary().size();
    boundary().resize(boundary().size() + geom.boundary().size());
    for(unsigned int i = 0; i < geom.boundary().size(); i++)
      boundary()[old_size + i] = geom.boundary()[i];

    //append and modify index info
    lines().insert(lines().end(),
                   geom.lines().begin(),
                   geom.lines().end());
    {
      lines_t::iterator i = lines().begin();
      advance(i, lines().size() - geom.lines().size());
      while(i != lines().end())
        {
          for(line_t::iterator j = i->begin();
              j != i->end();
              j++)
            *j += points().size() - geom.points().size();
          i++;
        }
    }
    
    tris().insert(tris().end(),
                  geom.tris().begin(),
                  geom.tris().end());
    {
      tris_t::iterator i = tris().begin();
      advance(i, tris().size() - geom.tris().size());
      while(i != tris().end())
        {
          for(tri_t::iterator j = i->begin();
              j != i->end();
              j++)
            *j += points().size() - geom.points().size();
          i++;
        }
    }

    quads().insert(quads().end(),
                   geom.quads().begin(),
                   geom.quads().end());
    {
      quads_t::iterator i = quads().begin();
      advance(i, quads().size() - geom.quads().size());
      while(i != quads().end())
        {
          for(quad_t::iterator j = i->begin();
              j != i->end();
              j++)
            *j += points().size() - geom.points().size();
          i++;
        }
    }

    return *this;
  }

  geometry geometry::tri_surface() const
  {
    geometry new_geom;

    new_geom._points = _points;
    new_geom._colors = _colors;

    //if we don't have boundary info, just insert the tris directly
    if(boundary().size() != points().size())
      {
        new_geom._tris = _tris;

        for(quads_t::const_iterator i = quads().begin();
            i != quads().end();
            i++)
          {
            tri_t quad_tris[2];
            quad_tris[0][0] = (*i)[0];
            quad_tris[0][1] = (*i)[1];
            quad_tris[0][2] = (*i)[3];

            quad_tris[1][0] = (*i)[1];
            quad_tris[1][1] = (*i)[2];
            quad_tris[1][2] = (*i)[3];

            new_geom.tris().push_back(quad_tris[0]);
            new_geom.tris().push_back(quad_tris[1]);
          }
      }
    else
      {
        for(tris_t::const_iterator i = tris().begin();
            i != tris().end();
            i++)
          if(boundary()[(*i)[0]] &&
             boundary()[(*i)[1]] &&
             boundary()[(*i)[2]])
            new_geom.tris().push_back(*i);

        for(quads_t::const_iterator i = quads().begin();
            i != quads().end();
            i++)
          if(boundary()[(*i)[0]] &&
             boundary()[(*i)[1]] &&
             boundary()[(*i)[2]] &&
             boundary()[(*i)[3]])
            {
              tri_t quad_tris[2];
              quad_tris[0][0] = (*i)[0];
              quad_tris[0][1] = (*i)[1];
              quad_tris[0][2] = (*i)[3];
                
              quad_tris[1][0] = (*i)[1];
              quad_tris[1][1] = (*i)[2];
              quad_tris[1][2] = (*i)[3];
                
              new_geom.tris().push_back(quad_tris[0]);
              new_geom.tris().push_back(quad_tris[1]);
            }
      }
      
    return new_geom;
  }

  geometry& geometry::calculate_surf_normals()
  {
    using namespace std;
    geometry tri_geom = this->tri_surface();
      
    //find all the triangles that each point is a part of
    map<point_t, list<unsigned int> > neighbor_tris;
    for(tris_t::const_iterator i = tri_geom.const_tris().begin();
        i != tri_geom.const_tris().end();
        i++)
      for(tri_t::const_iterator j = i->begin();
          j != i->end();
          j++)
        neighbor_tris[tri_geom.const_points()[*j]]
          .push_back(distance(tri_geom.const_tris().begin(),i));

    //now iterate over our points, and for each point push back a new normal vector
    normals().clear();
    for(points_t::const_iterator i = const_points().begin();
        i != const_points().end();
        i++)
      {
        vector_t norm = {{ 0.0, 0.0, 0.0 }};
        list<unsigned int> &neighbors = neighbor_tris[*i];
        for(list<unsigned int>::iterator j = neighbors.begin();
            j != neighbors.end();
            j++)
          {
            vector_t cur_tri_norm;
            const tri_t &tri = tri_geom.const_tris()[*j];
            point_t p[3] = { tri_geom.const_points()[tri[0]],
                             tri_geom.const_points()[tri[1]],
                             tri_geom.const_points()[tri[2]] };
            vector_t v1 = {{ p[1][0] - p[0][0],
                             p[1][1] - p[0][1],
                             p[1][2] - p[0][2] }};
            vector_t v2 = {{ p[2][0] - p[0][0], 
                             p[2][1] - p[0][1],
                             p[2][2] - p[0][2] }};
            cross(cur_tri_norm,v1,v2);
            normalize(cur_tri_norm);
            for(int k = 0; k < 3; k++)
              norm[k] += cur_tri_norm[k];
          }
        if(!neighbors.empty())
          for(int j = 0; j < 3; j++)
            norm[j] /= neighbors.size();
        normals().push_back(norm);
      }

    return *this;
  }
  
  geometry geometry::generate_wire_interior() const
  {
    geometry new_geom(*this);

    if(!boundary().empty() && 
       boundary().size() == points().size())
      {
        //tet mesh
        if(!tris().empty())
          {
            new_geom.tris().clear();
            for(tris_t::const_iterator j = tris().begin();
                j != tris().end();
                j++)
              {
                //add the tri face boundary lines if both verts
                //of the potential line to be added do not lie 
                //on the boundary
                if(!boundary()[(*j)[0]] || !boundary()[(*j)[1]])
                  {
                    line_t line = {{ (*j)[0], (*j)[1] }};
                    new_geom.lines().push_back(line);
                  }
                if(!boundary()[(*j)[1]] || !boundary()[(*j)[2]])
                  {
                    line_t line = {{ (*j)[1], (*j)[2] }};
                    new_geom.lines().push_back(line);
                  }
                if(!boundary()[(*j)[2]] || !boundary()[(*j)[0]])
                  {
                    line_t line = {{ (*j)[2], (*j)[0] }};
                    new_geom.lines().push_back(line);
                  }
                  
                //add the triangle to new_geom if all verts are on the boundary
                if(boundary()[(*j)[0]] && boundary()[(*j)[1]] && boundary()[(*j)[2]])
                  new_geom.tris().push_back(*j);
              }
          }
        //hex mesh
        if(!quads().empty())
          {
            new_geom.quads().clear();
            for(quads_t::const_iterator j = quads().begin();
                j != quads().end();
                j++)
              {
                //same as comment above except for quads...
                if(!boundary()[(*j)[0]] || !boundary()[(*j)[1]])
                  {
                    line_t line = {{ (*j)[0], (*j)[1] }};
                    new_geom.lines().push_back(line);
                  }
                if(!boundary()[(*j)[1]] || !boundary()[(*j)[2]])
                  {
                    line_t line = {{ (*j)[1], (*j)[2] }};
                    new_geom.lines().push_back(line);
                  }
                if(!boundary()[(*j)[2]] || !boundary()[(*j)[3]])
                  {
                    line_t line = {{ (*j)[2], (*j)[3] }};
                    new_geom.lines().push_back(line);
                  }
                if(!boundary()[(*j)[3]] || !boundary()[(*j)[0]])
                  {
                    line_t line = {{ (*j)[3], (*j)[0] }};
                    new_geom.lines().push_back(line);
                  }
                  
                if(boundary()[(*j)[0]] && boundary()[(*j)[1]] && 
                   boundary()[(*j)[2]] && boundary()[(*j)[3]])
                  new_geom.quads().push_back(*j);
              }
          }
      }

    return new_geom;
  }
  
  geometry& geometry::invert_normals()
  {
    for(normals_t::iterator i = normals().begin();
        i != normals().end();
        i++)
      for(int j = 0; j < 3; j++)
        (*i)[j] *= -1.0;
    return *this;
  }

  geometry& geometry::reorient()
  {
    using namespace std;

    //TODO: make this work for quads

    if(!tris().empty())
      {
        //find all the triangles that each point is a part of
        map<point_t, list<unsigned int> > neighbor_tris;
        for(tris_t::iterator i = tris().begin();
            i != tris().end();
            i++)
          for(tri_t::iterator j = i->begin();
              j != i->end();
              j++)
            neighbor_tris[points()[*j]]
              .push_back(distance(tris().begin(),i));
          
        /*
          IGNORE THIS
            
          for each triangle (i),
          for each point (j) in that triangle (i)
          for each triangle (k) that point (j) is a part of
          if the angle between the normal of k and i is obtuse,
          reverse the vertex order of k
        */
          
        //if we don't have normals, lets calculate them
        if(points().size() != normals().size())
          calculate_surf_normals();
          
        for(points_t::iterator i = points().begin();
            i != points().end();
            i++)
          {
            //skip non boundary verts because their normals are null
            if(!boundary()[distance(points().begin(),i)]) continue;
              
            list<unsigned int> &neighbors = neighbor_tris[*i];
            for(list<unsigned int>::iterator j = neighbors.begin();
                j != neighbors.end();
                j++)
              {
                tri_t &tri = tris()[*j];
                for(int k=0; k<3; k++)
                  if(dot(normals()[tri[k]],
                         normals()[distance(points().begin(),i)]) < 0)
                    for(int l=0; l<3; l++)
                      normals()[tri[k]][l] *= -1.0;
              }
          }
      }
      
    return *this;
  }

  geometry& geometry::clear()
  {
    return (*this = geometry());
  }

  geometry& geometry::project(const geometry& input)
  {
#ifdef CVC_GEOMETRY_ENABLE_PROJECT
    geometry ref = input.tri_surface();
    project_verts::project(points().begin(),
                           points().end(),
                           ref.const_points().begin(),
                           ref.const_points().end(),
                           ref.const_tris().begin(),
                           ref.const_tris().end());
#else
# warning geometry::project() disabled
#endif
    return *this;
  }

  void geometry::init_ptrs()
  {
    if(!_points)
      _points.reset(new points_t);
    if(!_boundary)
      _boundary.reset(new boundary_t);
    if(!_normals)
      _normals.reset(new normals_t);
    if(!_colors)
      _colors.reset(new colors_t);
    if(!_lines)
      _lines.reset(new lines_t);
    if(!_tris)
      _tris.reset(new tris_t);
    if(!_quads)
      _quads.reset(new quads_t);
  }

  void geometry::calc_extents() const
  {
    using namespace std;
    
    if(empty())
      {
        fill(_min.begin(),_min.end(),numeric_limits<scalar_t>::max());
        fill(_max.begin(),_max.end(),-numeric_limits<scalar_t>::max());
        return;
      }
    
    _min = _max = points()[0];

    for(points_t::const_iterator i = points().begin();
        i != points().end();
        i++)
      {
        if(_min[0] > (*i)[0]) _min[0] = (*i)[0];
        if(_min[1] > (*i)[1]) _min[1] = (*i)[1];
        if(_min[2] > (*i)[2]) _min[2] = (*i)[2];
        if(_max[0] < (*i)[0]) _max[0] = (*i)[0];
        if(_max[1] < (*i)[1]) _max[1] = (*i)[1];
        if(_max[2] < (*i)[2]) _max[2] = (*i)[2];
      }

    _extents_set = true;

    return;
  }

  template<class container_ptr_t>
  void make_unique(container_ptr_t& cp)
  {
    if(!cp.unique())
      {
        container_ptr_t tmp(cp);
        cp.reset(new typename container_ptr_t::element_type(*tmp));
      }
  }

  void geometry::pre_write(ARRAY_TYPE at)
  {
    //invalidate calculated extents
    _extents_set = false;

    switch(at)
      {
      case POINTS: make_unique(_points); break;
      case BOUNDARY: make_unique(_boundary); break;
      case NORMALS: make_unique(_normals); break;
      case COLORS: make_unique(_colors); break;
      case LINES: make_unique(_lines); break;
      case TRIS: make_unique(_tris); break;
      case QUADS: make_unique(_quads); break;
      }
  }

  geometry& geometry::read(const std::string& filename)
  {
    if(filename == ".CVC_BUNNY.raw")
      {
	*this = bunny();
	return *this;
      }

    std::string errors;
    boost::regex file_extension("^(.*)(\\.\\S*)$");
    boost::smatch what;

    if(boost::regex_match(filename, what, file_extension))
      {
        if(what[2].compare(".raw") == 0 ||
           what[2].compare(".rawn") == 0 ||
           what[2].compare(".rawnc") == 0 ||
           what[2].compare(".rawc") == 0 )
          read_raw(filename);  // raw data types
        else if(what[2].compare(".off") == 0) 
          read_off(filename);  // off files
        else 
          throw CVC_NAMESPACE::unsupported_geometry_file_type(std::string(BOOST_CURRENT_FUNCTION) + 
                                                              std::string(": Cannot read ") + filename);
      }

    return *this;
  }

  // 12/27/2013 -- Joe R. -- Adapted from old io.h header in VolumeRover
  void geometry::write(const std::string& filename) const
  {
    using namespace std;
    using namespace boost;
        
    ofstream outf(filename.c_str());
    if(!outf)
      throw write_error(string("Could not open ") + filename);
    
    unsigned int num_elems;
    if(boundary().size() != points().size()) //if we don't have boundary information
      {
        num_elems = 
          lines().size() > 0 ? lines().size() :
          tris().size() > 0 ? tris().size()   :
          quads().size() > 0 ? quads().size() :
          0;
      }
    else
      {
        num_elems =
          tris().size() > 0 ? tris().size()/4   :
          quads().size() > 0 ? quads().size()/6 :
          0;
      }
    outf << points().size() << " " << num_elems << endl;
    if(!outf)
      throw write_error("Could not write number of points or number of tris to file!");

    // arand: 4-21-2011
    // changed to check file extension
    // normals and colors are only printed if .rawn/.rawnc/.rawc
    // have been specified...
    bool haveNormals = (normals().size() == points().size());
    bool printNormals = (ends_with(filename,".rawn") || ends_with(filename,".rawnc"));
    bool haveColors = (colors().size() == points().size());
    bool printColors = (ends_with(filename,".rawc") || ends_with(filename,".rawnc"));

    if (printNormals && !haveNormals) {
      std::cout << "WARNING: file with normals requested but not available." << std::endl;
    }

    if (printColors && !haveColors) {
      std::cout << "WARNING: file with normals requested but not available." << std::endl;
    }
    
    for(typename points_t::const_iterator i = points().begin();
        i != points().end();
        i++)
      {
        outf << (*i)[0] << " " << (*i)[1] << " " << (*i)[2];

        if(haveNormals && printNormals)
          outf << " " 
               << normals()[distance(points().begin(),i)][0] << " " 
               << normals()[distance(points().begin(),i)][1] << " " 
               << normals()[distance(points().begin(),i)][2];
        if(haveColors && printColors)
          outf << " " 
               << colors()[distance(points().begin(),i)][0] << " " 
               << colors()[distance(points().begin(),i)][1] << " " 
               << colors()[distance(points().begin(),i)][2];
        if(boundary().size() == points().size())
          outf << " " << boundary()[distance(points().begin(),i)];
        outf << endl;
        if(!outf)
          throw write_error(str(format("Error writing vertex %1%") % distance(points().begin(),i)));
      }
    
    if(lines().size() != 0)
      {
        for(typename lines_t::const_iterator i = lines().begin(); i != lines().end(); i++)
          {
            typedef typename lines_t::value_type cell_type;
            for(typename cell_type::const_iterator j = i->begin(); j != i->end(); j++)
              {
                outf << *j;
                if(next(j) == i->end()) outf << endl;
                else outf << " ";
              }

            if(!outf)
              throw write_error(str(format("Error writing line %1%") % distance(lines().begin(),i)));
          }
      }
    else if(tris().size() != 0)
      {
        if(boundary().size() != points().size()) //if we don't have boundary info, don't treat these tris as part of a tet
          {
            for(typename tris_t::const_iterator i = tris().begin(); i != tris().end(); i++)
              {
                typedef typename tris_t::value_type cell_type;
                for(typename cell_type::const_iterator j = i->begin(); j != i->end(); j++)
                  {
                    outf << *j;
                    if(next(j) == i->end()) outf << endl;
                    else outf << " ";
                  }
                
                if(!outf)
                  throw write_error(str(format("Error writing triangle %1%") % distance(tris().begin(),i)));
              }
          }
        else
          {
            for(unsigned int i = 0; i < tris().size()/4; i++)
              {
                outf << tris()[4*i][0] << " " << tris()[4*i][1] << " " 
                     << tris()[4*i][2] << " " << tris()[4*i+1][2] << endl;
                if(!outf)
                  throw write_error(str(format("Error writing tetrahedron %1%") % i));
              }
          }
      }
    else if(quads().size() != 0)
      {
        if(boundary().size() != points().size()) //if we don't have boundary info, don't tread these quads as part of a hexa
          {
            for(typename quads_t::const_iterator i = quads().begin(); i != quads().end(); i++)
              {
                typedef typename quads_t::value_type cell_type;
                for(typename cell_type::const_iterator j = i->begin(); j != i->end(); j++)
                  {
                    outf << *j;
                    if(next(j) == i->end()) outf << endl;
                    else outf << " ";
                  }
                
                if(!outf)
                  throw write_error(str(format("Error writing quad %1%") % distance(quads().begin(),i)));
              }
          }
        else
          {
            for(unsigned int i = 0; i < quads().size()/6; i++)
              {
#if 0
                outf << quads()[6*i][0] << " " << quads()[6*i][1] << " "
                     << quads()[6*i][2] << " " << quads()[6*i][3] << " "
                     << quads()[6*i+1][1] << " " << quads()[6*i+1][0] << " "
                     << quads()[6*i+1][3] << " " << quads()[6*i+1][2] << endl;
#endif
                outf << quads()[6*i][0] << " " << quads()[6*i][1] << " "
                     << quads()[6*i][2] << " " << quads()[6*i][3] << " "
                     << quads()[6*i+1][0] << " " << quads()[6*i+1][1] << " "
                     << quads()[6*i+1][2] << " " << quads()[6*i+1][3] << endl;
                

                if(!outf)
                  throw write_error(str(format("Error writing hexahedron %1%") % i));           
              }
          }
      }
  }

  // arand: written 4-11-2011
  //        directly read a cvc-raw type file into the data structure
  geometry::geometry(const std::string& filename)
  {
    init_ptrs();
    read(filename);
  }

  // arand, 11-16-2011: added off reader
  //        this is not fully functional but it does handle the most
  //        common variants of off files...
  void geometry::read_off(const std::string & filename)
  {
    using namespace std;
    using namespace boost;

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

    try {
      num_verts = lexical_cast<unsigned int>(split_line[0]);
      num_elems = split_line.size() > 1 ?
        lexical_cast<unsigned int>(split_line[1]) : 0;
      
      for(unsigned int vt = 0; vt < num_verts; vt++) {
        
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

        switch(split_line.size()) {

        case 6: // colors
        case 7: // colors with alpha
          // arand: ignoring transparency...
          color_t color;
          for(int i = 3; i < 6; i++)
            color[i-3] = lexical_cast<double>(split_line[i]);
          colors().push_back(color);
        case 3: // no colors
          point_t point;
          for(int i = 0; i < 3; i++) {
            point[i] = lexical_cast<double>(split_line[i]);
          }
          points().push_back(point);

          break;
        default:
          throw read_error(str(format("Error reading file %1%, line %2%")
                               % filename
                               % line_num));
        }
      }

    }catch(std::exception& e) {
      throw read_error(str(format("Error reading file %1%, line %2%, contents: '%3%', reason: %4%")
                           % filename
                           % line_num
                           % line
                           % string(e.what())));
    }
        

    for(unsigned int tri = 0; tri < num_elems; tri++) {
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
      if (size < 3) {
        throw read_error(str(format("Error reading file %1%, line %2%")
                             % filename
                             % line_num));
      }

      // arand: just triangulate every surface so that
      //        we don't have special cases...
      for (int i=2; i<split_line.size()-1; i++) {
        tri_t tri;
        tri[0] = lexical_cast<unsigned int>(split_line[1]);
        tri[1] = lexical_cast<unsigned int>(split_line[i]);
        tri[2] = lexical_cast<unsigned int>(split_line[i+1]);
        tris().push_back(tri);
      }
    }
  }

  // 04/02/2010 -- creation.
  void geometry::read_raw(const std::string & filename)
  {
    using namespace std;
    using namespace boost;    

    ifstream inf(filename.c_str());
    if(!inf)
      throw read_error(string("Could not open ") + filename);
    
    unsigned int line_num = 0;
    string line;
    vector<string> split_line;
    unsigned int num_verts, num_elems;
    
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
    if(split_line.size() != 1 &&
       split_line.size() != 2)
      throw read_error(str(format("Not a cvc-raw file (wrong number of tokens: [%1%])")
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
              case 3: // raw
              case 4: // raw with boundary
                {
                  point_t point;
                  for(int i = 0; i < 3; i++)
                    point[i] = lexical_cast<double>(split_line[i]);
                  points().push_back(point);
                  if(split_line.size() == 4)
                    boundary().push_back(lexical_cast<int>(split_line[3]));
                }
                break;
              case 6: // rawn or rawc
              case 7: // rawn or rawc with boundary
                {
                  bool is_rawn = ends_with(filename,".rawn");
                  point_t point;
                  for(int i = 0; i < 3; i++)
                    point[i] = lexical_cast<double>(split_line[i]);
                  points().push_back(point);

                  if(is_rawn)
                    {
                      vector_t normal;
                      for(int i = 3; i < 6; i++)
                        normal[i-3] = lexical_cast<double>(split_line[i]);
                      normals().push_back(normal);
                    }
                  else
                    {
                      color_t color;
                      for(int i = 3; i < 6; i++)
                        color[i-3] = lexical_cast<double>(split_line[i]);
                      colors().push_back(color);
                    }

                  if(split_line.size() == 7)
                    boundary().push_back(lexical_cast<int>(split_line[6]));
                }
                break;
              case 9:  // rawnc
              case 10: // rawnc with boundary
                {
                  point_t point;
                  vector_t normal;
                  color_t color;

                  for(int i = 0; i < 3; i++)
                    point[i-0] = lexical_cast<double>(split_line[i]);
                  points().push_back(point);

                  for(int i = 3; i < 6; i++)
                    normal[i-3] = lexical_cast<double>(split_line[i]);
                  normals().push_back(normal);

                  for(int i = 6; i < 9; i++)
                    color[i-6] = lexical_cast<double>(split_line[i]);
                  colors().push_back(color);

                  if(split_line.size() == 10)
                    boundary().push_back(lexical_cast<int>(split_line[9]));
                }
                break;
              default:
                {
                  throw read_error(str(format("Not a cvc-raw file (wrong number of tokens: [%1%])")
                                       % split_line.size()));
                }
                break;
              }
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
            switch(split_line.size())
              {
              case 2: // lines
                {
                  line_t line;
                  for(unsigned int i = 0; i < split_line.size(); i++)
                    line[i] = lexical_cast<unsigned int>(split_line[i]);
                  lines().push_back(line);
                }
                break;
              case 3: // tris
                {
                  tri_t tri;
                  for(unsigned int i = 0; i < split_line.size(); i++)
                    tri[i] = lexical_cast<unsigned int>(split_line[i]);
                  tris().push_back(tri);
                }
                break;
              case 4: // quads or tetrahedrons
                {
                  //if we didn't collect boundary information earlier, then we have quads
                  if(boundary().size() != points().size())
                    {
                      quad_t quad;
                      for(unsigned int i = 0; i < split_line.size(); i++)
                        quad[i] = lexical_cast<unsigned int>(split_line[i]);
                      quads().push_back(quad);
                    }
                  else
                    {
                      unsigned int t[4];
                      for(unsigned int i = 0; i < 4; i++)
                        t[i] = lexical_cast<unsigned int>(split_line[i]);
                      
                      tri_t tet_tris[4];
                      tet_tris[0][0] = t[0]; tet_tris[0][1] = t[2]; tet_tris[0][2] = t[1];
                      tet_tris[1][0] = t[1]; tet_tris[1][1] = t[2]; tet_tris[1][2] = t[3];
                      tet_tris[2][0] = t[0]; tet_tris[2][1] = t[3]; tet_tris[2][2] = t[2];
                      tet_tris[3][0] = t[0]; tet_tris[3][1] = t[1]; tet_tris[3][2] = t[3];

                      for(unsigned int i = 0; i < 4; i++)
                        tris().push_back(tet_tris[i]);
                    }
                }
                break;
              case 8: // hexahedrons
                {
                  //if we don't have boundary information, then something is wrong!
                  if(boundary().size() != points().size())
                    throw read_error("Incorrect cvc-raw file: missing boundary info for hex verts");

                  unsigned int t[8];
                  for(unsigned int i = 0; i < 8; i++)
                    t[i] = lexical_cast<unsigned int>(split_line[i]);
                  // TODO: investigate hex and tet index ordering... somehow I had to do this different than geoframe
                  quad_t hex_quads[6];
#if 0
                  hex_quads[0][0] = t[0]; hex_quads[0][1] = t[3]; hex_quads[0][2] = t[2]; hex_quads[0][3] = t[1];
                  hex_quads[1][0] = t[4]; hex_quads[1][1] = t[5]; hex_quads[1][2] = t[6]; hex_quads[1][3] = t[7];
                  hex_quads[2][0] = t[0]; hex_quads[2][1] = t[4]; hex_quads[2][2] = t[7]; hex_quads[2][3] = t[3];
                  hex_quads[3][0] = t[1]; hex_quads[3][1] = t[2]; hex_quads[3][2] = t[6]; hex_quads[3][3] = t[5];
                  hex_quads[4][0] = t[0]; hex_quads[4][1] = t[1]; hex_quads[4][2] = t[5]; hex_quads[4][3] = t[4];
                  hex_quads[5][0] = t[2]; hex_quads[5][1] = t[3]; hex_quads[5][2] = t[7]; hex_quads[5][3] = t[6];
#endif

                  hex_quads[0][0] = t[0]; hex_quads[0][1] = t[1]; hex_quads[0][2] = t[2]; hex_quads[0][3] = t[3];
                  hex_quads[1][0] = t[4]; hex_quads[1][1] = t[5]; hex_quads[1][2] = t[6]; hex_quads[1][3] = t[7];
                  hex_quads[2][0] = t[0]; hex_quads[2][1] = t[4]; hex_quads[2][2] = t[7]; hex_quads[2][3] = t[3];
                  hex_quads[3][0] = t[1]; hex_quads[3][1] = t[2]; hex_quads[3][2] = t[6]; hex_quads[3][3] = t[5];
                  hex_quads[4][0] = t[0]; hex_quads[4][1] = t[1]; hex_quads[4][2] = t[5]; hex_quads[4][3] = t[4];
                  hex_quads[5][0] = t[2]; hex_quads[5][1] = t[3]; hex_quads[5][2] = t[7]; hex_quads[5][3] = t[6];

                  for(unsigned int i = 0; i < 6; i++)
                    quads().push_back(hex_quads[i]);
                }
                break;
              default:
                {
                  throw read_error(str(format("Not a cvc-raw file {num components found: %1%}") % split_line.size()));   
                }
                break;
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

    //adjust indices if they start from 1 rather than 0.  (i.e. search for a 0 in each index
    //list. If there are no zeros, decrement each index by one)
#ifdef CVC_GEOMETRY_CORRECT_INDEX_START
    {
      bool lines_has_zero = false;
      bool tris_has_zero = false;
      bool quads_has_zero = false;

      for(typename lines_t::const_iterator i = lines().begin();
          i != lines().end();
          i++)
        if(find(i->begin(), i->end(), 0) != i->end())
          {
            lines_has_zero = true;
            break;
          }

      for(typename tris_t::const_iterator i = tris().begin();
          i != tris().end();
          i++)
        if(find(i->begin(), i->end(), 0) != i->end())
          {
            tris_has_zero = true;
            break;
          }
      
      for(typename quads_t::const_iterator i = quads().begin();
          i != quads().end();
          i++)
        if(find(i->begin(), i->end(), 0) != i->end())
          {
            quads_has_zero = true;
            break;
          }

      if(!lines_has_zero)
        for(typename lines_t::iterator i = lines().begin();
            i != lines().end();
            i++)
          for_each(i->begin(), i->end(), --_1);

      if(!tris_has_zero)
        for(typename tris_t::iterator i = tris().begin();
            i != tris().end();
            i++)
          for_each(i->begin(), i->end(), --_1);

      if(!quads_has_zero)
        for(typename quads_t::iterator i = quads().begin();
            i != quads().end();
            i++)
          for_each(i->begin(), i->end(), --_1);
    }
#endif
  }

  // 06/01/2012 -- initial implementation
  geometry bunny()
  {
    geometry geom;
#ifdef CVC_GEOMETRY_ENABLE_BUNNY
    for(size_t i = 0; i < 34835; i++)
      {
        point_t point;
        vector_t normal;
        for(int j = 0; j < 3; j++)
          point[j] = BUNNY_VERTS[i*3+j];
        for(int j = 0; j < 3; j++)
          normal[j] = BUNNY_NORMS[i*3+j];
        geom.points().push_back(point);
        geom.normals().push_back(normal);
      }

    for(size_t i = 0; i < 69473; i++)
      {
        tri_t tri;
        for(int j = 0; j < 3; j++)
          tri[j] = BUNNY_TRIS[i*3+j];
        geom.tris().push_back(tri);
      }
#endif
    return geom;
  }
}
