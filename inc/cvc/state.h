/*
  Copyright 2012 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of libCVC.

  libCVC is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  libCVC is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id: State.h 5559 2012-05-11 21:43:22Z transfix $ */

#ifndef __CVC_STATE_H__
#define __CVC_STATE_H__

#include <cvc/namespace.h>
#include <cvc/types.h>
#include <cvc/exception.h>
#include <cvc/app.h>

#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>

namespace CVC_NAMESPACE
{
  CVC_DEF_EXCEPTION(network_error);
  CVC_DEF_EXCEPTION(xmlrpc_server_terminate);

  // ----------
  // cvc::state
  // ----------
  // Purpose: 
  //   Central program state manangement.  Provides a tree
  //   to which property values and arbitrary data can be attached.
  //   Written to be thread safe and to be used also as a thread
  //   messaging system.  With xmlrpc that messaging can extend to
  //   threads in other processes and nodes on the network.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.
  // 03/02/2012 -- Joe R. -- Added touch()
  // 03/15/2012 -- Joe R. -- Added initialized flag.
  // 03/16/2012 -- Joe R. -- Added reset(), ptree() and traverse()
  // 03/30/2012 -- Joe R. -- Added comment and hidden field.
  // 03/31/2012 -- Joe R. -- Added dataTypeName().
  class state
  {
  public:  
    typedef boost::shared_ptr<state> state_ptr;
    typedef std::map<std::string,state_ptr> child_map;
    typedef boost::function<void (std::string)> traversal_unary_func;

    static const std::string SEPARATOR;

    virtual ~state();

    // ***** Main API

    //Use instance() to grab a reference to the singleton application object.
    static state& instance();

    const std::string& name()   const { return _name;   }
    const state*       parent() const { return _parent; }

    //return's the parent's fullName
    std::string parentName() const { 
      std::string tmp;
      return parent() ? 
        (!parent()->name().empty() ?
         ((tmp=parent()->parentName()).empty()?
          parent()->name() : 
          tmp + SEPARATOR + parent()->name()) 
         : "")
        : "";
    }

    std::string fullName() const { 
      std::string pn = parentName();
      return pn.empty() ?
        name() :
        pn + SEPARATOR + name();
    }

    boost::posix_time::ptime lastMod();

    std::string value();
    std::string valueTypeName();
    std::vector<std::string> values(bool unique = false); //shortcut for comma separated values in value()
    state& value(const std::string& v, bool setValueType = true);
    template <class T> T value() { return boost::lexical_cast<T>(value()); }
    template <class T> state& value(const T& v) {
      {
        boost::mutex::scoped_lock lock(_mutex);
        _valueTypeName = cvcapp.dataTypeName<T>();
      }
      return value(boost::lexical_cast<std::string>(v),false);
    }
    signal valueChanged;

    boost::any data();
    state& data(const boost::any&);

    template<class T>
    T data()
    {
      return boost::any_cast<T>(data());
    }

    template<class T>
    bool isData()
    {
      try
	{
	  T val = data<T>();
	}
      catch(std::exception& e)
	{
	  return false;
	}
      return true;
    }

    std::string dataTypeName();

    signal dataChanged;

    state& operator()(const std::string& childname = std::string());
    std::vector<std::string> children(const std::string& re = std::string());
    size_t numChildren();
    map_change_signal childChanged;

    operator std::string(){ return value(); }
    
    signal destroyed;

    void touch();

    bool initialized() const { return _initialized; }

    //like propertyData from CVC::App
    template<class T>
    std::vector<T> valueData(bool uniqueElements = false)
    {
      using namespace std;
      using namespace boost;
      using namespace boost::algorithm;
      vector<string> vals = values(uniqueElements);
      vector<T> ret_data;
      BOOST_FOREACH(string dkey, vals)
        {
          trim(dkey);
          if(CVC_NAMESPACE::state::instance()(dkey).isData<T>())
            ret_data.push_back(CVC_NAMESPACE::state::instance()(dkey).data<T>());
        }
      return ret_data;
    }

    void reset();

    //converting to and from a boost property tree.  Useful for saving and restoring state.
    boost::property_tree::ptree ptree();
    operator boost::property_tree::ptree(){ return ptree(); }
    void ptree(const boost::property_tree::ptree&);

    void save(const std::string& filename);
    void restore(const std::string& filename);

    void traverse(traversal_unary_func func, const std::string& re = std::string());
    signal traverseEnter;
    signal traverseExit;

    std::string comment();
    state& comment(const std::string& c);
    signal commentChanged;

    bool hidden();
    state& hidden(bool h);
    signal hiddenChanged;

  protected:
    state(const std::string& n = std::string(),
          const state* p = NULL);

    void notifyParent(const std::string& childname);
    void notifyXmlRpc();

    boost::mutex                     _mutex;
    boost::posix_time::ptime         _lastMod;

    std::string                      _name;
    const state*                     _parent;

    std::string                      _value;
    std::string                      _valueTypeName;
    boost::any                       _data;
    std::string                      _comment;
    bool                             _hidden;
    child_map                        _children;

    bool                             _initialized;
    
    static state_ptr                 instancePtr();
    static state_ptr                 _instance;
    static boost::mutex              _instanceMutex;
  private:
    state(const state&);
  };
}

//Shorthand to access the cvc::state object from anywhere
#define cvcstate CVC_NAMESPACE::state::instance()

#endif
