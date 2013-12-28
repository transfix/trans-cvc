/*
  Copyright 2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of libcvc.

  libcvc is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  libcvc is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id: App.h 5626 2012-05-23 22:07:11Z arand $ */

#ifndef __CVC_APP_H__
#define __CVC_APP_H__

#include <cvc/namespace.h>
#include <cvc/types.h>
#include <cvc/config.h>


#ifdef WIN32
#include <pthread.h> // arand, fix mingw problem
#endif

#include <boost/thread.hpp>
#include <boost/signals2.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/current_function.hpp>
#include <boost/format.hpp>
#include <boost/utility.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/property_tree/ptree.hpp>

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include <typeinfo>

namespace CVC_NAMESPACE
{
  // --------
  // cvc::app
  // --------
  // Purpose: 
  //   Main CVC Application class.  Used to keep track of all data objects
  //   used by this application, as well as the properties of the applicaton.
  //   The idea is that each CVC application with have their own application class
  //   that inherits this class.  The creator of the App subclasses should provide
  //   a static instance method to return a reference to the singleton object.
  //   Another possible design choice would be for CVC applications to simply
  //   use the existing App class and just use the data_map as a container for everything:
  //   UI objects such as the application's main window as well as regular data objects
  //   like Volumes and Geometry.
  // ---- Change History ----
  // 01/30/2011 -- Joe R. -- Initial implementation.
  // 02/06/2011 -- Joe R. -- Need to move the signal calls outside of the lock
  //                         scopes because it would cause deadlocks if slots
  //                         tried accessing the maps themselves.
  // 02/10/2011 -- Joe R. -- Moved most inline code to App.cpp.
  // 03/20/2011 -- Joe R. -- Added isData<> for convenience.
  // 04/06/2011 -- Joe R. -- Added data<T>() and log()
  // 04/18/2011 -- Joe R. -- Added propertyData<T>()
  // 04/23/2011 -- Joe R. -- Added thread progress
  // 04/29/2011 -- Joe R. -- Added listify and registerDataType
  // 06/17/2011 -- Joe R. -- Copied Exceptions, Types, and a few other things from VolMagick
  //                         Also added a new lookup table for mapping C++ types to DataType enums
  // 07/15/2011 -- Joe R. -- Added mutex_map stuff.  Used to easily define mutexes on files.
  // 07/29/2011 -- Joe R. -- Added another listProperty and added threadInfo stuff.
  // 09/09/2011 -- Joe R. -- Removing _instanceMutex static member because we're now using a mutex
  //                         wrapped in a helper class in app.cpp for the same usage.
  // 10/07/2011 -- Joe R. -- Added dataType(std::string)
  // 10/09/2011 -- Joe R. -- Added a new argument to startThread: 'wait'
  // 12/16/2011 -- Joe R. -- Added save/load property map functions.
  // 02/24/2012 -- Joe R. -- Moved wait_for_threads() to app.
  // 03/31/3012 -- Joe R. -- Added boost::any dataTypeName().
  class app
  {
  public:
    typedef boost::shared_ptr<app> app_ptr;

    //virtual ~app(); // arand: why virtual?
    ~app();

    // ***** Main API

    //Use instance() to grab a reference to the singleton application object.
    static app& instance();

    //Regular data access
    data_map data();
    boost::any data(const std::string& key);
    void data(const std::string& key, const boost::any& value);
    void data(const data_map& map);

    //Returns a nice name for the type of data referenced by the key
    //if it has been registered.
    std::string dataTypeName(const std::string& key);

    //Returns a nice name for the type of data if it has been registered
    template<class T>
    std::string dataTypeName()
    {
      boost::this_thread::interruption_point();
      boost::mutex::scoped_lock lock(_dataMutex);
      std::string rawname = typeid(T).name();
      if(_dataTypeNames.find(rawname)!=_dataTypeNames.end())
        return _dataTypeNames[rawname];
      else return rawname;
    }

    //Returns a nice name for the type of boost::any data if it has been registered
    std::string dataTypeName(const boost::any& d);

    //Returns the enum of the type of data with key.  If not found or
    //type not registered with enum, it will return Undefined.
    data_type dataType(const std::string& key);

    //Returns a data_type for the C++ type if one is available
    template<class T>
    data_type dataType()
    {
      boost::this_thread::interruption_point();
      boost::mutex::scoped_lock lock(_dataMutex);
      std::string rawname = typeid(T).name();
      if(_dataTypeEnum.find(rawname)!=_dataTypeEnum.end())
        return _dataTypeEnum[rawname];
      else return Undefined;
    }

    template<class T>
    T data(const std::string& key)
    {
      return boost::any_cast<T>(data(key));
    }

    template<class T>
    bool isData(const std::string& key)
    {
      try
        {
          T val = data<T>(key);
        }
      catch(std::exception& e)
        {
          return false;
        }
      return true;
    }

    //Return a vector of keys containing each datum of type T
    template<class T>
    std::vector<std::string> data()
    {
      std::vector<std::string> keys;
      data_map map = data();
      BOOST_FOREACH(data_map::value_type val, map)
        if(isData<T>(val.first))
          keys.push_back(val.first);
      return keys;
    }

    //Return a vector of objects of type T given an input container of strings
    template<class T>
    std::vector<T> data(const std::vector<std::string>& keys)
    {
      std::vector<T> ret;
      data_map map = data();
      BOOST_FOREACH(std::string key, keys)
        if(isData<T>(key))
          ret.push_back(data<T>(key));
      return ret;
    }

    //For each key, copy a corresponding datum in the vector to the datamap
    template<class Object>
    void data(const std::vector<std::string>& keys, const std::vector<Object>& v)
    {
      if(keys.empty() || v.empty()) return;
      for(size_t i = 0; 
          i < keys.size() && i < v.size();
          i++)
        data(keys[i],v[i]);
    }

    //Duplicate the value across all keys in the vector
    template<class T>
    void data(const std::vector<std::string>& keys, const T& value)
    {
      BOOST_FOREACH(std::string key, keys)
        data(key, value);
    }

    //Given a list of data keys as a comma separated list contained in a string,
    //returns a vector of objects of type T, each corresponding to a key in the list
    template<class T>
    std::vector<T> listData(const std::string& keylist)
    {
      using namespace std;
      using namespace boost;
      using namespace boost::algorithm;
      string separators = properties("system.list_separators");
      if(separators.empty()) separators=",";
      vector<string> keys;
      split(keys, keylist, is_any_of(separators));
      BOOST_FOREACH(string& key, keys) trim(key);
      return data<T>(keys);
    }

    std::vector<std::string> listify(const std::string& keylist);
    std::string listify(const std::vector<std::string>& keys);
    
    data_reader_collection dataReaders();
    void dataReaders(const data_reader_collection& dlc);
    data_reader dataReader(data_reader_collection::size_type idx);
    void dataReader(const data_reader& dl);

    bool readData(const std::string& path);

    property_map properties();
    std::string properties(const std::string& key);
    void properties(const std::string& key, const std::string& val);
    void properties(const property_map& map);
    void addProperties(const property_map& map);
    bool hasProperty(const std::string& key);

    //Given an input property key, output's the value for that key as
    //a vector of strings.  This is useful if the value of your property is a comma separated list
    std::vector<std::string> listProperty(const std::string& key, 
                                          bool uniqueElements = false);

    //Returns a collection of properties of type T.  Convenience function for above.
    // 07/29/2011 - Joe. R. -- initial implementation
    template<class T>
    std::vector<T> listProperty(const std::string& key,
                                bool uniqueElements = false)
    {
      using namespace std;
      using namespace boost;
      vector<string> vals =
        listProperty(key, uniqueElements);
      vector<T> ret_data;
      BOOST_FOREACH(string dkey, vals)
        {
          trim(dkey);
          ret_data.push_back(lexical_cast<T>(dkey));
        }
      return ret_data;
    }
    
    void listPropertyAppend(const std::string& key, const std::string& val);
    void listPropertyRemove(const std::string& key, const std::string& val);

    //shortcut for adding any type data to properties via
    //boost lexical_cast
    template<class T>
    void properties(const std::string& key, const T& val)
    {
      properties(key,boost::lexical_cast<std::string>(val));
    }

    //shortcut for retrieving any type via lexical_cast
    // 12/02/2011 -- transfix -- checking for property existence before attempting lexical_cast
    template<class T>
    T properties(const std::string& key)
    {
      if(hasProperty(key))
        return boost::lexical_cast<T>(properties(key));
      else
        return T();
    }

    //Shortcut for getting a vector of objects of type T from
    //the datamap using data keys stored in a property map value as a list.
    //Lists are just strings separated by any character in the value of the
    //system.list_separators property.
    template<class T>
    std::vector<T> propertyData(const std::string& propKey,
                                bool uniqueElements = false)
    {
      using namespace std;
      using namespace boost;
      using namespace boost::algorithm;
      vector<string> vals =
        listProperty(propKey, uniqueElements);
      vector<T> ret_data;
      BOOST_FOREACH(string dkey, vals)
        {
          trim(dkey);
          if(isData<T>(dkey))
            ret_data.push_back(data<T>(dkey));
        }
      return ret_data;
    }

    void readPropertyMap(const std::string& path);
    void writePropertyMap(const std::string& path);

    // ***** Thread API
    thread_map threads();
    thread_ptr threads(const std::string& key);
    void threads(const std::string& key, const thread_ptr& val);
    void threads(const thread_map& map);
    bool hasThread(const std::string& key);
    double threadProgress(const std::string& key = std::string());
    void threadProgress(double progress); //0.0 - 1.0
    void threadProgress(const std::string& key, double progress);
    void finishThreadProgress(const std::string& key = std::string());
    std::string threadKey(); //returns the thread key for this thread
    void removeThread(const std::string& key);
    std::string uniqueThreadKey(const std::string& hint = std::string());

    //set a string to associate with the thread to state it's current activity
    void threadInfo(const std::string& key, const std::string& infostr);
    std::string threadInfo(const std::string& key = std::string());
    void thisThreadInfo(const std::string& infostr)
    {
      threadInfo(std::string(),infostr);
    }
    std::string thisThreadInfo()
    {
      return threadInfo();
    }

    //T is a class with operator()
    // 10/09/2011 -- transfix -- added new argument 'wait'
    template<class T>
    void startThread(const std::string& key, const T& t, bool wait = true)
    {
      //If waiting and an existing thread with this key is running,
      //stop the existing running thread with this key and wait for
      //it to actually stop before starting a new one.  Else just use
      //a unique key.
      if(wait && hasThread(key))
        {
          thread_ptr t = threads(key);
          t->interrupt(); //initiate thread quit
          t->join();      //wait for it to quit
        }

      threads(
        wait ? key : app::instance().uniqueThreadKey(key),
        CVC_NAMESPACE::thread_ptr(new boost::thread(t))
      );
    }

    //Used to easily manage saving/restoring thread info as we
    //traverse a threads stack.
    class thread_info
    {
    public:
      thread_info(const std::string& info = "running")
        {
          _origInfo = app::instance().thisThreadInfo();
          _origProgress = app::instance().threadProgress();
          app::instance().thisThreadInfo(info);
        }
      ~thread_info()
        {
          app::instance().thisThreadInfo(_origInfo);
          app::instance().threadProgress(_origProgress);
        }
    private:
      std::string _origInfo;
      double _origProgress;
    };

    //Instantiate one of these on the shallowest point of your 
    //thread's stack
    class thread_feedback
    {
    public:
      thread_feedback(const std::string& info = "running")
        : _threadInfo(info)
        {
          app::instance().threadProgress(0.0);
        }

      ~thread_feedback()
        {
          app::instance().finishThreadProgress();
          app::instance().removeThread(app::instance().threadKey());
        }
    private:
      thread_info _threadInfo;
    };

    //output
    void log(unsigned int level, const std::string& buf);

    template<class T>
    void registerDataType(const std::string& datatypename)
    {
      boost::this_thread::interruption_point();
      boost::mutex::scoped_lock lock(_dataMutex);
      _dataTypeNames[typeid(T).name()] = datatypename;
    }

    template<class T>
    void registerDataType(data_type dt)
    {
      boost::this_thread::interruption_point();
      boost::mutex::scoped_lock lock(_dataMutex);
      _dataTypeEnum[typeid(T).name()] = dt;
    }

#ifndef registerDataType
#define registerDataType( type )  registerDataType<type>( #type )
#else
#warning registerDataType already defined!
#endif

    mutex_ptr mutex(const std::string& name);

    //info is an arbitrary string that can be used to describe
    //who currently has a lock on the mutex
    void mutexInfo(const std::string& name, const std::string& in);
    std::string mutexInfo(const std::string& name);

    // --------------------
    // scoped_lock
    // --------------------
    // Purpose:
    //   Class used to lock a mutex on the mutex_map.  The
    //   Specified mutex will be locked as long as this object
    //   is in scope.
    // ---- Change History ----
    // 07/15/2011 -- Joe R. -- Initial implementation.
    class scoped_lock
    {
    public:
      scoped_lock(const std::string& name,
                  const std::string& info = std::string()) : 
        _scopedLock(*app::instance().mutex(name)),
        _name(name)
          {
            //prepend the thread key
            app::instance().mutexInfo(name,
                                      app::instance().threadKey() + ": " + info);
          }
      ~scoped_lock()
        {
          app::instance().mutexInfo(_name,"");
        }
    private:
      boost::mutex::scoped_lock _scopedLock;
      std::string _name;
    };

    //Connect to these signals to monitor this application's state.
    map_change_signal dataChanged;
    map_change_signal propertiesChanged;
    map_change_signal threadsChanged;
    map_change_signal mutexesChanged;

    void wait() { wait_for_threads(); } //non-static for convenience

  protected:
    app();

    void propertyTreeTraverse(const boost::property_tree::ptree& pt,
                              const std::string& parentkey = std::string());

    data_map               _data;
    data_type_name_map     _dataTypeNames;
    data_type_enum_map     _dataTypeEnum; //mapping of C++ types to the data_type enums
    boost::mutex           _dataMutex;
    property_map           _properties;
    boost::mutex           _propertiesMutex;
    thread_map             _threads;
    thread_progress_map    _threadProgress;
    thread_key_map         _threadKeys;
    thread_info_map        _threadInfo;
    void                   updateThreadKeys();
    boost::mutex           _threadsMutex;
    boost::mutex           _logMutex;

    //Collection of functions that can read data from the filesystem (or URI??)
    //and stick it in the datamap.  A file loading operation will iterate
    //through the vector and call each function, passing the requested
    //filename.  The first one to return true will break the iteration,
    //signifying that it was successful in reading the data off disk and
    //sticking it in the map.
    data_reader_collection _dataReaders;
    boost::mutex           _dataReadersMutex;

    mutex_map              _mutexMap;
    boost::mutex           _mutexMapMutex;

    static app_ptr         instancePtr();
    static app_ptr         _instance;
    static boost::mutex    _instanceMutex;

    static void wait_for_threads();
  private:
    app(const app&);
  };

  typedef app::thread_info     thread_info;
  typedef app::thread_feedback thread_feedback;
  typedef app::scoped_lock     scoped_lock;
}

//Shorthand to access the app object from anywhere
#define cvcapp CVC_NAMESPACE::app::instance()

#endif
