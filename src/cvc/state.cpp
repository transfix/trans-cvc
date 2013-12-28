/*
  Copyright 2012 The University of Texas at Austin

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

#include <cvc/state.h>
#include <cvc/app.h>
#include <cvc/exception.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/current_function.hpp>
#include <boost/regex.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/asio.hpp>

#ifndef CVC_STATE_XML_PROPERTY_TREE
#include <boost/property_tree/info_parser.hpp>
#else
#include <boost/property_tree/xml_parser.hpp>
#endif

#ifdef USING_XMLRPC
#include <xmlrpc/XmlRpc.h>
#endif

#include <algorithm>
#include <set>

#ifdef USING_XMLRPC
namespace
{
  // -----------------
  // getLocalIpAddress
  // -----------------
  // Purpose: 
  //   This function is kind of a hack way to get the default interface's
  //   local ip address.  From http://bit.ly/ADIcC1
  // ---- Change History ----
  // 02/24/2012 -- Joe R. -- Initial implementation.
  std::string getLocalIPAddress()
  {
    using namespace boost::asio;

    ip::address addr;
    try
      {
        io_service netService;
        ip::udp::resolver   resolver(netService);
        ip::udp::resolver::query query(ip::udp::v4(), "cvcweb.ices.utexas.edu", "");
        ip::udp::resolver::iterator endpoints = resolver.resolve(query);
        ip::udp::endpoint ep = *endpoints;
        ip::udp::socket socket(netService);
        socket.connect(ep);
        addr = socket.local_endpoint().address();
      } 
    catch (std::exception& e)
      {
        throw CVC_NAMESPACE::network_error(e.what());
      }
    return addr.to_string();
  } 

  // --------------------
  // notify_xmlrpc_thread
  // --------------------
  // Purpose: 
  //   Sends the value of the state object specified in 'which'
  //   to the xmlrpc server running at host:port.
  // ---- Change History ----
  // 02/20/2012 -- Joe R. -- Initial implementation.
  // 03/09/2012 -- Joe R. -- Using stateName string instead of direct State ptr
  class notify_xmlrpc_thread
  {
  public:
    notify_xmlrpc_thread(const std::string& threadName,
                         const std::string& host, int port,
                         const std::string& stateName)
      : _threadName(threadName), _host(host), 
        _port(port), _stateName(stateName) {}
    
    notify_xmlrpc_thread(const notify_xmlrpc_thread& t)
      : _threadName(t._threadName), _host(t._host), 
        _port(t._port), _stateName(t._stateName) {}
    
    notify_xmlrpc_thread& operator=(const notify_xmlrpc_thread& t)
    {
      _threadName = t._threadName;
      _host = t._host;
      _port = t._port;
      _stateName = t._stateName;
      return *this;
    }
    
    void operator()()
    {
      using namespace boost;
      using namespace std;
      
      CVC_NAMESPACE::thread_feedback feedback;

      XmlRpc::XmlRpcClient c(_host.c_str(),_port);
      XmlRpc::XmlRpcValue params,result;

      params[0] = _stateName;
      params[1] = cvcstate(_stateName).value();
      params[2] = posix_time::to_simple_string(cvcstate(_stateName).lastMod());

      for(int i = 0; i < 3; i++)
        cvcapp.log(6,str(format("%s :: params[%d] = %s\n")
                         % BOOST_CURRENT_FUNCTION
                         % i
                         % string(params[i])));

      c.execute("cvcstate_set_value",params,result);
    }

    const std::string& threadName() const { return _threadName; }
    const std::string& stateName() const  { return _stateName; }

  protected:
    std::string _threadName;
    std::string _host;
    int _port;
    std::string _stateName;
  };

  // --------------------------
  // notify_xmlrpc_thread_setup
  // --------------------------
  // Purpose: 
  //   The initial thread spawned by a value change.  This will in turn
  //   spawn 1 thread for each host specified in __system.xmlrpc.hosts.
  // ---- Change History ----
  // 02/20/2012 -- Joe R. -- Initial implementation.
  // 03/09/2012 -- Joe R. -- Using stateName string instead of direct State ptr
  // 04/21/2012 -- Joe R. -- Instead of synching everyone, do it one way.
  class notify_xmlrpc_thread_setup
  {
  public:
    notify_xmlrpc_thread_setup(const std::string stateName) : _stateName(stateName) {}

    void operator()()
    {
      using namespace std;
      using namespace boost;
      using namespace boost::algorithm;

      CVC_NAMESPACE::thread_feedback feedback;
    
      //if no hosts have been set, don't do anything.
      if(cvcstate("__system.xmlrpc.hosts").value().empty())
        {
          cvcapp.log(3,str(format("%s :: no hosts listed in __system.xmlrpc.hosts\n")
                           % BOOST_CURRENT_FUNCTION));
          return;
        }

      //If this object is under any hierarchy specified in notify_states, forward over xmlrpc.
      {
        vector<string> vals = cvcstate("__system.xmlrpc.notify_states").values(true);

        vector<string> parts;
        string fn = _stateName;
        split(parts,fn,is_any_of(CVC_NAMESPACE::state::SEPARATOR));
        if(parts.size()<=1) return;
        bool filter = true;
        BOOST_FOREACH(string name, vals)
          {
            if(parts[0] == name)
              filter = false; //if it is one of the notify_states, call xmlrpc
          }
        if(filter) return;
      }

      vector<string> hosts = cvcstate("__system.xmlrpc.hosts").values(true);
      BOOST_FOREACH(string host, hosts)
        {
          vector<string> parts;
          split(parts,host,is_any_of(":"));
          if(parts.empty()) continue;
          string hostname = parts[0];
          int port = cvcstate("__system.xmlrpc.port").value<int>(); //use the port for this process
          if(parts.size()>1)
            port = lexical_cast<int>(parts[1]);
          string threadName = "notify_xmlrpc_thread_"+_stateName+"_"+host;

          cvcapp.log(6,str(format("%s :: hostname = %s, port = %d, name = %s\n")
                           % BOOST_CURRENT_FUNCTION
                           % hostname
                           % port
                           % _stateName));

          //Stick the thread on the datamap.  There is another thread that will actually launch
          //these threads some time later.  Doing this, we won't flood the network with xmlrpc
          //requests if several state changes happen in quick succession.
          cvcapp.data(threadName, notify_xmlrpc_thread(threadName,hostname,port,_stateName));
        }
    }
  protected:
    std::string _stateName;
  };

  // -----------------------------
  // process_notify_xmlrpc_threads
  // -----------------------------
  // Purpose: 
  //   Starts all the notify threads that are on the data map.
  // ---- Change History ----
  // 03/10/2012 -- Joe R. -- creation.
  class process_notify_xmlrpc_threads
  {
  public:
    void operator()()
    {
      while(1)
        {
          //Sleep for 200ms before each iteration.
          {
            CVC_NAMESPACE::thread_info ti("sleeping");
            boost::xtime xt;
            boost::xtime_get( &xt, boost::TIME_UTC );
            xt.nsec += 1000000000 / 5;
            boost::thread::sleep( xt );
          }
          
          std::vector<std::string> keys = cvcapp.data<notify_xmlrpc_thread>();
          std::vector<notify_xmlrpc_thread> threads = 
            cvcapp.data<notify_xmlrpc_thread>(keys);
          BOOST_FOREACH(notify_xmlrpc_thread thread, threads)
            {
              cvcapp.startThread(thread.threadName(),thread,false);
              cvcapp.data(thread.threadName(),boost::any()); //erase from the datamap
            }
        }
    }
  };

#define XMLRPC_METHOD_PROTOTYPE(name, description)                      \
  class name : public XmlRpc::XmlRpcServerMethod                        \
  {                                                                     \
  public:                                                               \
    name(XmlRpc::XmlRpcServer *s) :                                     \
      XmlRpc::XmlRpcServerMethod(#name, s) {}                           \
    void execute(XmlRpc::XmlRpcValue &params,                           \
                 XmlRpc::XmlRpcValue &result);                          \
    std::string help() { return std::string(description); }             \
  };

#define XMLRPC_METHOD_DEFINITION(name)                                    \
  void xmlrpc_server_thread::name::execute(XmlRpc::XmlRpcValue &params,   \
                                           XmlRpc::XmlRpcValue &result)

  // --------------------
  // xmlrpc_server_thread
  // --------------------
  // Purpose: 
  //   The thread that manages the XmlRpcServer instance.
  // ---- Change History ----
  // 02/20/2012 -- Joe R. -- Initial implementation.
  // 02/24/2012 -- Joe R. -- Moving default initilization here to avoid deadlock
  // 03/02/2012 -- Joe R. -- Running a thread to sync up with other hosts.
  // 03/10/2012 -- Joe R. -- Starting process_notify_xmlrpc_threads.
  class xmlrpc_server_thread
  {
  public:
    xmlrpc_server_thread() {}

    void operator()()
    {
      CVC_NAMESPACE::thread_feedback feedback;

      //instantiate the server and its methods.
      XmlRpc::XmlRpcServer s;
      cvcstate_set_value set_value(&s);
      cvcstate_get_value get_value(&s);
      cvcstate_get_state_names get_state_names(&s);
      cvcstate_terminate terminate(&s);

      if(cvcstate("__system.xmlrpc.port").value().empty())
        cvcstate("__system.xmlrpc.port")
          .value(int(23196))
          .comment("The port used by the xmlrpc server.");
      if(cvcstate("__system.xmlrpc.hosts").value().empty())
        cvcstate("__system.xmlrpc.hosts")
          .value("localhost:23196") //loopback test for now
          .comment("Comma separated list of host:port combinations used to broadcast "
                   "node changes via notify_xmlrpc_thread.");
      if(cvcstate("__system.xmlrpc.notify_states").value().empty())
        cvcstate("__system.xmlrpc.notify_states")
          .comment("Comma separated list of nodes to broadcast notification for.");

      int port = cvcstate("__system.xmlrpc.port").value<int>();
      std::string portstr = cvcstate("__system.xmlrpc.port");

      cvcapp.startThread("process_notify_xmlrpc_threads",process_notify_xmlrpc_threads(),false);

      try
        {
          std::string host = boost::asio::ip::host_name();
          std::string ipaddr = getLocalIPAddress();

          //Useful info to have
          cvcstate("__system.xmlrpc.hostname")
            .value(host)
            .comment("The hostname of the host running the xmlrpc server thread.");
          cvcstate("__system.xmlrpc.ipaddr")
            .value(ipaddr)
            .comment("The ip address bound by the xmlrpc server.");
 
          //Start the server, and run it indefinitely.
          //For some reason, time_from_string and boost_regex creashes if the main thread is waiting in atexit().
          //So, make sure main() has a cvcapp.wait_for_threads() call at the end.
          XmlRpc::setVerbosity(0);
          s.bindAndListen(port);
          s.enableIntrospection(true);
          s.work(-1.0);
        }
      catch(std::exception& e)
        {
          using namespace boost;
          cvcapp.log(1,str(format("%s :: xmlrpc_server_thread shutting down: %s\n")
                           % BOOST_CURRENT_FUNCTION % e.what()));
        }
    }

  private:
    //our exported methods
    XMLRPC_METHOD_PROTOTYPE(cvcstate_set_value, "Sets a state object's value");
    XMLRPC_METHOD_PROTOTYPE(cvcstate_get_value, "Gets a state object's value");
    XMLRPC_METHOD_PROTOTYPE(cvcstate_get_state_names, "Get a list of root's children using a PERL regular expression");
    XMLRPC_METHOD_PROTOTYPE(cvcstate_terminate, "Quits the server");
  };

  XMLRPC_METHOD_DEFINITION(cvcstate_set_value)
  {
    using namespace std;
    using namespace boost;
    using namespace boost::posix_time;

    string fullStateName = params[0];
    string stateval = params[1];
    ptime modtime = time_from_string(params[2]);

    for(int i = 0; i < 3; i++)
      cvcapp.log(6,str(format("%s :: params[%d] = %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % i
                       % string(params[i])));

    //search for a child with this state name
    vector<string> children = cvcstate().children(fullStateName);

    //if the object doesn't exist, or if the incoming value is newer, set it.
    if(children.empty() ||
       modtime > cvcstate(fullStateName).lastMod())
      {
        cvcstate(fullStateName).value(stateval);

        std::vector<std::string> children = cvcstate().children();
        BOOST_FOREACH(std::string child, children)
          cvcapp.log(4,str(format("%s :: %s = %s\n")
                           % BOOST_CURRENT_FUNCTION
                           % child
                           % cvcstate(child).value()));
      }
  }

  XMLRPC_METHOD_DEFINITION(cvcstate_get_value)
  {
    using namespace std;
    using namespace boost;
    using namespace boost::posix_time;

    string fullStateName = params[0];
    result[0] = cvcstate(fullStateName).value();
    result[1] = to_simple_string(cvcstate(fullStateName).lastMod());

    cvcapp.log(6,str(format("%s :: fullStateName = %s\n")
                     % BOOST_CURRENT_FUNCTION
                     % fullStateName));
  }

  XMLRPC_METHOD_DEFINITION(cvcstate_get_state_names)
  {
    using namespace std;
    using namespace boost;

    vector<string> ret = cvcstate().children(params[0]);
    for(size_t i = 0; i < ret.size(); i++)
      result[i] = ret[i];

    cvcapp.log(6,str(format("%s :: cvcstate_get_state_names(%s): num results %d\n")
                     % BOOST_CURRENT_FUNCTION
                     % string(params[0])
                     % ret.size()));

  }

  XMLRPC_METHOD_DEFINITION(cvcstate_terminate)
  {
    throw CVC_NAMESPACE::xmlrpc_server_terminate("Quitting...");
  }
}
#endif //USING_XMLRPC

namespace CVC_NAMESPACE
{
  const std::string state::SEPARATOR(".");
  state::state_ptr state::_instance;
  boost::mutex state::_instanceMutex;

  // ------------
  // state::state
  // ------------
  // Purpose: 
  //   Constructor for a state object.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.
  // 02/20/2012 -- Joe R. -- Adding notifyXmlRpc slot.
  // 03/02/2012 -- Joe R. -- Setting last mod to minimum date by default.
  // 03/15/2012 -- Joe R. -- Added initialized flag.
  // 03/30/2012 -- Joe R. -- Added hidden flag.
  state::state(const std::string& n, const state* p) :
    _name(n),
    _parent(p),
    _lastMod(boost::posix_time::min_date_time),
    _hidden(false),
    _initialized(false)
  {
    //This slot propagates child changes up to parents
    childChanged.connect(
      map_change_signal::slot_type(
        &state::notifyParent, this, _1
      )
    );

#ifdef USING_XMLRPC
    valueChanged.connect(
      signal::slot_type(
        &state::notifyXmlRpc, this
      )
    );
#endif
  }

  // -------------
  // state::~state
  // -------------
  // Purpose: 
  //   Destructor.  Just signals that this object has been destroyed.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.
  state::~state()
  {
    destroyed();
  }

  // ------------------
  // state::instancePtr
  // ------------------
  // Purpose: 
  //   Returns a pointer to the root state object singleton. Currently
  //   stores the root state object on the cvcapp datamap, though this
  //   might not always be the case.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.  
  state::state_ptr state::instancePtr()
  {
    boost::mutex::scoped_lock lock(_instanceMutex);

    if(!_instance)
      {
        //Keep our static instance in the data map!
        const std::string statekey("__state");
        state_ptr ptr(new state);
        cvcapp.data(statekey,ptr);

#ifdef USING_XMLRPC
        //Create a new XMLRPC thread to handle IPC
        cvcapp.startThread("xmlrpc_server_thread",xmlrpc_server_thread(),false);
#endif

        _instance = cvcapp.data<state_ptr>(statekey);
      }
    return _instance;
  }

  // ---------------
  // state::instance
  // ---------------
  // Purpose: 
  //   Returns a reference to the singleton root object.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.  
  state& state::instance()
  {
    return *instancePtr();
  }

  // --------------
  // state::lastMod
  // --------------
  // Purpose: 
  //   Returns the time this object was last modified.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.  
  boost::posix_time::ptime state::lastMod()
  {
    boost::this_thread::interruption_point();
    boost::mutex::scoped_lock lock(_mutex);
    return _lastMod;
  }

  // ------------
  // state::value
  // ------------
  // Purpose: 
  //   Returns the string value of this object.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.  
  std::string state::value()
  {
    boost::this_thread::interruption_point();
    boost::mutex::scoped_lock lock(_mutex);
    return _value;
  }

  // --------------------
  // state::valueTypeName
  // --------------------
  // Purpose: 
  //   Returns the type of the value as a string.
  // ---- Change History ----
  // 03/31/2012 -- Joe R. -- Initial implementation.  
  std::string state::valueTypeName()
  {
    boost::this_thread::interruption_point();
    boost::mutex::scoped_lock lock(_mutex);
    return _valueTypeName;
  }

  // --------------
  // state::comment
  // --------------
  // Purpose: 
  //   Returns the string comment for this object.
  // ---- Change History ----
  // 03/30/2012 -- Joe R. -- Initial implementation.  
  std::string state::comment()
  {
    boost::this_thread::interruption_point();
    boost::mutex::scoped_lock lock(_mutex);
    return _comment;
  }

  // --------------
  // state::comment
  // --------------
  // Purpose: 
  //   Sets a comment for this state object, useful at runtime for the user.
  // ---- Change History ----
  // 03/30/2012 -- Joe R. -- Initial implementation.  
  state& state::comment(const std::string& c)
  {
    boost::this_thread::interruption_point();
    if(comment() == c) return *this; //do nothing if equal
    
    {
      boost::mutex::scoped_lock lock(_mutex);
      _comment = c;
      _lastMod = boost::posix_time::microsec_clock::universal_time();
      _initialized = true;      
    }

    commentChanged();
    if(parent()) parent()->childChanged(name());
    return *this;
  }

  // -------------
  // state::hidden
  // -------------
  // Purpose: 
  //   Returns the hidden flag for this object.
  // ---- Change History ----
  // 03/30/2012 -- Joe R. -- Initial implementation.  
  bool state::hidden()
  {
    boost::this_thread::interruption_point();
    boost::mutex::scoped_lock lock(_mutex);
    return _hidden;
  }

  // -------------
  // state::hidden
  // -------------
  // Purpose: 
  //   Sets a hidden flag for this state object, useful to hide internal API
  //   state objects that users shouldn't change.
  // ---- Change History ----
  // 03/30/2012 -- Joe R. -- Initial implementation.  
  state& state::hidden(bool h)
  {
    boost::this_thread::interruption_point();
    if(hidden() == h) return *this; //do nothing if equal

    {
      boost::mutex::scoped_lock lock(_mutex);
      _hidden = h;
      _lastMod = boost::posix_time::microsec_clock::universal_time();
      _initialized = true;      
    }

    hiddenChanged();
    if(parent()) parent()->childChanged(name());
    return *this;
  }

  // -------------
  // state::values
  // -------------
  // Purpose: 
  //   Returns a vector of strings if the value of the object
  //   is comma separated.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.  
  std::vector<std::string> state::values(bool unique)
  {
    boost::this_thread::interruption_point();
    boost::mutex::scoped_lock lock(_mutex);
    
    using namespace std;
    using namespace boost;
    using namespace boost::algorithm;

    vector<string> vals;
    if(_value.empty()) return vals;

    string valstr = _value;
    split(vals,valstr,is_any_of(","));
    BOOST_FOREACH(string& val, vals) trim(val);
    if(unique)
      {
        set<string> vals_set;
        copy(vals.begin(), vals.end(), 
             inserter(vals_set, vals_set.begin()));
        vals.resize(vals_set.size());
        copy(vals_set.begin(), vals_set.end(),
             vals.begin());
      }
    return vals;
  }

  // ------------
  // state::value
  // ------------
  // Purpose: 
  //   Sets the value of this object.  Returns a reference to this
  //   to make it possible to add this to a chain of commands.   
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.  
  // 03/15/2012 -- Joe R. -- Added initialized flag.
  state& state::value(const std::string& v, bool setValueType)
  {
    boost::this_thread::interruption_point();
    if(value() == v) return *this; //do nothing if equal

    {
      boost::mutex::scoped_lock lock(_mutex);
      _value = v;
      _lastMod = boost::posix_time::microsec_clock::universal_time();
      _initialized = true;

      if(setValueType)
        _valueTypeName = cvcapp.dataTypeName<std::string>();
    }

    valueChanged();
    if(parent()) parent()->childChanged(name());
    return *this;
  }

  // ------------
  // state::touch
  // ------------
  // Purpose: 
  //   Triggers signals as if this state obj changed.
  // ---- Change History ----
  // 03/02/2012 -- Joe R. -- Initial implementation.
  void state::touch()
  {
    boost::this_thread::interruption_point();
    {
      boost::mutex::scoped_lock lock(_mutex);
      _lastMod = boost::posix_time::microsec_clock::universal_time();
    }
    valueChanged();
    dataChanged();
    if(parent()) parent()->childChanged(name());
  }

  // ------------
  // state::reset
  // ------------
  // Purpose: 
  //   Sets value and data to default state, and does the same for
  //   all children.
  // ---- Change History ----
  // 03/16/2012 -- Joe R. -- Creation.
  // 03/30/2012 -- Joe R. -- Resetting comment.
  void state::reset()
  {
    boost::this_thread::interruption_point();
    {
      boost::mutex::scoped_lock lock(_mutex);
      _value = std::string();
      _valueTypeName = std::string();
      _data = boost::any();
      _comment = std::string();
      _hidden = false;
      _initialized = false;
      BOOST_FOREACH(child_map::value_type val, _children)
        val.second->reset();
    }
    touch();
  }

  // ------------
  // state::ptree
  // ------------
  // Purpose: 
  //   Returns a property tree describing this state object 
  //   and it's children.  Only stores string values, not 'data'.
  // ---- Change History ----
  // 03/16/2012 -- Joe R. -- Creation.
  boost::property_tree::ptree state::ptree()
  {
    using namespace boost;
    property_tree::ptree pt;

    boost::this_thread::interruption_point();
    cvcapp.log(1,boost::str(boost::format("%s :: %s = %s\n")
                            % BOOST_CURRENT_FUNCTION
                            % fullName()
                            % value()));

    pt.push_back(property_tree::ptree::value_type(fullName(),
                                                  property_tree::ptree(value())));
    
    {
      boost::mutex::scoped_lock lock(_mutex);
      BOOST_FOREACH(child_map::value_type val, _children)
        {
          boost::this_thread::interruption_point();
          property_tree::ptree child_pt = val.second->ptree();
          pt.push_back(property_tree::ptree::value_type(val.second->fullName(), 
                                                        child_pt));
        }
    }

    return pt;
  }  

  // ------------
  // state::ptree
  // ------------
  // Purpose: 
  //  Sets this state object and its children based on an incoming
  //  property tree.
  // ---- Change History ----
  // 03/16/2012 -- Joe R. -- Creation.
  void state::ptree(const boost::property_tree::ptree& pt)
  {
    using namespace boost;
    BOOST_FOREACH(const property_tree::ptree::value_type &v, pt)
      (*this)(v.first).value(v.second.get_value<std::string>());
  }  

  // -----------
  // state::save
  // -----------
  // Purpose: 
  //  Saves this state object and its children to the specified filename.
  // ---- Change History ----
  // 03/16/2012 -- Joe R. -- Creation.
  void state::save(const std::string& filename)
  {
#ifndef CVC_STATE_XML_PROPERTY_TREE
    write_info(filename, ptree());
#else
    write_xml(filename, ptree());
#endif
  }

  // --------------
  // state::restore
  // --------------
  // Purpose: 
  //  Restores this state object and its children from the specified filename.
  // ---- Change History ----
  // 03/16/2012 -- Joe R. -- Creation.
  void state::restore(const std::string& filename)
  {
    using namespace boost;
    property_tree::ptree pt;
#ifndef CVC_STATE_XML_PROPERTY_TREE
    read_info(filename, pt);
#else
    read_xml(filename, pt);
#endif
    ptree(pt);
  }

  // ---------------
  // state::traverse
  // ---------------
  // Purpose: 
  //  Traverses the state tree, calling func for this and each child.
  //  Use 're' to filter what children get visited.
  // ---- Change History ----
  // 03/16/2012 -- Joe R. -- Creation.
  // 04/15/2012 -- Joe R. -- Triggering enter/exit signals.
  void state::traverse(traversal_unary_func func, const std::string& re)
  {
    traverseEnter();
    func(fullName());
    std::vector<std::string> ch = children(re);
    BOOST_FOREACH(std::string c, ch)
      cvcstate(c).traverse(func,re);
    traverseExit();
  }
  
  // -----------
  // state::data
  // -----------
  // Purpose: 
  //   Returns the data of this object.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.  
  boost::any state::data()
  {
    boost::this_thread::interruption_point();
    boost::mutex::scoped_lock lock(_mutex);
    return _data;
  }

  // -----------
  // state::data
  // -----------
  // Purpose: 
  //   Sets this object's data.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.
  // 03/15/2012 -- Joe R. -- Added initialized flag.
  // 04/20/2012 -- Joe R. -- Returning reference to this.
  state& state::data(const boost::any& d)
  {
    boost::this_thread::interruption_point();
    {
      boost::mutex::scoped_lock lock(_mutex);
      _data = d;
      _lastMod = boost::posix_time::microsec_clock::universal_time();
      _initialized = true;
    }
    dataChanged();
    if(parent()) parent()->childChanged(name());
    return *this;
  }

  // -------------------
  // state::dataTypeName
  // -------------------
  // Purpose: 
  //   Returns a string representing the type of the data.
  // ---- Change History ----
  // 03/31/2012 -- Joe R. -- Initial implementation.
  std::string state::dataTypeName()
  {
    return cvcapp.dataTypeName(data());
  }

  // -----------------
  // state::operator()
  // -----------------
  // Purpose: 
  //   Used for child object lookups.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.
  // 03/15/2012 -- Joe R. -- Added initialized flag.
  state& state::operator()(const std::string& childname)
  {
    using namespace std;
    using namespace boost::algorithm;

    boost::this_thread::interruption_point();

    vector<string> keys;
    split(keys, childname, is_any_of(SEPARATOR));
    if(keys.empty()) return *this;
    BOOST_FOREACH(string& key, keys) trim(key);
    //Ignore beginning empty keys
    while(!keys.empty() && keys.front().empty())
      keys.erase(keys.begin());
    if(keys.empty()) return *this;

    string nearest = keys.front();
    keys.erase(keys.begin());
    {
      boost::mutex::scoped_lock lock(_mutex);
      //If we have the child state in our map, take out its part of the
      //keys vector and recursively call its operator().
      //If not, create a new one.
      if(_children.find(nearest)!=_children.end() &&
         _children[nearest])
        return (*_children[nearest])(join(keys,SEPARATOR));
      else
        {
          state_ptr s(new state(nearest,this));
          _children[nearest] = s;
          _lastMod = boost::posix_time::microsec_clock::universal_time();
          _initialized = true;
          return (*_children[nearest])(join(keys,SEPARATOR));
        }
    }
  }

  // ---------------
  // state::children
  // ---------------
  // Purpose: 
  //   Returns a vector of children state object names. Filters children by
  //   a regular expression if regex isn't empty.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.
  // 02/24/2012 -- Joe R. -- Adding regex support.
  std::vector<std::string> state::children(const std::string& re)
  {
    boost::this_thread::interruption_point();
    boost::mutex::scoped_lock lock(_mutex);
    std::vector<std::string> ret;
    BOOST_FOREACH(child_map::value_type val, _children)
      {
        if(!re.empty())
          {
            boost::regex expression(re.c_str());
            boost::cmatch what;

            cvcapp.log(6,boost::str(boost::format("%s :: check match %s\n")
                                    % BOOST_CURRENT_FUNCTION
                                    % val.second->fullName()));

            if(boost::regex_match(val.second->fullName().c_str(),what,expression))
              {
                cvcapp.log(6,boost::str(boost::format("%s :: matched! %s\n")
                                        % BOOST_CURRENT_FUNCTION
                                        % val.second->fullName()));

                ret.push_back(val.second->fullName());
              }
          }
        else
          ret.push_back(val.second->fullName());

        //Get any matches from this state's children if any.
        std::vector<std::string> childret = 
          val.second->children(re);
        ret.insert(ret.end(), childret.begin(), childret.end());
      }
    return ret;
  }

  // ------------------
  // state::numChildren
  // ------------------
  // Purpose: 
  //   Returns the number of children.
  // ---- Change History ----
  // 04/06/2012 -- Joe R. -- Initial implementation.
  size_t state::numChildren()
  {
    boost::this_thread::interruption_point();
    boost::mutex::scoped_lock lock(_mutex);
    return _children.size();
  }

  // -------------------
  // state::notifyParent
  // -------------------
  // Purpose: 
  //   Used to propagate child change signals up the tree to the root node.
  //   Because of this, every change to the entire tree will trigger the root node's
  //   childChanged signal.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.
  void state::notifyParent(const std::string& childname)
  {
    boost::this_thread::interruption_point();
    if(parent()) parent()->childChanged(name() + SEPARATOR + childname);
  }

  // -------------------
  // state::notifyXmlRpc
  // -------------------
  // Purpose: 
  //   Used to propagate this node's changes to any network host that is listed
  //   in the __system.xmlrpc.hosts state object's value.  This spawns threads, so
  //   it will not block while the RPC call is being performed.
  // ---- Change History ----
  // 02/18/2012 -- Joe R. -- Initial implementation.
  // 02/20/2012 -- Joe R. -- Moved to its own thread to avoid possible deadlocks since
  //                          we don't yet have read/RW mutexes in use yet.
  void state::notifyXmlRpc()
  {
#ifdef USING_XMLRPC
    cvcapp.startThread("notify_xmlrpc_thread_setup",notify_xmlrpc_thread_setup(fullName()),false);
#endif
  }
}
