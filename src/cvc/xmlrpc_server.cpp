#include <cvc/app.h>
#include <cvc/state.h>
#include <cvc/utility.h>

#include <xmlrpc/XmlRpc.h>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/asio.hpp>

namespace CVC_NAMESPACE
{
  CVC_DEF_EXCEPTION(xmlrpc_server_error);
  CVC_DEF_EXCEPTION(xmlrpc_server_terminate);

  // --------------------
  // notify_xmlrpc_thread
  // --------------------
  // Purpose: 
  //   Sends the value of the state object specified in 'which'
  //   to the xmlrpc server running at host:port.
  // ---- Change History ----
  // 02/20/2012 -- Joe R. -- Creation.
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
  // 02/20/2012 -- Joe R. -- Creation.
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
  // 02/20/2012 -- Joe R. -- Creation.
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

      //document the api
      cvcstate("__system.xmlrpc.port")
	.comment("The port used by the xmlrpc server.");
      cvcstate("__system.xmlrpc.hosts")
	.comment("Comma separated list of host:port combinations used to broadcast "
		 "node changes via notify_xmlrpc_thread.");
      cvcstate("__system.xmlrpc.notify_states")
	.comment("Comma separated list of nodes to broadcast notification for.");

#if 0
      if(cvcstate("__system.xmlrpc.hosts").value().empty())
        cvcstate("__system.xmlrpc.hosts")
          .value("localhost:23196"); //loopback test for now
#endif

      //cvcapp.startThread("process_notify_xmlrpc_threads",process_notify_xmlrpc_threads(),false);

      try
        {
	  using namespace boost;
          std::string host = asio::ip::host_name();
          std::string ipaddr = CVC_NAMESPACE::get_local_ip_address();

	  int port = -1;
	  try
	    {
	      port = cvcstate("__system.xmlrpc.port").value<int>();
	    }
	  catch(bad_lexical_cast&)
	    {
	      throw xmlrpc_server_error("invalid port");
	    }
	  std::string portstr = cvcstate("__system.xmlrpc.port");

          //Useful info to have
          cvcstate("__system.xmlrpc.hostname")
            .value(host)
            .comment("The hostname of the host running the xmlrpc server thread.");
          cvcstate("__system.xmlrpc.ipaddr")
            .value(ipaddr)
            .comment("The ip address bound by the xmlrpc server.");
 
	  cvcapp.log(1,str(format("%s :: %s\n")
			   % BOOST_CURRENT_FUNCTION
			   % cvcstate("__system").obj()));

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
    throw xmlrpc_server_terminate("Quitting...");
  }
}

namespace
{
  class xmlrpc_server_thread_init
  {
  public:

    static void monitor()
    {
      try
	{
	  if(!cvcstate("__system.xmlrpc").value().empty() &&
	     boost::lexical_cast<int>(cvcstate("__system.xmlrpc").value()))
	    {
	      //Create a new XMLRPC thread to handle IPC
	      cvcapp.startThread("xmlrpc_server_thread",CVC_NAMESPACE::xmlrpc_server_thread(),false);
	    }
	  else
	    {
	      if(cvcapp.hasThread("xmlrpc_server_thread"))
		cvcapp.threads("xmlrpc_server_thread")->interrupt();
	    }
	}
      catch(boost::bad_lexical_cast&)
	{
	  cvcapp.log(3,boost::str(boost::format("%s :: error parsing __system.xmlrpc\n") % BOOST_CURRENT_FUNCTION));
	}
    }

    //sets a monitor function to observe the value of __system.xmlrpc.
    //If it is set to anything that evaluates to true, the xmlrpc server thread will be started.
    //If it is set to false, the running xmlrpc server will be terminated.
    static void init()
    {
      cvcstate("__system.xmlrpc").valueChanged.connect(monitor);
      monitor();
    }

    xmlrpc_server_thread_init()
    {
      CVC_NAMESPACE::state::on_startup(init);
    }
  } static_init;
}
