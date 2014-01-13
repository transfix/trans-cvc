#include <cvc/app.h>
#include <cvc/state.h>
#include <cvc/utility.h>

#include <xmlrpc/XmlRpc.h>

#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/asio.hpp>

namespace CVC_NAMESPACE
{
  // --------------------
  // xmlrpc_client_thread
  // --------------------
  // Purpose: 
  //   Performs a client xmlrpc call to the specified method.
  // ---- Change History ----
  // 01/12/2014 -- Joe R. -- Creation.
  class xmlrpc_client_thread
  {
  public:
    xmlrpc_client_thread(const std::string& host, int port,
                         const std::string& method_name, const XmlRpc::XmlRpcValue& params,
                         boost::posix_time::ptime mod_time = boost::posix_time::microsec_clock::universal_time())
      : _host(host), 
        _port(port), _method_name(method_name), 
        _params(params), _mod_time(mod_time) 
    {
      _thread_name = boost::str(boost::format("xmlrpc_client_thread_%s_%d_%s")
                                % host % port % method_name);
    }
    
    xmlrpc_client_thread(const xmlrpc_client_thread& t)
      : _thread_name(t._thread_name), _host(t._host), 
        _port(t._port), _method_name(t._method_name), 
        _params(t._params), _mod_time(t._mod_time) {}
    
    xmlrpc_client_thread& operator=(const xmlrpc_client_thread& t)
    {
      _thread_name = t._thread_name;
      _host = t._host;
      _port = t._port;
      _method_name = t._method_name;
      _params = t._params;
      _mod_time = t._mod_time;
      return *this;
    }
    
    void operator()()
    {
      using namespace boost;
      using namespace std;

      using namespace boost;
      using namespace std;
      
      CVC_NAMESPACE::thread_feedback feedback;

      XmlRpc::XmlRpcClient c(_host.c_str(),_port);
      XmlRpc::XmlRpcValue params,result;

      params = _params;

      //params[0] = _params;

      //this might be ignored if the server has a newer value
      //      params[1] = posix_time::to_simple_string(_mod_time);


      for(int i = 0; i < 3; i++)
        cvcapp.log(6,str(format("%s :: params[%d] = %s\n")
                         % BOOST_CURRENT_FUNCTION
                         % i
                         % string(params[i])));

      c.execute(_method_name.c_str(), params, result);

      cvcapp.data(thread_name() + "_result", result);
    }

    const std::string& thread_name() const { return _thread_name; }

  protected:
    std::string _thread_name;
    std::string _host;
    int _port;
    std::string _method_name;
    XmlRpc::XmlRpcValue _params;
    boost::posix_time::ptime _mod_time;
  };

  // --------------------------
  // queue_xmlrpc_client_thread
  // --------------------------
  // Purpose: 
  //   Queues up a thread on the data map to set a remote state value at a later time.
  //   process_xmlrpc_client_threads is in charge of actually launching the thread.
  // ---- Change History ----
  // 01/13/2014 -- Joe R. -- Creation.
  class queue_xmlrpc_client_thread
  {
  public:
    queue_xmlrpc_client_thread(const xmlrpc_client_thread& xct)
      : _xct(xct) {}
    
    queue_xmlrpc_client_thread(const queue_xmlrpc_client_thread& q)
      : _xct(q._xct) {}

    queue_xmlrpc_client_thread& operator=(const queue_xmlrpc_client_thread& rhs)
    {
      _xct = rhs._xct;
    }

    void operator()()
    {
      CVC_NAMESPACE::thread_feedback feedback;

      //Stick the thread on the datamap.  There is another thread that will actually launch
      //these threads some time later.  Doing this, we won't flood the network with xmlrpc
      //requests if data gets old before it is sent out over the network.
      cvcapp.data(_xct.thread_name(), _xct);
    }
  private:
    xmlrpc_client_thread _xct;
  };

  // -----------------------------
  // process_xmlrpc_client_threads
  // -----------------------------
  // Purpose: 
  //   Starts all the xmlrpc client threads that are on the data map.
  // ---- Change History ----
  // 01/13/2014 -- Joe R. -- creation.
  class process_xmlrpc_client_threads
  {
  public:
    void operator()()
    {
      CVC_NAMESPACE::thread_feedback feedback;

      while(1)
        {
          //Sleep for 200ms before each iteration.
          cvcapp.sleep(cvcstate("__system.xmlrpc_client.interval").value<double>());
          
          std::vector<std::string> keys = cvcapp.data<xmlrpc_client_thread>();
          std::vector<xmlrpc_client_thread> threads = 
            cvcapp.data<xmlrpc_client_thread>(keys);
          BOOST_FOREACH(xmlrpc_client_thread& thread, threads)
            {
              if(cvcapp.hasThread(thread.thread_name()))
                continue;
              cvcapp.startThread(thread.thread_name(),thread);
              cvcapp.data(thread.thread_name(),boost::any()); //erase from the datamap
            }
        }
    }
  };
}

namespace
{
  class process_xmlrpc_client_threads_init
  {
  public:

    static void init()
    {
      cvcstate("__system.xmlrpc_client.interval").value(double(200.0));
      cvcapp.startThread("process_xmlrpc_client_threads",
			 CVC_NAMESPACE::process_xmlrpc_client_threads());
    }

    process_xmlrpc_client_threads_init()
    {
      CVC_NAMESPACE::state::on_startup(init);
    }
  } static_init;
}
