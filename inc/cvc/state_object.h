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

/* $Id: StateObject.h 5883 2012-07-20 19:52:38Z transfix $ */

#include <cvc/app.h>
#include <cvc/state.h>

namespace CVC_NAMESPACE
{
  // -----------------
  // cvc::state_object
  // -----------------
  // Purpose: 
  //   This class is for making it more convenient to use the cvc::state heirarchy
  //   to store class member data for your object.
  //
  //   Use it like this:
  //
  //      class new_object : public cvc::state_object<new_object>
  //      {
  //        ...
  //      protected:
  //        virtual void handleStateChanged(const std::string& childState);
  //      };
  //
  //   Then later, you can do stuff like this:
  //
  //      new_object *p = ...
  //      p->state("member_variable").value<int>(1234);
  //
  //   The advantage of this being that you can easily monitor changes to an object's
  //   state as long as it is using the state graph to store it's data.  You can also
  //   change an object's state easily and it should respond to each change in a meaningful
  //   way.
  //
  // ---- Change History ----
  // 05/27/2012 -- Joe R. -- Creation.
  template <class This> //This should be the type of the inheriting class
  class state_object
  {
  public:
    state_object() 
    { 
      cvcapp.registerDataType(This);

      //watch this object's state
      _stateConnection = getState().childChanged.connect(
         map_change_signal::slot_type(
           &state_object<This>::stateChanged, this, _1
         )
      );
    }

    ~state_object() { _stateConnection.disconnect(); }

    //Use this to easily get the name of this viewer's state object.
    std::string stateName(const std::string& childState = std::string()) const { 
      std::string viewer_root = cvcapp.dataTypeName<This>()+
        CVC::State::SEPARATOR+
        boost::lexical_cast<std::string>(this);
      return 
        !childState.empty() ? 
        viewer_root + CVC::State::SEPARATOR + childState :
        viewer_root;
    }

    //Shortcut for accessing the state corresponding to an instance of this
    //class or it's children.
    state& getState(const std::string& s = std::string()) const {
      return cvcstate(stateName(s));
    }

  protected:
    boost::signals2::connection _stateConnection;

    //Classes that are state_objects should implement this function for themselves.
    //Note: each call happens in its own thread.
    virtual void handleStateChanged(const std::string& childState)
    {
      cvcapp.log(2,str(boost::format("%s :: state changed: %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % childState));
    }

  private:
    //Responding to state changes.  Every change will launch a new thread and will
    //immediately return, therefore this call is non-blocking.
    void stateChanged(const std::string& childState)
    {
      cvcapp.startThread(stateName(childState) + "_stateChanged",
                         boost::bind(&state_object<This>::handleStateChanged, 
                                     boost::ref(*this),
                                     childState));
    }
  };
}
