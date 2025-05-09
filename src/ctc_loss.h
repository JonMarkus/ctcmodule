/*
 *  ctc_loss.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef CTC_LOSS_H
#define CTC_LOSS_H

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

// Includes from sli:
#include "dictdatum.h"

namespace ctc
{
using namespace nest;

  void register_ctc_loss( const std::string& name );

/** @BeginDocumentation
Name: ctc_loss - Custom-made CTC loss for a eprop model

Description:
Coming

Parameters:

Remarks:
This is highly experimental!
 
 See ctc_readout_neuron for communication scheme with that neuron.

References:

Sends: LearningSignalConnectionEvent

Receives: LearningSignalConnectionEvent

Author:
Hans Ekkehard Plesser, Jon Markus Berg

*/
class ctc_loss : public nest::ArchivingNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model
   * manager.
   */
  ctc_loss();

  /**
   * The copy constructor is used to create model copies and instances of the
   * model.
   * @node The copy constructor needs to initialize the parameters and the
   * state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c pre_run_hook().
   */
  ctc_loss( const ctc_loss& );

  /**
   * Import sets of overloaded virtual functions.
   * This is necessary to ensure proper overload and overriding resolution.
   * @see http://www.gotw.ca/gotw/005.htm.
   */
  using nest::Node::handle;
  using nest::Node::handles_test_event;
  using Node::sends_secondary_event;
  

  /**
   * @defgroup mynest_handle Functions handling incoming events.
   * We tell nest that we can handle incoming events of various types by
   * defining @c handle() and @c connect_sender() for the given event.
   * @{
   */
  void handle( nest::InstantaneousRateConnectionEvent& ) override;
  void handle( nest::DataLoggingRequest& ) override; //! allow recording with multimeter

  size_t handles_test_event( nest::InstantaneousRateConnectionEvent&, size_t ) override;
  size_t handles_test_event( nest::DataLoggingRequest&, size_t ) override;
  /** @} */

  void
    sends_secondary_event( nest::FlexibleDataEvent& ) override
  {
  }

  void get_status( DictionaryDatum& ) const override;
  void set_status( const DictionaryDatum&) override;

private:
  //! Reset internal buffers of neuron.
  void init_buffers_() override;

  //! Initialize auxiliary quantities, leave parameters and state untouched.
  void pre_run_hook() override;

  //! Take neuron through given time interval
  void update( nest::Time const&, const long, const long ) override;

  // The next two classes need to be friends to access the State_ class/member
  friend class nest::RecordablesMap< ctc_loss >;
  friend class nest::UniversalDataLogger< ctc_loss >;

  /**
   * Free parameters of the neuron.
   *
   * These are the parameters that can be set by the user through @c SetStatus.
   * They are initialized from the model prototype when the node is created.
   * Parameters do not change during calls to @c update().
   *
   * @note Parameters_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If Parameters_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time .
   *         You may also want to define the assignment operator.
   *       - If Parameters_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
   */
  struct Parameters_
  {
    double w_stream;   //!< weighting parameter for ctc stream
    std::vector<long> target;  //! ctc target without blank
    long n_steps;    //! how many time steps to predict
    long n_target;   //! length of target
    long loss_delay; //! delay of the loss signal
    
    //! Initialize parameters to their default values.
    Parameters_();

    //! Store parameter values in dictionary.
    void get( DictionaryDatum& ) const;

    //! Set parameter values from dictionary.
    void set( const DictionaryDatum&);
  };

  /**
   * Dynamic state of the neuron.
   *
   * These are the state variables that are advanced in time by calls to
   * @c update(). In many models, some or all of them can be set by the user
   * through @c SetStatus. The state variables are initialized from the model
   * prototype when the node is created.
   *
   * @note State_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If State_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time .
   *         You may also want to define the assignment operator.
   *       - If State_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
   */
  struct State_
  {
    
    std::vector<std::vector<double>> stream;    //! steering stream
    std::vector<double> last_state;  //! state vector of the most recent sequence point
    int sequence_point;   //! current sequence point

    /**
     * Construct new default State_ instance based on values in Parameters_.
     * This c'tor is called by the no-argument c'tor of the neuron model. It
     * takes a reference to the parameters instance of the model, so that the
     * state can be initialized in accordance with parameters, e.g.,
     * initializing the membrane potential with the resting potential.
     */
    State_( const Parameters_& );

    /** Store state values in dictionary. */
    void get( DictionaryDatum& ) const;

    /**
     * Set membrane potential from dictionary.
     * @note Receives Parameters_ so it can test that the new membrane potential
     *       is below threshold.
     */
    void set( const DictionaryDatum&, const Parameters_& );
  };

  /**
   * Buffers of the neuron.
   * Ususally buffers for incoming spikes and data logged for analog recorders.
   * Buffers must be initialized by @c init_buffers_(), which is called before
   * @c pre_run_hook() on the first call to @c Simulate after the start of NEST,
   * ResetKernel.
   * @node Buffers_ needs neither constructor, copy constructor or assignment
   *       operator, since it is initialized by @c init_nodes_(). If Buffers_
   *       has members that cannot destroy themselves, Buffers_ will need a
   *       destructor.
   */
  struct Buffers_
  {
    Buffers_( ctc_loss& );
    Buffers_( const Buffers_&, ctc_loss& );

    size_t num_inputs_;  //! Number of neurons providing input
    //! Buffer for incoming predictions.
    //! Outer dim: sending neuron, inner dim: time
    std::vector<std::vector<double>> p_symbol_;

    //! Logger for all analog data
    nest::UniversalDataLogger< ctc_loss > logger_;
  };

  /**
   * Internal variables of the neuron.
   * These variables must be initialized by @c pre_run_hook, which is called before
   * the first call to @c update() upon each call to @c Simulate.
   * @node Variables_ needs neither constructor, copy constructor or assignment
   *       operator, since it is initialized by @c pre_run_hook(). If Variables_
   *       has members that cannot destroy themselves, Variables_ will need a
   *       destructor.
   */
  struct Variables_
  {
  };

  /**
   * @defgroup Access functions for UniversalDataLogger.
   * @{
   */
  /** @} */

  /**
   * @defgroup pif_members Member variables of neuron model.
   * Each model neuron should have precisely the following four data members,
   * which are one instance each of the parameters, state, buffers and variables
   * structures. Experience indicates that the state and variables member should
   * be next to each other to achieve good efficiency (caching).
   * @note Devices require one additional data member, an instance of the @c
   *       Device child class they belong to.
   * @{
   */
  Parameters_ P_; //!< Free parameters.
  State_ S_;      //!< Dynamic state.
  Variables_ V_;  //!< Internal Variables
  Buffers_ B_;    //!< Buffers.

  //! Mapping of recordables names to access functions
  static nest::RecordablesMap< ctc_loss > recordablesMap_;

  /** @} */
};

inline size_t
ctc::ctc_loss::handles_test_event( nest::InstantaneousRateConnectionEvent&, size_t receptor_type )
{
  if ( not B_.p_symbol_.empty() )
  {
    throw nest::IllegalConnection( "Cannot connect to CTC after simulation started" );
  }
  
  if ( receptor_type < 1 )
  {
    throw nest::UnknownReceptorType( receptor_type, get_name() );
  }
  
  if ( receptor_type != B_.num_inputs_ + 1 )
  {
    throw nest::IllegalConnection("must connect to CTC with strictly subsequent rports");
  }
  ++B_.num_inputs_;
  
  return receptor_type;
}

inline size_t
ctc::ctc_loss::handles_test_event( nest::DataLoggingRequest& dlr, size_t receptor_type )
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c DataLoggingRequest on port 0.
  // The function also tells the built-in UniversalDataLogger that this node
  // is recorded from and that it thus needs to collect data during simulation.
  if ( receptor_type != 0 )
  {
    throw nest::UnknownReceptorType( receptor_type, get_name() );
  }

  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
ctc_loss::get_status( DictionaryDatum& d ) const
{
  // get our own parameter and state data
  P_.get( d );
  S_.get( d );

  // get information managed by parent class
  ArchivingNode::get_status( d );

  ( *d )[ nest::names::recordables ] = recordablesMap_.get_list();
}

inline void
ctc_loss::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d);         // throws if BadProperty
  State_ stmp = S_;      // temporary copy in case of errors
  stmp.set( d, ptmp );   // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  ArchivingNode::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}

} // namespace

#endif /* #ifndef CTC_LOSS_H */
