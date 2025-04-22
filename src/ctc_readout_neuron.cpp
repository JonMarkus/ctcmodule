/*
 *  ctc_readout_neuron.cpp
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

// nest models
#include "ctc_readout_neuron.h"

// C++
#include <limits>

// libnestutil
#include "dict_util.h"
#include "numerics.h"

// nestkernel
#include "exceptions.h"
#include "kernel_manager.h"
#include "nest_impl.h"
#include "universal_data_logger_impl.h"

// sli
#include "dictutils.h"

void
ctc::register_ctc_readout_neuron( const std::string& name )
{
  register_node_model< ctc_readout_neuron >( name );
}

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< ctc::ctc_readout_neuron > ctc::ctc_readout_neuron::recordablesMap_;

namespace nest {
template <>
void
RecordablesMap< ctc::ctc_readout_neuron >::create()
{
  insert_( names::eprop_history_duration, &ctc::ctc_readout_neuron::get_eprop_history_duration );
  insert_( names::error_signal, &ctc::ctc_readout_neuron::get_error_signal_ );
  insert_( names::readout_signal, &ctc::ctc_readout_neuron::get_readout_signal_ );
  insert_( names::readout_signal_unnorm, &ctc::ctc_readout_neuron::get_readout_signal_unnorm_ );
  insert_( "p_symbol", &ctc::ctc_readout_neuron::get_p_symbol_ );
  insert_( names::target_signal, &ctc::ctc_readout_neuron::get_target_signal_ );
  insert_( names::V_m, &ctc::ctc_readout_neuron::get_v_m_ );
}
}

/* ----------------------------------------------------------------
 * Default constructors for parameters, state, and buffers
 * ---------------------------------------------------------------- */

ctc::ctc_readout_neuron::Parameters_::Parameters_()
  : C_m_( 250.0 )
  , E_L_( 0.0 )
  , I_e_( 0.0 )
  , loss_( "mean_squared_error" )
  , regular_spike_arrival_( true )
  , tau_m_( 10.0 )
  , V_min_( -std::numeric_limits< double >::max() )
  , training_mode(1)
{
}

ctc::ctc_readout_neuron::State_::State_()
  : error_signal_( 0.0 )
  , readout_signal_( 0.0 )
  , readout_signal_unnorm_( 0.0 )
  , p_symbol_( 0.0 )
  , target_signal_( 0.0 )
  , i_in_( 0.0 )
  , v_m_( 0.0 )
  , z_in_( 0.0 )
{
}

ctc::ctc_readout_neuron::Buffers_::Buffers_( ctc::ctc_readout_neuron& n )
  : logger_( n )
{
}

ctc::ctc_readout_neuron::Buffers_::Buffers_( const Buffers_&, ctc::ctc_readout_neuron& n )
  : logger_( n )
{
}

/* ----------------------------------------------------------------
 * Getter and setter functions for parameters and state
 * ---------------------------------------------------------------- */

void
ctc::ctc_readout_neuron::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::C_m, C_m_ );
  def< double >( d, names::E_L, E_L_ );
  def< double >( d, names::I_e, I_e_ );
  def< std::string >( d, names::loss, loss_ );
  def< bool >( d, names::regular_spike_arrival, regular_spike_arrival_ );
  def< double >( d, names::tau_m, tau_m_ );
  def< double >( d, names::V_min, V_min_ + E_L_ );
  def< bool >( d, "training_mode", training_mode );
}

double
ctc::ctc_readout_neuron::Parameters_::set( const DictionaryDatum& d, Node* node )
{
  // if leak potential is changed, adjust all variables defined relative to it
  const double ELold = E_L_;
  updateValueParam< double >( d, names::E_L, E_L_, node );
  const double delta_EL = E_L_ - ELold;

  V_min_ -= updateValueParam< double >( d, names::V_min, V_min_, node ) ? E_L_ : delta_EL;

  updateValueParam< double >( d, names::C_m, C_m_, node );
  updateValueParam< double >( d, names::I_e, I_e_, node );
  updateValueParam< std::string >( d, names::loss, loss_, node );
  updateValueParam< bool >( d, names::regular_spike_arrival, regular_spike_arrival_, node );
  updateValueParam< double >( d, names::tau_m, tau_m_, node );
  updateValueParam< bool >( d, "training_mode", training_mode, node );

  if ( C_m_ <= 0 )
  {
    throw BadProperty( "Membrane capacitance C_m > 0 required." );
  }

  if ( loss_ != "mean_squared_error" and loss_ != "cross_entropy" )
  {
    throw BadProperty( "Loss function loss from [\"mean_squared_error\", \"cross_entropy\"] required." );
  }

  if ( tau_m_ <= 0 )
  {
    throw BadProperty( "Membrane time constant tau_m > 0 required." );
  }

  return delta_EL;
}

void
ctc::ctc_readout_neuron::State_::get( DictionaryDatum& d, const Parameters_& p ) const
{
  def< double >( d, names::V_m, v_m_ + p.E_L_ );
  def< double >( d, names::error_signal, error_signal_ );
  def< double >( d, names::readout_signal, readout_signal_ );
  def< double >( d, names::readout_signal_unnorm, readout_signal_unnorm_ );
  def< double >( d, "p_symbol", p_symbol_ );
  def< double >( d, names::target_signal, target_signal_ );
}

void
ctc::ctc_readout_neuron::State_::set( const DictionaryDatum& d, const Parameters_& p, double delta_EL, Node* node )
{
  v_m_ -= updateValueParam< double >( d, names::V_m, v_m_, node ) ? p.E_L_ : delta_EL;
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

ctc::ctc_readout_neuron::ctc_readout_neuron()
  : EpropArchivingNodeReadout()
  , P_()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
}

ctc::ctc_readout_neuron::ctc_readout_neuron( const ctc_readout_neuron& n )
  : EpropArchivingNodeReadout( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
ctc::ctc_readout_neuron::init_buffers_()
{
  B_.normalization_rate_ = 0;
  B_.spikes_.clear();   // includes resize
  B_.currents_.clear(); // includes resize
  B_.logger_.reset();   // includes resize
  
  B_.ctc_loss_.resize( kernel().connection_manager.get_min_delay(), 0.0 );
}

void
ctc::ctc_readout_neuron::pre_run_hook()
{
  B_.logger_.init(); // ensures initialization in case multimeter connected after Simulate

  if ( P_.loss_ == "mean_squared_error" )
  {
    compute_error_signal = &ctc_readout_neuron::compute_error_signal_mean_squared_error;
    V_.signal_to_other_readouts_ = false;
  }
  else if ( P_.loss_ == "cross_entropy" )
  {
    compute_error_signal = &ctc_readout_neuron::compute_error_signal_cross_entropy;
    V_.signal_to_other_readouts_ = true;
  }

  const double dt = Time::get_resolution().get_ms();

  V_.P_v_m_ = std::exp( -dt / P_.tau_m_ );
  V_.P_i_in_ = P_.tau_m_ / P_.C_m_ * ( 1.0 - V_.P_v_m_ );
  V_.P_z_in_ = P_.regular_spike_arrival_ ? 1.0 : 1.0 - V_.P_v_m_;
}


/* ----------------------------------------------------------------
 * Update function
 * ---------------------------------------------------------------- */

void
ctc::ctc_readout_neuron::update( Time const& origin, const long from, const long to )
{
  const long update_interval = kernel().simulation_manager.get_eprop_update_interval().get_steps();
  const long learning_window = kernel().simulation_manager.get_eprop_learning_window().get_steps();
  const bool with_reset = kernel().simulation_manager.get_eprop_reset_neurons_on_update();
  const long shift = get_shift();

  const size_t buffer_size = kernel().connection_manager.get_min_delay();

  std::vector< double > error_signal_buffer( buffer_size, 0.0 );
  std::vector< double > readout_signal_unnorm_buffer( buffer_size, 0.0 );
  std::vector< double > p_symbol_buffer( buffer_size, 0.0 );
  
  for ( long lag = from; lag < to; ++lag )
  {
    const long t = origin.get_steps() + lag;
    const long interval_step = ( t - shift ) % update_interval;
    const long interval_step_signals = ( t - shift - delay_out_norm_ ) % update_interval;

    if ( interval_step == 0 )
    {
      erase_used_eprop_history();

      if ( with_reset )
      {
        S_.v_m_ = 0.0;
      }
    }

    S_.z_in_ = B_.spikes_.get_value( lag );

    S_.v_m_ = V_.P_i_in_ * S_.i_in_ + V_.P_z_in_ * S_.z_in_ + V_.P_v_m_ * S_.v_m_;
    S_.v_m_ = std::max( S_.v_m_, P_.V_min_ );

    // S_.error_signal_ = B_.ctc_loss_[ lag ];

    // error signal calculation should somehow use B_.ctc_loss_.
    // ( this->*compute_error_signal )( lag );
    // const double norm_rate = B_.normalization_rate_ + S_.readout_signal_unnorm_;
    // S_.readout_signal_unnorm_ = S_.v_m_ + P_.E_L_;
    S_.readout_signal_unnorm_ = std::exp( S_.v_m_ + P_.E_L_ );
    S_.readout_signal_ = S_.readout_signal_unnorm_ ;

    // write dummy value for testing, build in lag so we can see lags arrive correctly
    // p_symbol_buffer[ lag ] = S_.readout_signal_;

    if(interval_step + 1 == update_interval - learning_window)
    {
      // send signal to loss 
     p_symbol_buffer[ lag ] = S_.readout_signal_;
     
     S_.error_signal_ = 0.0;


    }

    else if ( interval_step_signals < update_interval - learning_window )
    {
      S_.target_signal_ = 0.0;
      S_.readout_signal_ = 0.0;
      S_.error_signal_ = 0.0;
    }
    else
    {
      
      S_.readout_signal_ = 0.0;
      if(P_.training_mode)
      {
        S_.error_signal_ = B_.ctc_loss_[ lag ];
      }
      else
      {
        S_.error_signal_ = 0.0;
      }
    }

    B_.normalization_rate_ = 0.0;

    // print out dummy values for testing. This is stuff received from ctc_readout
    // std::cerr << std::fixed << std::setprecision(0);
    // std::cerr << "readout update loop: t = " << origin.get_steps() << ", gid = " << get_node_id()
    //   << ", loss received = ";
    // for ( const auto&l : B_.ctc_loss_ )
    // {
    //   std::cerr << l << " ";
    // }
    // std::cerr << std::endl;
    // end dummy printing
    
    if ( V_.signal_to_other_readouts_ )
    {
      readout_signal_unnorm_buffer[ lag ] = S_.readout_signal_unnorm_;
    }

    error_signal_buffer[ lag ] = S_.error_signal_;
    


    append_new_eprop_history_entry( t );
    write_error_signal_to_history( t, S_.error_signal_ );

    S_.i_in_ = B_.currents_.get_value( lag ) + P_.I_e_;

    B_.logger_.record_data( t );
  }

  LearningSignalConnectionEvent error_signal_event;
  error_signal_event.set_coeffarray( error_signal_buffer );
  kernel().event_delivery_manager.send_secondary( *this, error_signal_event );

  if ( V_.signal_to_other_readouts_ )
  {
    // time is one time step longer than the final interval_step to enable sending the
    // unnormalized readout signal one time step in advance so that it is available
    // in the next times step for computing the normalized readout signal
    DelayedRateConnectionEvent readout_signal_unnorm_event;
    readout_signal_unnorm_event.set_coeffarray( readout_signal_unnorm_buffer );
    kernel().event_delivery_manager.send_secondary( *this, readout_signal_unnorm_event );
  }
  
  // send p_symbol to CTC
  InstantaneousRateConnectionEvent p_symbol_event;
  p_symbol_event.set_coeffarray( p_symbol_buffer );
  kernel().event_delivery_manager.send_secondary( *this, p_symbol_event );

  // print out dummy values for testing. This is stuff received from ctc_readout
  // std::cerr << "readout update end: t =  " << origin.get_steps() << ", gid = " << get_node_id()
  //   << ", p_symbol sending = ";
  // for ( const auto&p : p_symbol_buffer )
  // {
  //   std::cerr << p << " ";
  // }
  // std::cerr << std::endl;
  // end dummy printing

  return;
}

/* ----------------------------------------------------------------
 * Error signal functions
 * ---------------------------------------------------------------- */

void
ctc::ctc_readout_neuron::compute_error_signal_mean_squared_error( const long lag )
{
  S_.readout_signal_ = S_.readout_signal_unnorm_;
  S_.readout_signal_unnorm_ = S_.v_m_ + P_.E_L_;
  S_.error_signal_ = S_.readout_signal_ - S_.target_signal_;
}

void
ctc::ctc_readout_neuron::compute_error_signal_cross_entropy( const long lag )
{
  const double norm_rate = B_.normalization_rate_ + S_.readout_signal_unnorm_;
  S_.readout_signal_ = S_.readout_signal_unnorm_ / norm_rate;
  S_.readout_signal_unnorm_ = std::exp( S_.v_m_ + P_.E_L_ );
  S_.error_signal_ = S_.readout_signal_ - S_.target_signal_;
}

/* ----------------------------------------------------------------
 * Event handling functions
 * ---------------------------------------------------------------- */

void
ctc::ctc_readout_neuron::handle( DelayedRateConnectionEvent& e )
{
  const size_t rport = e.get_rport();
  assert( rport < SUP_RATE_RECEPTOR );

  auto it = e.begin();
  assert( it != e.end() );

  const double signal = e.get_weight() * e.get_coeffvalue( it );
  if ( rport == READOUT_SIG )
  {
    B_.normalization_rate_ += signal;
  }
  else if ( rport == TARGET_SIG )
  {
    S_.target_signal_ = signal;
  }

  assert( it == e.end() );
}

void
ctc::ctc_readout_neuron::handle( FlexibleDataEvent& e )
{
  // Extract the block that is relevant for this neuron.
  const size_t rport = e.get_rport();
  const size_t min_delay = kernel().connection_manager.get_min_delay();
  
  // make sure we have min_delay entries per rport
  assert( rport > 0 );
  assert( rport * min_delay <= e.size() );

  // Advance to section of buffer with data for this readout neuron
  // Take into account that iterator for e is based on underlying uint type
  auto it = e.begin()
    + number_of_uints_covered< decltype(B_.ctc_loss_)::value_type >() * (rport-1) * min_delay;
  assert( it != e.end() );

  for ( size_t i = 0 ; i < min_delay and it != e.end() ; ++i )
  {
    // get_coeffvalue() advances it
    B_.ctc_loss_.at( i ) = e.get_coeffvalue( it );
  }
}


void
ctc::ctc_readout_neuron::handle( SpikeEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  B_.spikes_.add_value(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), e.get_weight() * e.get_multiplicity() );
}

void
ctc::ctc_readout_neuron::handle( CurrentEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  B_.currents_.add_value(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), e.get_weight() * e.get_current() );
}

void
ctc::ctc_readout_neuron::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

double
ctc::ctc_readout_neuron::compute_gradient( std::vector< long >& presyn_isis,
  const long,
  const long t_previous_trigger_spike,
  const double kappa,
  const bool average_gradient )
{
  auto eprop_hist_it = get_eprop_history( t_previous_trigger_spike );

  double grad = 0.0;  // gradient value to be calculated
  double L = 0.0;     // error signal
  double z = 0.0;     // spiking variable
  double z_bar = 0.0; // low-pass filtered spiking variable

  for ( long presyn_isi : presyn_isis )
  {
    z = 1.0; // set spiking variable to 1 for each incoming spike

    for ( long t = 0; t < presyn_isi; ++t )
    {
      assert( eprop_hist_it != eprop_history_.end() );

      L = eprop_hist_it->error_signal_;

      z_bar = V_.P_v_m_ * z_bar + V_.P_z_in_ * z;
      grad += L * z_bar;
      z = 0.0; // set spiking variable to 0 between spikes

      ++eprop_hist_it;
    }
  }
  presyn_isis.clear();

  const long learning_window = kernel().simulation_manager.get_eprop_learning_window().get_steps();
  if ( average_gradient )
  {
    grad /= learning_window;
  }

  return grad;
}

