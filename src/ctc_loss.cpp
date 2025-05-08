/*
 *  ctc_loss.cpp
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

#include "ctc_loss.h"

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "nest_impl.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "sharedptrdatum.h"

using namespace nest;

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< ctc::ctc_loss > ctc::ctc_loss::recordablesMap_;

void
ctc::register_ctc_loss( const std::string& name )
{
  register_node_model< ctc_loss >( name );
}

namespace nest
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< ctc::ctc_loss >::create()
{
  // use standard names whereever you can for consistency!
  // insert_( names::V_m, &ctc::ctc_loss::get_V_m_ );
}
}


/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

ctc::ctc_loss::Parameters_::Parameters_()
  : w_stream(1)
  , target({0, 1, 0, 2, 0})
  , n_steps(5000)
  , n_target(target.size())
  , loss_delay(3)
{
}

ctc::ctc_loss::State_::State_( const Parameters_& p )
: stream(p.target.size(), std::vector<double>(p.n_steps, 0.0))
, last_state(std::vector<double>(p.target.size(), 0.0))
, sequence_point(0)
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
ctc::ctc_loss::Parameters_::get( DictionaryDatum& d ) const
{
  def<double>(d,  "w_stream", w_stream );
  def<std::vector<long>>(d, "target", target);
  def<long>(d, "n_steps", n_steps);
}

void
ctc::ctc_loss::Parameters_::set( const DictionaryDatum& d)
{
  // add check that num input ch cannot be changed after making connections or simulating
  updateValue< double >( d, "w_stream", w_stream );
  updateValue< long >( d, "n_steps", n_steps);
  updateValue< std::vector<long> >(d, "target", target );
  
 
}

void
ctc::ctc_loss::State_::get( DictionaryDatum& d ) const
{
  // def< std::vector<std::vector<double>> >( d, "stream", stream );
  // def< std::vector<double>>(d, "last_state", last_state);
  // def<long>(d, "sequence_point", sequence_point);
}

void
ctc::ctc_loss::State_::set( const DictionaryDatum& d, const Parameters_& p )
{
}

ctc::ctc_loss::Buffers_::Buffers_( ctc_loss& n )
  : num_inputs_( 0 )
  , logger_( n )
{
}

ctc::ctc_loss::Buffers_::Buffers_( const Buffers_& b, ctc_loss& n )
  : num_inputs_(b.num_inputs_)
  , logger_( n )
{
}


/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

ctc::ctc_loss::ctc_loss()
  : ArchivingNode()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
}

ctc::ctc_loss::ctc_loss( const ctc_loss& n )
  : ArchivingNode( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
ctc::ctc_loss::init_buffers_()
{
  B_.p_symbol_.clear();
  B_.logger_.reset();  // includes resize
}

void
ctc::ctc_loss::pre_run_hook()
{
  S_.sequence_point = 0;
  P_.n_target = P_.target.size();

  B_.p_symbol_.resize( B_.num_inputs_ );
  for ( auto& ps : B_.p_symbol_ )
  {
    ps.resize( kernel().connection_manager.get_min_delay(), 0 );
  }

  P_.n_target = P_.target.size();


  if(P_.n_steps < P_.n_target/2){
    throw BadParameterValue("n_steps is to small for this target");
  }

  // Initialize the alpha matrix
  std::vector<std::vector<double>> alpha(P_.n_target, std::vector<double>(P_.n_steps, 0.0));
  alpha[0][0] = 0.5;
  alpha[1][0] = 0.5;

  // Compute alpha values
  for (int j = 1; j < P_.n_steps; j++) {
      for (int i = 0; i < P_.n_target; i++) {
          if (i == 0) {
              alpha[i][j] = alpha[i][j - 1];
          } else if (i == 1) {
              alpha[i][j] = alpha[i][j - 1] + alpha[i - 1][j - 1];
          } else if (P_.target[i] == P_.target[i - 2]) {
              alpha[i][j] = alpha[i][j - 1] + alpha[i - 1][j - 1];
          } else {
              alpha[i][j] = alpha[i][j - 1] + alpha[i - 1][j - 1] + alpha[i - 2][j - 1];
          }
      }
  }

  // Initialize the beta matrix
  std::vector<std::vector<double>> beta(P_.n_target, std::vector<double>(P_.n_steps, 0.0));
  beta[P_.n_target - 1][P_.n_steps - 1] = 0.5;
  beta[P_.n_target - 2][P_.n_steps - 1] = 0.5;

  // Compute beta values
  for (int j = P_.n_steps - 2; j >= 0; j--) {
      for (int i = P_.n_target - 1; i >= 0; i--) {
          if (i == P_.n_target - 1) {
              beta[i][j] = beta[i][j + 1];
          } else if (i == P_.n_target - 2) {
              beta[i][j] = beta[i][j + 1] + beta[i + 1][j + 1];
          } else if (P_.target[i] == P_.target[i + 2]) {
              beta[i][j] = beta[i][j + 1] + beta[i + 1][j + 1];
          } else {
              beta[i][j] = beta[i][j + 1] + beta[i + 1][j + 1] + beta[i + 2][j + 1];
          }
      }
  }

  S_.stream.clear();
  S_.stream.resize(P_.n_target, std::vector<double>(P_.n_steps, 0.0));
  

  // Compute the stream matrix by multiplying alpha and beta element-wise
  for (int i = 0; i < P_.n_target; i++) {
      for (int j = 0; j < P_.n_steps; j++) {
        S_.stream.at(i).at(j) = alpha[i][j] * beta[i][j];
      }
  }

  // Normalize each column of the stream matrix
  for (int j = 0; j < P_.n_steps; j++) {
      double colSum = 0.0;
      for (int i = 0; i < P_.n_target; i++) {
          colSum += S_.stream[i][j];
      }
      if (colSum > 0.0) { // Avoid division by zero
          for (int i = 0; i < P_.n_target; i++) {
            S_.stream.at(i).at(j) /= colSum;
          }
      }
  }


  B_.logger_.init();
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
ctc::ctc_loss::update( Time const& slice_origin, const long from_step, const long to_step )
{

  // time zero: slice_origin.get_steps() == 0 and from_step == 0

  // outer dimension: inputs/targets, inner dimension: time, but flattened as needed by set_coeffarray
  const size_t min_delay = kernel().connection_manager.get_min_delay();
  std::vector< double > loss_buffer( B_.num_inputs_ * min_delay, 0 );

  const long update_interval = kernel().simulation_manager.get_eprop_update_interval().get_steps();
  const long learning_window = kernel().simulation_manager.get_eprop_learning_window().get_steps();


  for ( long lag = from_step; lag < to_step; ++lag )
  {
    const long t = slice_origin.get_steps() + lag;
    const long interval_step_signals = ( t - P_.loss_delay ) % update_interval;

    if( interval_step_signals  == update_interval - learning_window )
    {
    
    

    std::vector<double> prediction;
    for ( const auto& pa : B_.p_symbol_ )
    {
        prediction.push_back(pa[lag]);
    }

    // no prdiction can be less than zero
    for (auto& p : prediction) 
    {
        if (p < 0) 
        {
            p = 0;
        }
    }

    // make prediction to probabilities
    double sum_pred = std::accumulate(prediction.begin(), prediction.end(), 0.0);
    if (sum_pred != 0) 
    {
      for (auto& p : prediction) 
      {
          p /= sum_pred;
      }
    }
    else
    {
      // warning("Prediction sum is zero");
      std::cerr << "Prediction sum is zero" << std::endl;
      for (auto& p : prediction) 
      {
          p = 1.0 / prediction.size();
      }
      sum_pred = 1.0;
    }


    // Forward state
    std::vector<double> forward_state(P_.n_target, 0.0);

    // Step target
    std::vector<double> step_target(prediction.size(), 0.0);

    if (S_.sequence_point == 0) 
    {
        forward_state[0] = prediction[P_.target[0]];
        forward_state[1] = prediction[P_.target[1]];
    } 
    else 
    {
        std::vector<double> node_sums(prediction.size(), 0.0);
        for (int i = 0; i < P_.n_target; ++i) 
        {
            int e = P_.target[i];
            if (i == 0) 
            {
                forward_state[i] = S_.last_state[i];
            } 
            else if (i == 1) 
            {
                forward_state[i] = S_.last_state[i] + S_.last_state[i - 1];
            } 
            else if (P_.target[i] == P_.target[i - 2]) 
            {
                forward_state[i] = S_.last_state[i] + S_.last_state[i - 1];
            } 
            else 
            {
                forward_state[i] = S_.last_state[i] + S_.last_state[i - 1] + S_.last_state[i - 2];
            }
            node_sums[e] += forward_state[i];
        }

        for (int i = 0; i < P_.n_target; ++i) 
        {
            int e = P_.target[i];
            if (node_sums[e] > 0) 
            {
                forward_state[i] = forward_state[i] * forward_state[i] / node_sums[e] * prediction[e];
            }
        }
    }
    // Normalize forward_state
    double sum_forward = std::accumulate(forward_state.begin(), forward_state.end(), 0.0);
    for (auto& f : forward_state) 
    {
        f /= sum_forward;
    }
    

    // Compute ctc_state
    std::vector<double> ctc_state(P_.n_target, 0.0);
    for (int i = 0; i < P_.n_target; ++i) 
    {
      ctc_state[i] = forward_state[i] * pow(S_.stream.at(i).at(S_.sequence_point), P_.w_stream);
    }
 

    // Normalize ctc_state
    double sum_ctc = std::accumulate(ctc_state.begin(), ctc_state.end(), 0.0);

    if (sum_ctc == 0)
    {
      // warning("CTC sum is zero");
      std::cerr << "CTC sum is " << sum_ctc << std::endl;
      for(int i = 0; i < P_.n_target; ++i )
      {
        ctc_state[i] = S_.stream[i][S_.sequence_point];
      }
  }
    else
    {

        for (auto& c : ctc_state) 
        {
            c /= sum_ctc;
        }
    } 

    
    // Update S_.last_state
    S_.last_state = forward_state;

    // Update step target
    for (int i = 0; i < P_.n_target; ++i) 
    {
        step_target[P_.target[i]] += ctc_state[i];
    }
    
    loss_buffer.at(0 + lag) = 0;

     // filling loss buffer here with dummy values for testing
     for ( size_t i = 1 ; i < B_.num_inputs_ ; ++i )
     {
       loss_buffer.at( i * min_delay + lag ) = (prediction.at( i ) - step_target.at( i )) * sum_pred;
     }

    S_.sequence_point ++;

    B_.logger_.record_data( slice_origin.get_steps() + lag );
    }
  }
  

  
  // send loss back to readout neuron, blocked by neuron
  FlexibleDataEvent loss_event;
  loss_event.set_coeffarray( loss_buffer );
  kernel().event_delivery_manager.send_secondary( *this, loss_event );
}

void
ctc::ctc_loss::handle( InstantaneousRateConnectionEvent& e )
{
  auto it_event = e.begin();
  const long channel = e.get_rport();

   for ( size_t i = 0 ; it_event != e.end() ; ++i )
  {
    const double prediction = e.get_coeffvalue( it_event ); // get_coeffvalue advances iterator

    assert( 0 < channel and channel <= B_.num_inputs_ );
    B_.p_symbol_.at( channel - 1 ).at( i ) = prediction;
  }
}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void
ctc::ctc_loss::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e ); // the logger does this for us
}
