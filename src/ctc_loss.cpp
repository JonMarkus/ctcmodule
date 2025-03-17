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
  : num_input_channels(1)
  , output_channel(1)
  , w_stream(1)
  , target()
{
}

ctc::ctc_loss::State_::State_( const Parameters_& p )
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
ctc::ctc_loss::Parameters_::get( DictionaryDatum& d ) const
{
  def<long>(d, "num_input_channels", num_input_channels );
  def<long>(d, "output_channel", output_channel );
  def<double>(d,  "w_stream", w_stream );
  def<std::string>(d, "target", target);
}

void
ctc::ctc_loss::Parameters_::set( const DictionaryDatum& d )
{
  // add check that num input ch cannot be changed after making connections or simulating
  
  updateValue< long >( d, "num_input_channels", num_input_channels );
  updateValue< long >( d, "output_channel", output_channel );
  updateValue< double >( d, "w_stream", w_stream );
  updateValue< std::string >(d, "target", target );
  
  if ( num_input_channels < 1 )
  {
    throw BadProperty("num_input_channels > 0 required");
  }

  if ( output_channel < 1 or output_channel > num_input_channels )
  {
    throw BadProperty("Bad output channel");
  }
}

void
ctc::ctc_loss::State_::get( DictionaryDatum& d ) const
{
}

void
ctc::ctc_loss::State_::set( const DictionaryDatum& d, const Parameters_& p )
{
}

ctc::ctc_loss::Buffers_::Buffers_( ctc_loss& n )
  : logger_( n )
{
}

ctc::ctc_loss::Buffers_::Buffers_( const Buffers_&, ctc_loss& n )
  : logger_( n )
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
  B_.predictions.clear();   // includes resize
  B_.logger_.reset();  // includes resize
}

void
ctc::ctc_loss::pre_run_hook()
{
  B_.predictions.resize( P_.num_input_channels );
  B_.logger_.init();
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
ctc::ctc_loss::update( Time const& slice_origin, const long from_step, const long to_step )
{
  // time zero: slice_origin.get_steps() == 0 and from_step == 0
  
  for ( long lag = from_step; lag < to_step; ++lag )
  {
    // order is important in this loop, since we have to use the old values
    // (those upon entry to the loop) on right hand sides everywhere

    // log membrane potential
    B_.logger_.record_data( slice_origin.get_steps() + lag );
  }
}

void
ctc::ctc_loss::handle( LearningSignalConnectionEvent& e )
{
   for ( auto it_event = e.begin(); it_event != e.end(); )
  {
    const long channel = e.get_rport();
    const double prediction = e.get_coeffvalue( it_event ); // get_coeffvalue advances iterator

    assert( 0 < channel and channel <= P_.num_input_channels );
    B_.predictions.at( channel - 1 ).add_value(
					       e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), prediction );
  }

}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void
ctc::ctc_loss::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e ); // the logger does this for us
}
