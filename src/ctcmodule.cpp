/*
 *  ctcmodule.cpp
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

// include headers with your own stuff
#include "ctc_loss.h"

// Includes from NEST
#include "nest_extension_interface.h"

namespace ctc
{
  class CtcModule : public nest::NESTExtensionInterface
  {
  public:
    CtcModule() {}
    virtual ~CtcModule() {}

    void initialize() override;
  };
}

// Define module instance outside of namespace to avoid name-mangling problems
ctc::CtcModule ctcmodule_LTX_module;

void ctc::CtcModule::initialize()
{
  /* Register a neuron or device model.
   */
  ctc::register_ctc_loss( "ctc_loss" );
}

