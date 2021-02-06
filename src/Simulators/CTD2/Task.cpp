//***************************************************************************
// Copyright 2007-2016 Universidade do Porto - Faculdade de Engenharia      *
// Laboratório de Sistemas e Tecnologia Subaquática (LSTS)                  *
//***************************************************************************
// This file is part of DUNE: Unified Navigation Environment.               *
//                                                                          *
// Commercial Licence Usage                                                 *
// Licencees holding valid commercial DUNE licences may use this file in    *
// accordance with the commercial licence agreement provided with the       *
// Software or, alternatively, in accordance with the terms contained in a  *
// written agreement between you and Universidade do Porto. For licensing   *
// terms, conditions, and further information contact lsts@fe.up.pt.        *
//                                                                          *
// European Union Public Licence - EUPL v.1.1 Usage                         *
// Alternatively, this file may be used under the terms of the EUPL,        *
// Version 1.1 only (the "Licence"), appearing in the file LICENCE.md       *
// included in the packaging of this file. You may not use this work        *
// except in compliance with the Licence. Unless required by applicable     *
// law or agreed to in writing, software distributed under the Licence is   *
// distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF     *
// ANY KIND, either express or implied. See the Licence for the specific    *
// language governing permissions and limitations at                        *
// http://ec.europa.eu/idabc/eupl.html.                                     *
//***************************************************************************
// Author: Pedro Calado                                                     *
//***************************************************************************

// DUNE headers.
#include <DUNE/DUNE.hpp>
#include <math.h>

using DUNE_NAMESPACES;

namespace Simulators
{
  namespace CTD2
  {
    //! %Task arguments.
    struct Arguments
    {
      //! Standard deviation of temperature measurements.
      double std_dev_temp;
      //! Mean temperature value.
      float mean_temp;
      //! Standard deviation of conductivity measurements.
      double std_dev_cond;
      //! Mean conductivity value.
      float mean_cond;
      //! Standard deviation of depth measurements.
      double std_dev_depth;
      //! PRNG type.
      std::string prng_type;
      //! PRNG seed.
      int prng_seed;
      //! Time of year profile
      //std::string season;
      //! X-value of estimated state
      double x_pos;
      //! Y-value of estimated state
      double y_pos;
      //! y = ax^2 + bx + c
      double a_val;
      double b_val;
      double c_val;
      //inside outside plume
      bool is_inside;

      //! Width of fuzzy transition
      double fuzzy_trans;
      //! Temperature difference
      double temp_diff;
      //! Conductivity difference
      double cond_diff;
      //! fuzzy gradient
      float gradient_f;
      //! distance from front
      float dist_from_front;
      //! Sigmoid solution
      float sig_sol;
      //! Position noise
      float pos_noise_std_dev;

    };

    //! %SVS simulator task.
    struct Task: public Tasks::Periodic
    {
      //! Temperature.
      IMC::Temperature m_temp;
      //! Current sound speed.
      IMC::SoundSpeed m_sspeed;
      IMC::Conductivity m_cond;
      IMC::Salinity m_salinity;
      IMC::Depth m_depth;
      //! Simulated state.
      IMC::SimulatedState m_sstate;
      //! PRNG handle.
      Random::Generator* m_prng;
      //! Task arguments.
      Arguments m_args;

      Task(const std::string& name, Tasks::Context& ctx):
        Tasks::Periodic(name, ctx),
        m_prng(nullptr)
      {
        // Retrieve configuration values.
        param("Standard Deviation - Temperature", m_args.std_dev_temp)
        .defaultValue("0.02");

        param("Mean Value - Temperature", m_args.mean_temp)
        .defaultValue("14.0");

        param("Standard Deviation - Conductivity", m_args.std_dev_cond)
        .defaultValue("1.0");

        param("Mean Value - Conductivity", m_args.mean_cond)
        .defaultValue("3.0");

        param("Standard Deviation - Depth", m_args.std_dev_depth)
        .defaultValue("0.1");

        param("PRNG Type", m_args.prng_type)
        .defaultValue(Random::Factory::c_default);

        param("PRNG Seed", m_args.prng_seed)
        .defaultValue("-1");

        param("a",m_args.a_val)
        .defaultValue("0.0");

        param("b",m_args.b_val)
        .defaultValue("0.0");

        param("c",m_args.c_val)
        .defaultValue("0.0");

        param("is the vehicle starting inside the plume?",m_args.is_inside)
        .defaultValue("true");

        param("Width of fuzzy transition", m_args.fuzzy_trans)
        .defaultValue("10.0");

        param("temperature diffrence", m_args.temp_diff)
        .defaultValue("-2.1");

        param("Conductivity diffrence", m_args.cond_diff)
        .defaultValue("0.2");

        param("Position noise", m_args.pos_noise_std_dev)
        .defaultValue("2.5");
        //param("CTD season", m_args.season)
        //.defaultValue("spring")

        // Register consumers.
        bind<IMC::SimulatedState>(this);
      }

      //! Initialize resources.
      void
      onResourceInitialization() override
      {
        requestDeactivation();
      }

      //! Acquire resources.
      void
      onResourceAcquisition() override
      {
        m_prng = Random::Factory::create(m_args.prng_type, m_args.prng_seed);
      }

      //! Release resources.
      void
      onResourceRelease() override
      {
        Memory::clear(m_prng);
      }

      void
      consume(const IMC::SimulatedState* msg)
      {
        if (!isActive())
        {
          setEntityState(IMC::EntityState::ESTA_NORMAL, Status::CODE_ACTIVE);
          requestActivation();
        }

        m_args.x_pos = (msg->x)+m_args.pos_noise_std_dev*m_prng->gaussian();
        m_args.y_pos = (msg->y)+m_args.pos_noise_std_dev*m_prng->gaussian();

/*
        if (m_args.y_pos >= m_args.a_val*m_args.x_pos*m_args.x_pos + m_args.b_val*m_args.x_pos + m_args.c_val) {
          m_args.is_inside = true;
        }else{
          m_args.is_inside = false;
        }
*/
        if(m_args.fuzzy_trans >  abs(m_args.y_pos - m_args.a_val*m_args.x_pos*m_args.x_pos + m_args.b_val*m_args.x_pos + m_args.c_val)){
          m_args.gradient_f = std::max(((m_args.y_pos - m_args.a_val*m_args.x_pos*m_args.x_pos + m_args.b_val*m_args.x_pos + m_args.c_val) / m_args.fuzzy_trans),0.0);
        }else{
          m_args.gradient_f = 1.0;
        }
        m_args.gradient_f= m_args.y_pos - (m_args.a_val*m_args.x_pos*m_args.x_pos + m_args.b_val*m_args.x_pos + m_args.c_val);
        m_args.sig_sol = (exp(m_args.gradient_f/m_args.fuzzy_trans)/(exp(m_args.gradient_f/m_args.fuzzy_trans)+1));

        m_sstate = *msg;
      }

      void
      task() override
      {
        // Return if task is not active.
        if (!isActive())
          return;

        m_temp.setTimeStamp();

        m_depth.setTimeStamp(m_temp.getTimeStamp());
        m_depth.value = std::max(m_sstate.z + m_prng->gaussian() * m_args.std_dev_depth, 0.0);

        m_cond.setTimeStamp(m_temp.getTimeStamp());
        m_cond.value = (m_args.mean_cond + m_prng->gaussian() * m_args.std_dev_cond) + 0.1*log10(m_depth.value+10.0) + m_args.cond_diff*m_args.sig_sol;
        //std::cout << "\n\n" << m_cond.value << std::endl;

        m_temp.value = m_args.mean_temp + (m_prng->gaussian() * m_args.std_dev_temp) + (m_depth.value*10.0*exp(-m_depth.value*0.1))*0.05 + m_args.temp_diff*m_args.sig_sol;
        //std::cout << m_temp.value << std::endl;
        //std::cout << m_depth.value << std::endl;
        //std::cout << m_args.sig_sol << std::endl;

        // Compute pressure.
        double pressure = (m_depth.value * c_gravity * c_seawater_density + c_sea_level_pressure) / c_pascal_per_bar;

        m_salinity.setTimeStamp(m_temp.getTimeStamp());
        m_salinity.value = UNESCO1983::computeSalinity(m_cond.value, pressure, m_temp.value);

        m_sspeed.setTimeStamp(m_temp.getTimeStamp());
        m_sspeed.value = (m_salinity.value < 0.0) ? -1.0 : UNESCO1983::computeSoundSpeed(m_salinity.value, pressure, m_temp.value);

        dispatch(m_temp, DF_KEEP_TIME);
        dispatch(m_cond, DF_KEEP_TIME);
        dispatch(m_depth, DF_KEEP_TIME);
        dispatch(m_salinity, DF_KEEP_TIME);
        dispatch(m_sspeed, DF_KEEP_TIME);
      }
    };
  }
}

DUNE_TASK
