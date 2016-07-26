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
// Author: PGonçalves                                                       *
//***************************************************************************

// ISO C++ 98 headers.
#include <queue>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string>

// DUNE headers.
#include <DUNE/DUNE.hpp>

//OpenCV headers
#include <opencv2/opencv.hpp>

// Local headers.
#include "IpCamCap.hpp"
#include "OperationCV.hpp"
#include "StereoMatch.hpp"

namespace Vision
{
  namespace Tracking3d
  {
    using DUNE_NAMESPACES;
    using namespace cv;

    static const unsigned int c_sleep_time = 1500;
    static const unsigned int c_number_of_pre_frames = 40;
    static const unsigned int c_max_cam = 2;

    struct Arguments
    {
      //! IpCam1 URL
      std::string url_ipcam1;
      //! IpCam2 URL
      std::string url_ipcam2;
      //! IpCam Name
      std::string name_ipcam[c_max_cam];
      //! Intrinsic values of IPCam 1.
      std::vector<double> intrinsicCam1;
      //! Intrinsic values of IPCam 2.
      std::vector<double> intrinsicCam2;
      //! Distortion values of IPCam 1.
      std::vector<double> distortionCam1;
      //! Distortion values of IPCam 2.
      std::vector<double> distortionCam2;
      //! Position of marks in Pixels (X) CAM 1.
      std::vector<double> positionPixelsX1;
      //! Position of marks in Pixels (Y) CAM 1.
      std::vector<double> positionPixelsY1;
      //! Position of marks in Pixels (X) CAM 2.
      std::vector<double> positionPixelsX2;
      //! Position of marks in Pixels (Y) CAM 2.
      std::vector<double> positionPixelsY2;
      //! Position of marks in meters (X).
      std::vector<double> positionMetersX;
      //! Position of marks in meters (Y).
      std::vector<double> positionMetersY;
      //! Position of marks in meters (Z).
      std::vector<double> positionMetersZ;
      //! OffSet of Y axis in meters
      float offSetY;
      //! Size of Tpl match
      int tpl_size;
      //! Windows search size
      int window_search_size;
      //! Frames to Refresh
      int frames_to_refresh;
    };

    struct Task : public DUNE::Tasks::Task
    {
      //! Temperature messages.
      Arguments m_args;
      IMC::GetImageCoords m_getImgCoord[c_max_cam];
      IMC::GetWorldCoordinates m_getworldcoord;
      IpCamCap* m_cap1;
      IpCamCap* m_cap2;
      OperationCV* m_operation1;
      OperationCV* m_operation2;
      StereoMatch* m_stereo_match;
      //! buffer for frame of Cam1
      IplImage* m_frameCam1;
      //! buffer for frame of Cam2
      IplImage* m_frameCam2;
      //! Init tpl values
      bool m_initValuesTpl;
      //! State of tracking Cam1
      bool m_isTrackingCam1;
      //! State of tracking Cam2
      bool m_isTrackingCam2;
      //! Read timestamp.
      double m_tstamp;

      //! Constructor.
      //! @param[in] name task name.
      //! @param[in] ctx context.
      Task(const std::string& name, Tasks::Context& ctx) :
        DUNE::Tasks::Task(name, ctx), m_cap1(NULL), m_cap2(NULL),
        m_tstamp(0)
      {
        param("IpCam1 - URL", m_args.url_ipcam1)
        .description("IpCam1 Addresses");

        param("IpCam2 - URL", m_args.url_ipcam2)
        .description("IpCam2 Addresses");

        param("IpCam1 - Name", m_args.name_ipcam[0])
        .defaultValue("Cam1")
        .description("IpCam1 Name");

        param("IpCam2 - Name", m_args.name_ipcam[1])
        .defaultValue("Cma2")
        .description("IpCam2 Name");

        param("IpCam1 - Intrinsic Matrix", m_args.intrinsicCam1)
        .description("Intrinsic values of IPCam 1");

        param("IpCam2 - Intrinsic Matrix", m_args.intrinsicCam2)
        .description("Intrinsic values of IPCam 2");

        param("IpCam1 - Distortion Vector", m_args.distortionCam1)
        .description("Distortion values of IPCam 1");

        param("IpCam2 - Distortion Vector", m_args.distortionCam2)
        .description("Distortion values of IPCam 2");

        param("Position IpCam1 - Pixels X", m_args.positionPixelsX1)
        .description("Position of marks in Pixels");

        param("Position IpCam1 - Pixels Y", m_args.positionPixelsY1)
        .description("Position of marks in Pixels");

        param("Position IpCam2 - Pixels X", m_args.positionPixelsX2)
        .description("Position of marks in Pixels");

        param("Position IpCam2 - Pixels Y", m_args.positionPixelsY2)
        .description("Position of marks in Pixels");

        param("Position - Meters X", m_args.positionMetersX)
        .description("Position of marks in Meters");

        param("Position - Meters Y", m_args.positionMetersY)
        .description("Position of marks in Meters");

        param("Position - Meters Z", m_args.positionMetersZ)
        .description("Position of marks in Meters");

        param("OffSet Y", m_args.offSetY)
        .description("OffSet of Y in Meters");

        param("Tpl Size", m_args.tpl_size)
        .defaultValue("50")
        .description("Size of TPL match");

        param("Window Search Size", m_args.window_search_size)
        .defaultValue("160")
        .description("Size of Window Search Size");

        param("Frames to Refresh", m_args.frames_to_refresh)
        .defaultValue("30")
        .description("Number of frames necessary to auto refresh TPL");

        bind<IMC::SetImageCoords>(this);
      }

      //! Reserve entity identifiers.
      void
      onEntityReservation(void)
      {
        unsigned eid = 0;

        for (unsigned i = 0; i < c_max_cam; ++i)
          {
            if (m_args.name_ipcam[i].empty())
            continue;

            try
            {
              eid = resolveEntity(m_args.name_ipcam[i]);
            }
            catch (Entities::EntityDataBase::NonexistentLabel& e)
            {
              (void)e;
              eid = reserveEntity(m_args.name_ipcam[i]);
            }

            m_getImgCoord[i].setSourceEntity(eid);
          }
        }

      //! Initialize resources.
      void
      onResourceInitialization(void)
      {
        m_cap1 = new IpCamCap(this, m_args.url_ipcam1);
        m_cap2 = new IpCamCap(this, m_args.url_ipcam2);
        m_operation1 = new OperationCV(this, m_args.url_ipcam1, m_args.tpl_size, m_args.window_search_size, m_args.frames_to_refresh);
        m_operation2 = new OperationCV(this, m_args.url_ipcam2, m_args.tpl_size, m_args.window_search_size, m_args.frames_to_refresh);
        setEntityState(IMC::EntityState::ESTA_NORMAL, Status::CODE_IDLE);
        m_initValuesTpl = false;
        m_isTrackingCam1 = false;
        m_stereo_match = new StereoMatch(this);
        inf("Running stereo calibration");
        m_stereo_match->loadParametersForStereo(m_args.intrinsicCam1, m_args.distortionCam1, m_args.intrinsicCam2,
                                                m_args.distortionCam2, m_args.positionPixelsX1, m_args.positionPixelsY1,
                                                m_args.positionPixelsX2, m_args.positionPixelsY2, m_args.positionMetersX,
                                                m_args.positionMetersY, m_args.positionMetersZ);
        m_cap1->start();
        Delay::waitMsec(c_sleep_time);
        m_cap2->start();
        preLoadFrame(c_number_of_pre_frames);
      }

      //! Release resources.
      void
      onResourceRelease(void)
      {
        if (m_cap1 != NULL)
        {
          m_cap1->stopAndJoin();
          delete m_cap1;
          m_cap1 = NULL;
        }
        if (m_cap2 != NULL)
        {
          m_cap2->stopAndJoin();
          delete m_cap2;
          m_cap2 = NULL;
        }
      }

      //! Consume Message
      void
      consume(const IMC::SetImageCoords* msg)
      {
        if (msg->camid == 1)
        {
          if (m_frameCam1 != NULL)
          {
            if (m_initValuesTpl)
              m_operation1->setNewTPL(msg->x, msg->y, m_frameCam1, m_args.name_ipcam[0]);
          }
        }
        else
        {
          if (m_frameCam2 != NULL)
          {
            if (m_initValuesTpl)
              m_operation2->setNewTPL(msg->x, msg->y, m_frameCam2, m_args.name_ipcam[1]);
          }
        }
      }

      void
      preLoadFrame(int ntimes)
      {
        int t = 0;
        while (t < ntimes && !stopping())
        {
          m_frameCam1 = m_cap1->capFrame();
          m_frameCam2 = m_cap2->capFrame();
          if (!m_cap1->isConnected() || !m_cap2->isConnected())
            setEntityState(IMC::EntityState::ESTA_ERROR, Status::CODE_COM_ERROR);
          if (m_frameCam1 != NULL && m_frameCam2 != NULL)
            t++;
        }
        if (!stopping())
        {
          m_operation1->inicTplTest(m_frameCam1);
          m_operation2->inicTplTest(m_frameCam2);
          m_initValuesTpl = true;
        }
      }

      void
      getPosition(void)
      {
        m_tstamp = Clock::getSinceEpoch();
        m_frameCam1 = m_cap1->capFrame();
        m_frameCam2 = m_cap2->capFrame();

        if (m_frameCam1 != NULL && m_frameCam2 != NULL)
        {
          setEntityState(IMC::EntityState::ESTA_NORMAL, Status::CODE_ACTIVE);
          m_isTrackingCam1 = m_operation1->trackObject(m_frameCam1, m_args.name_ipcam[0]);
          m_isTrackingCam2 = m_operation2->trackObject(m_frameCam2, m_args.name_ipcam[1]);

          for (unsigned int i = 0; i < c_max_cam; i++)
          {
            m_getImgCoord[i].setTimeStamp(m_tstamp);
            m_getImgCoord[i].camid = i + 1;
            if (i == 0)
            {
              m_getImgCoord[i].x = m_operation1->m_coordImage.x;
              m_getImgCoord[i].y = m_operation1->m_coordImage.y;
            }
            else if (i == 1)
            {
              m_getImgCoord[i].x = m_operation2->m_coordImage.x;
              m_getImgCoord[i].y = m_operation2->m_coordImage.y;
            }

            dispatch(m_getImgCoord[i], DF_KEEP_TIME);
          }

          if (m_isTrackingCam1 && m_isTrackingCam2)
          {
            m_getworldcoord.setSourceEntity(getEntityId());

            if (m_stereo_match->getRealCoord(m_operation1->m_coordImage.x, m_operation1->m_coordImage.y, m_operation2->m_coordImage.x, m_operation2->m_coordImage.y))
            {
              m_getworldcoord.x = m_stereo_match->m_real_coord.x;
              m_getworldcoord.y = m_stereo_match->m_real_coord.y - m_args.offSetY;
              m_getworldcoord.z = m_stereo_match->m_real_coord.z;
            }
            else
            {
              m_getworldcoord.x = 0;
              m_getworldcoord.y = 0;
              m_getworldcoord.z = 0;
              m_getworldcoord.tracking = false;
            }
              m_getworldcoord.tracking = true;
              dispatch(m_getworldcoord);
          }
          else
          {
            m_getworldcoord.setSourceEntity(getEntityId());
            m_getworldcoord.tracking = false;
            dispatch(m_getworldcoord);
          }
        }
        else
        {
          setEntityState(IMC::EntityState::ESTA_ERROR, Status::CODE_IO_ERROR);
          throw RestartNeeded(DTR("null frame"), 1, true);
        }
      }

      //! Main loop.
      void
      onMain(void)
      {
        while (!stopping())
        {
          getPosition();
          waitForMessages(0.01);
          m_frameCam1 = NULL;
          m_frameCam2 = NULL;
        }
      }
    };
  }
}

DUNE_TASK
