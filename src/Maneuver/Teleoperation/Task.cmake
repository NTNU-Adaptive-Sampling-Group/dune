if(LITE)
  if(Maneuver.Teleoperation)
    set(TASK_ENABLED TRUE)
  else(Maneuver.Teleoperation)
    set(TASK_ENABLED FALSE)
  endif(Maneuver.Teleoperation)
endif(LITE)
