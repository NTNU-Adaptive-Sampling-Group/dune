if(LITE)
  if(Maneuver.CommsRelay)
    set(TASK_ENABLED TRUE)
  else(Maneuver.CommsRelay)
    set(TASK_ENABLED FALSE)
  endif(Maneuver.CommsRelay)
endif(LITE)
