if(LITE)
  if(Sensors.SW100)
    set(TASK_ENABLED TRUE)
  else(Sensors.SW100)
    set(TASK_ENABLED FALSE)
  endif(Sensors.SW100)
endif(LITE)
