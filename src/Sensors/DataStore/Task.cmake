if(LITE)
  if(Sensors.DataStore)
    set(TASK_ENABLED TRUE)
  else(Sensors.DataStore)
    set(TASK_ENABLED FALSE)
  endif(Sensors.DataStore)
endif(LITE)
