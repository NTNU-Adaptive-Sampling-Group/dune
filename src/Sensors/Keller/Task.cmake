if(LITE)
  if(Sensors.Keller)
    set(TASK_ENABLED TRUE)
  else(Sensors.Keller)
    set(TASK_ENABLED FALSE)
  endif(Sensors.Keller)
endif(LITE)
