if(LITE)
  if(Sensors.MLBL)
    set(TASK_ENABLED TRUE)
  else(Sensors.MLBL)
    set(TASK_ENABLED FALSE)
  endif(Sensors.MLBL)
endif(LITE)
