if(LITE)
  if(Sensors.Imagenex852)
    set(TASK_ENABLED TRUE)
  else(Sensors.Imagenex852)
    set(TASK_ENABLED FALSE)
  endif(Sensors.Imagenex852)
endif(LITE)
