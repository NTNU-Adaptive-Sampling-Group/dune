if(LITE)
  if(EmulatedGPS)
    set(TASK_ENABLED TRUE)
  else(EmulatedGPS)
    set(TASK_ENABLED FALSE)
  endif(EmulatedGPS)
endif(LITE)
