if(LITE)
  if(Control.UAV.LOS)
    set(TASK_ENABLED TRUE)
  else(Control.UAV.LOS)
    set(TASK_ENABLED FALSE)
  endif(Control.UAV.LOS)
endif(LITE)
