if(LITE)
  if(Control.AUV.Speed)
    set(TASK_ENABLED TRUE)
  else(Control.AUV.Speed)
    set(TASK_ENABLED FALSE)
  endif(Control.AUV.Speed)
endif(LITE)
