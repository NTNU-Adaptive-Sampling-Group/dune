if(LITE)
  if(Control.ROV.Depth)
    set(TASK_ENABLED TRUE)
  else(Control.ROV.Depth)
    set(TASK_ENABLED FALSE)
  endif(Control.ROV.Depth)
endif(LITE)
