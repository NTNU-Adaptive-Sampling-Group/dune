if(LITE)
  if(Monitors.Emergency)
    set(TASK_ENABLED TRUE)
  else(Monitors.Emergency)
    set(TASK_ENABLED FALSE)
  endif(Monitors.Emergency)
endif(LITE)
