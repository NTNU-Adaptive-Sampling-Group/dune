if(LITE)
  if(Navigation.General.ROV)
    set(TASK_ENABLED TRUE)
  else(Navigation.General.ROV)
    set(TASK_ENABLED FALSE)
  endif(Navigation.General.ROV)
endif(LITE)
