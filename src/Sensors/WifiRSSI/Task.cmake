if(LITE)
  if(WifiRSSI)
    set(TASK_ENABLED TRUE)
  else(WifiRSSI)
    set(TASK_ENABLED FALSE)
  endif(WifiRSSI)
endif(LITE)
