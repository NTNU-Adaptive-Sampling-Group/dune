if(LITE)
  if(TREX)
    set(TASK_ENABLED TRUE)
    # This task is for Linux and Mac OS X only.
    if(NOT DUNE_OS_LINUX AND NOT DUNE_OS_DARWIN)
      set(TASK_ENABLED FALSE)
    endif(NOT DUNE_OS_LINUX AND NOT DUNE_OS_DARWIN)
  else(TREX)
    set(TASK_ENABLED FALSE)
  endif(TREX)
else(LITE)
  # This task is for Linux and Mac OS X only.
  if(NOT DUNE_OS_LINUX AND NOT DUNE_OS_DARWIN)
    set(TASK_ENABLED FALSE)
  endif(NOT DUNE_OS_LINUX AND NOT DUNE_OS_DARWIN)
endif(LITE)
