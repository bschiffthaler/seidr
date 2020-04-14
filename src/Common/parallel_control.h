#pragma once

#if defined(SEIDR_PSTL)
#include <tbb/global_control.h>
#include <tbb/task_scheduler_init.h>

void
set_pstl_threads(int target);

#else

void
set_pstl_threads(int target);

#endif