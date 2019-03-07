#pragma once

#if defined(SEIDR_PSTL)
#include <tbb/task_scheduler_init.h>
#endif

void set_pstl_threads(int target);