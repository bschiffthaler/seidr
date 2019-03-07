#include <common.h>
#include <parallel_control.h>
#include <BSlogger.h>

#if defined(SEIDR_PSTL)
#include <tbb/task_scheduler_init.h>
#include <omp.h>
#endif

void set_pstl_threads(int target)
{
  logger log(std::cerr, "set_pstl_threads");
#ifdef SEIDR_PSTL
  assert_in_range<int>(target, 1, GET_MAX_PSTL_THREADS(),
                       "--threads");
  SET_NUM_PSTL_THREADS(target);
  omp_set_num_threads(target);
#else
  if (target > 1)
  {
    log(LOG_WARN) << "More than one thread selected, but seidr wasn't "
                  << "compiled for parallel sorting. "
                  << "Flag will have no effect.\n";
  }
#endif
}