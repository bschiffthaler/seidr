#include <BSlogger.hpp>
#include <common.h>
#include <parallel_control.h>

#if defined(SEIDR_PSTL)
#include <omp.h>
#include <tbb/global_control.h>
#include <tbb/task_scheduler_init.h>

static tbb::global_control global_control(
  tbb::global_control::max_allowed_parallelism,
  tbb::task_scheduler_init::default_num_threads());

void
set_pstl_threads(int target)
{
  logger log(std::cerr, "set_pstl_threads");
#ifdef SEIDR_PSTL
  assert_in_range<int>(target, 1, GET_MAX_PSTL_THREADS(), "--threads");
  tbb::task_scheduler_init init(target);
  omp_set_num_threads(target);
#else
  if (target > 1) {
    log(LOG_WARN) << "More than one thread selected, but seidr wasn't "
                  << "compiled for parallel sorting. "
                  << "Flag will have no effect.\n";
  }
#endif
}

#else
void
set_pstl_threads(int target) // NOLINT(clang-diagnostic-unused-parameter)
{}
#endif