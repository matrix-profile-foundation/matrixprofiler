#include <Rcpp.h>
#include <RcppParallel.h>

#if RCPP_PARALLEL_USE_TBB // TBB support turned on

#include <string>
#include <exception>
#include <tbb/task_scheduler_init.h>

using namespace Rcpp;

void setThreadOptions(int numThreads, int stackSize) {

  static tbb::task_scheduler_init *s_pTaskScheduler = NULL;

  try {
    if (!s_pTaskScheduler) {
      s_pTaskScheduler = new tbb::task_scheduler_init(numThreads, stackSize);
    } else {
      s_pTaskScheduler->terminate();
      s_pTaskScheduler->initialize(numThreads, stackSize);
    }
  } catch (const std::exception &e) {
    std::cout << "Error loading TBB: " << e.what() << std::endl;
  } catch (...) {
    std::cout << "Error loading TBB: (Unknown error)" << std::endl;
  }

  return;
}

int defaultNumThreads() {
  int value = tbb::task_scheduler_init::default_num_threads();
  return value;
}

#else // TBB support not turned on

#include <tthread/tinythread.h>

void setThreadOptions(int numThreads, int stackSize) {
  return;
}

int defaultNumThreads() {
  int value = tthread::thread::hardware_concurrency();
  return value;
}

#endif
