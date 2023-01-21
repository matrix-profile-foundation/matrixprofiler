#ifndef __RCPP_PARALLEL_TINYTHREAD_FIX__
#define __RCPP_PARALLEL_TINYTHREAD_FIX__

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

namespace RcppParallel2 {
namespace { // NOLINT: anonymous namespace
// Function to calculate the ranges for a given input
inline std::vector<IndexRange> splitInputRange(const IndexRange &range, std::size_t grainSize) {

  // determine max number of threads
  std::size_t threads = tthread::thread::hardware_concurrency();
  char *numThreads = ::getenv("RCPP_PARALLEL_NUM_THREADS"); // NOLINT(concurrency-mt-unsafe)
  if (numThreads != nullptr) {
    long const parsedThreads = ::strtol(numThreads, nullptr, 10);
    if (parsedThreads > 0) {
      threads = parsedThreads;
    }
  }

  // compute grainSize (including enforcing requested minimum)
  std::size_t const length = range.end() - range.begin();
  if (threads == 1) {
    grainSize = length;
  } else if ((length % threads) == 0) { // perfect division
    grainSize = std::max(length / threads, grainSize);
  } else {                              // imperfect division, divide by threads - 1
    grainSize = std::max(length / (threads - 1), grainSize);
  }
  // allocate ranges
  std::vector<IndexRange> ranges;
  std::size_t begin = range.begin();
  std::size_t end = begin;
  while (begin < range.end()) {
    if ((range.end() - (begin + grainSize)) < grainSize) {
      end = range.end();
    } else {
      end = std::min(begin + grainSize, range.end());
    }

    ranges.emplace_back(begin, end);
    begin = end;
  }

  // return ranges
  return ranges;
}

// Execute the Worker over the IndexRange in parallel
inline void ttParallelFor(std::size_t begin, std::size_t end, Worker &worker, std::size_t grainSize = 1) {
  // split the work
  IndexRange const inputRange(begin, end);
  std::vector<IndexRange> const ranges = RcppParallel2::splitInputRange(inputRange, grainSize);

  // create threads
  std::vector<tthread::thread *> threads;
  threads.reserve(ranges.size());

  for (IndexRange const range : ranges) {
    threads.push_back(new tthread::thread(workerThread, new Work(range, worker)));
  }

  // join and delete them
  for (tthread::thread *thread : threads) {
    thread->join();
    delete thread;
  }
}

} // anonymous namespace
} // namespace RcppParallel2

#endif // __RCPP_PARALLEL_TINYTHREAD_FIX__
