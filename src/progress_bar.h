// This is a personal academic project. Dear PVS-Studio, please check it.

#ifndef MATRIXPROFILER_PROGRESS_BAR_H
#define MATRIXPROFILER_PROGRESS_BAR_H

#include <RcppThread.h>
#include <memory>
#include <cstddef>

namespace matrixprofiler {

class ProgressBar {
private:
  std::unique_ptr<RcppThread::ProgressBar> bar_;

public:
  explicit ProgressBar(std::size_t total, bool show = true, std::size_t print_every = 1) {
    if (show && total > 0) {
      bar_ = std::unique_ptr<RcppThread::ProgressBar>(
        new RcppThread::ProgressBar(total, print_every)
      );
    }
  }

  void increment(std::size_t amount = 1) {
    if (bar_) {
      for (std::size_t i = 0; i < amount; ++i) {
        (*bar_)++;
      }
    }
  }

  std::size_t operator++() {
    if (bar_) {
      return ++(*bar_);
    }
    return 0;
  }

  std::size_t operator++(int) {
    if (bar_) {
      return (*bar_)++;
    }
    return 0;
  }
};

} // namespace matrixprofiler

#endif // MATRIXPROFILER_PROGRESS_BAR_H
