#ifndef INCLUDE_INCLUDE_LOG_H_
#define INCLUDE_INCLUDE_LOG_H_

#include "datastructures.h"
#include <chrono>
#include <string>
#include <vector>

namespace Log {

enum class Algorithm {
  SIMPLIFICATION_SIMPLE_MANHATTAN,
  SIMPLIFICATION_SIMPLE_EUCLIDEAN,
  SIMPLIFICATION_SIMPLE_IMPLICIT_EUCLIDEAN,
  SIMPLIFICATION_SIMPLE_CHEBYSHEV,
  SIMPLIFICATION_SIMPLE_IMPLICIT_MINKOWSKI,

  SIMPLIFICATION_ADVANCED_MANHATTAN,
  SIMPLIFICATION_ADVANCED_EUCLIDEAN,
  SIMPLIFICATION_ADVANCED_IMPLICIT_EUCLIDEAN,
  SIMPLIFICATION_ADVANCED_CHEBYSHEV,
  SIMPLIFICATION_ADVANCED_IMPLICIT_MINKOWSKI,
};

class PerformanceLogger final {
private:
  std::string const header;
  Algorithm const alogrithm;
  std::chrono::duration<double> total_time;
  std::chrono::duration<double> min_time;
  std::chrono::duration<double> max_time;
  size_t data_count;
  size_t min_dim;
  size_t max_dim;
  size_t min_points;
  size_t max_points;

public:
  PerformanceLogger(Algorithm, std::string);

  void add_data(DataStructures::Polyline &, std::chrono::duration<double>);
  void emit();
};

} // namespace Log
#endif // INCLUDE_INCLUDE_LOG_H_
