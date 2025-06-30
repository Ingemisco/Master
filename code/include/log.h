#ifndef INCLUDE_INCLUDE_LOG_H_
#define INCLUDE_INCLUDE_LOG_H_

#include "datastructures.h"
#include "simplification.h"
#include <chrono>
#include <string>
#include <vector>

namespace Log {

extern std::string measurement_directory;

enum class Algorithm {
  SIMPLIFICATION_SIMPLE_MANHATTAN,
  SIMPLIFICATION_SIMPLE_EUCLIDEAN,
  SIMPLIFICATION_SIMPLE_IMPLICIT_EUCLIDEAN,
  SIMPLIFICATION_SIMPLE_SEMIEXPLICIT_EUCLIDEAN,
  SIMPLIFICATION_SIMPLE_CHEBYSHEV,
  SIMPLIFICATION_SIMPLE_IMPLICIT_MINKOWSKI,

  SIMPLIFICATION_ADVANCED_MANHATTAN,
  SIMPLIFICATION_ADVANCED_EUCLIDEAN,
  SIMPLIFICATION_ADVANCED_IMPLICIT_EUCLIDEAN,
  SIMPLIFICATION_ADVANCED_SEMIEXPLICIT_EUCLIDEAN,
  SIMPLIFICATION_ADVANCED_CHEBYSHEV,
  SIMPLIFICATION_ADVANCED_IMPLICIT_MINKOWSKI,
};

typedef std::chrono::duration<double> Time;

class PerformanceLogger final {
private:
  std::string const header;
  Algorithm const algorithm;
  std::vector<Time> times;
  std::vector<size_t> simplification_sizes;
  std::vector<std::string> case_names;
  size_t dimension = 0;
  size_t point_count = 0;

public:
  PerformanceLogger(Algorithm, std::string);

  void add_data(DataStructures::Polyline &, Simplification::Simplification &,
                std::chrono::duration<double>, std::string);
  void emit();
};

} // namespace Log
#endif // INCLUDE_INCLUDE_LOG_H_
