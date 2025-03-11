#include "datastructures.h"
#include "log.h"
#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace Log {
PerformanceLogger::PerformanceLogger(Algorithm algorithm, std::string header)
    : header(header), alogrithm(algorithm), total_time(0), min_time(DBL_MAX),
      max_time(0) {}

void PerformanceLogger::add_data(DataStructures::Polyline &polyline,
                                 std::chrono::duration<double> elapsed_time) {

  this->total_time += elapsed_time;
  this->data_count++;

  this->max_dim = std::max(this->max_dim, polyline.dimension);
  this->min_dim = std::min(this->min_dim, polyline.dimension);
  this->max_points = std::max(this->max_points, polyline.point_count);
  this->min_points = std::min(this->min_points, polyline.point_count);
  this->min_time = std::min(this->min_time, elapsed_time);
  this->max_time = std::max(this->max_time, elapsed_time);
}

static std::string algorithm_name(Algorithm algorithm) {
  switch (algorithm) {
  case Algorithm::SIMPLIFICATION_SIMPLE_MANHATTAN:
    return "Manhattan Simplification (Simple)";
  case Algorithm::SIMPLIFICATION_SIMPLE_EUCLIDEAN:
    return "Euclidean Simplification (Simple)";
  case Algorithm::SIMPLIFICATION_SIMPLE_IMPLICIT_EUCLIDEAN:
    return "Euclidean Simplification (Simple, Implicit)";
  case Algorithm::SIMPLIFICATION_SIMPLE_CHEBYSHEV:
    return "Chebyshev Simplification (Simple)";
  case Algorithm::SIMPLIFICATION_SIMPLE_IMPLICIT_MINKOWSKI:
    return "Minkowski Simplification (Simple, Implicit)";

  case Algorithm::SIMPLIFICATION_ADVANCED_MANHATTAN:
    return "Manhattan Simplification (Advanced)";
  case Algorithm::SIMPLIFICATION_ADVANCED_EUCLIDEAN:
    return "Euclidean Simplification (Advanced)";
  case Algorithm::SIMPLIFICATION_ADVANCED_IMPLICIT_EUCLIDEAN:
    return "Euclidean Simplification (Advanced, Implicit)";
  case Algorithm::SIMPLIFICATION_ADVANCED_CHEBYSHEV:
    return "Chebyshev Simplification (Advanced)";
  case Algorithm::SIMPLIFICATION_ADVANCED_IMPLICIT_MINKOWSKI:
    return "Minkowski Simplification (Advanced, Implicit)";
  }
}

void PerformanceLogger::emit() {
  if (this->data_count == 0) {
    throw std::runtime_error("No data has been entered.");
  }
  if (!std::filesystem::exists("measurements")) {
    std::filesystem::create_directory("measurements");
  }
  if (!std::filesystem::is_directory("measurements")) {
    throw std::runtime_error("File called 'measurements' exists. Cannot create "
                             "a directory of same name.");
  }

  // Get the current time
  auto now = std::chrono::system_clock::now();
  auto now_time_t = std::chrono::system_clock::to_time_t(now);
  auto now_tm = *std::localtime(&now_time_t);

  // Format the timestamp
  std::ostringstream oss;
  oss << std::put_time(&now_tm, "%Y-%m-%d-%H-%M-%S-") << rand();
  auto filename = "measurements/log-" + oss.str() + ".json";

  std::ofstream out(filename);

  out << "{\n";
  out << "  \"title\": \"" << this->header << "\",\n";
  out << "  \"algorithm\": \"" << algorithm_name(this->alogrithm) << "\",\n";
  out << "  \"min\": " << this->min_time.count() << ",\n";
  out << "  \"max\": " << this->max_time.count() << ",\n";
  out << "  \"avg\": " << (this->total_time / this->data_count).count()
      << ",\n";
  out << "  \"min_point_count\": " << this->min_points << ",\n";
  out << "  \"max_point_count\": " << this->max_points << ",\n";
  out << "  \"min_dimension_count\": " << this->min_dim << ",\n";
  out << "  \"max_dimension_count\": " << this->max_dim << "\n";
  out << "}";
  out.close();
}

} // namespace Log
