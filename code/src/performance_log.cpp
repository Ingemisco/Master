#include "datastructures.h"
#include "log.h"
#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <random>
#include <stdexcept>

namespace Log {
std::string measurement_directory = "";

PerformanceLogger::PerformanceLogger(Algorithm algorithm, std::string header)
    : header(header), algorithm(algorithm) {}

void PerformanceLogger::add_data(DataStructures::Polyline &polyline,
                                 Simplification::Simplification &simplification,
                                 std::chrono::duration<double> elapsed_time,
                                 std::string name) {

  this->times.push_back(elapsed_time);
  this->simplification_sizes.push_back(simplification->size());
  this->dimension = std::max(this->dimension, polyline.dimension);
  this->point_count = std::max(this->point_count, polyline.point_count);
  this->case_names.push_back(name);
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

	default:
		__builtin_unreachable();
  }
}

void PerformanceLogger::emit() {
  if (this->times.size() == 0) {
    throw std::runtime_error("No data has been entered.");
  }
  if (!std::filesystem::exists("measurements")) {
    std::filesystem::create_directory("measurements");
  }
  if (!std::filesystem::is_directory("measurements")) {
    throw std::runtime_error("File called 'measurements' exists. Cannot create "
                             "a directory of same name.");
  }

  std::string dir = "measurements";
  if (measurement_directory != "") {
    dir +=
        "/" + std::filesystem::path(measurement_directory).filename().string();
    if (!std::filesystem::exists(dir)) {
      std::filesystem::create_directory(dir);
    }
    if (!std::filesystem::is_directory(dir)) {
      throw std::runtime_error("File called '" + dir +
                               "' exists. Cannot create "
                               "a directory of same name.");
    }
  }

  // Get the current time
  auto now = std::chrono::system_clock::now();
  auto now_time_t = std::chrono::system_clock::to_time_t(now);
  auto now_tm = *std::localtime(&now_time_t);

  // Format the timestamp
  std::ostringstream oss;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, 10000000);
  oss << std::put_time(&now_tm, "%Y-%m-%d-%H-%M-%S-") << dis(gen);
  auto filename = dir + "/log-" + oss.str() + ".json";

  std::ofstream out(filename);
  auto minmax = std::ranges::minmax_element(this->times);
  auto avg =
      std::accumulate(this->times.begin(), this->times.end(), Time(0.0)) /
      this->times.size();

  out << "{\n";
  out << "  \"title\": \"" << this->header << "\",\n";
  out << "  \"algorithm\": \"" << algorithm_name(this->algorithm) << "\",\n";
  out << "  \"min\": " << minmax.min->count() << ",\n";
  out << "  \"max\": " << minmax.max->count() << ",\n";
  out << "  \"avg\": " << avg.count() << ",\n";
  out << "  \"dimension\": " << this->dimension << ",\n";
  out << "  \"point_count\": " << this->point_count << ",\n";

  out << "  \"data\": [\n";
  for (unsigned int i = 0; i < this->case_names.size(); i++) {
    out << "    { \"file\": \"" << this->case_names[i]
        << "\", \"time\": " << this->times[i].count()
        << ", \"simplification_size\": " << this->simplification_sizes[i]
        << " }";
    if (i < this->case_names.size() - 1) {
      out << ",";
    }

    out << "\n";
  }
  out << "  ]\n";

  out << "}";
  out.close();
}

} // namespace Log
