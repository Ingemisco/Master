#include "global.h"
#include "log.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <random>
#include <stdexcept>

namespace Log {
std::string measurement_directory = "";

void PerformanceLogger::add_data(size_t simplification_size, std::chrono::duration<double> elapsed_time, std::string name) {
	auto &data_set = this->data_sets.back();
  data_set.times.push_back(elapsed_time);
  data_set.simplification_sizes.push_back(simplification_size);
  data_set.case_names.push_back(name);
}

std::string algorithm_name(Algorithm algorithm) {
  switch (algorithm) {
  case Algorithm::SIMPLIFICATION_SIMPLE_MANHATTAN:
    return "Manhattan Simplification (Simple)";
  case Algorithm::SIMPLIFICATION_SIMPLE_EUCLIDEAN:
    return "Euclidean Simplification (Simple)";
  case Algorithm::SIMPLIFICATION_SIMPLE_IMPLICIT_EUCLIDEAN:
    return "Euclidean Simplification (Simple, Implicit)";
  case Algorithm::SIMPLIFICATION_SIMPLE_SEMIEXPLICIT_EUCLIDEAN:
    return "Euclidean Simplification (Simple, Semiexplicit)";
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
  case Algorithm::SIMPLIFICATION_ADVANCED_SEMIEXPLICIT_EUCLIDEAN:
    return "Euclidean Simplification (Advanced, Semiexplicit)";
  case Algorithm::SIMPLIFICATION_ADVANCED_CHEBYSHEV:
    return "Chebyshev Simplification (Advanced)";
  case Algorithm::SIMPLIFICATION_ADVANCED_IMPLICIT_MINKOWSKI:
    return "Minkowski Simplification (Advanced, Implicit)";
  case Algorithm::SIMPLIFICATION_IMAI_IRI_EUCLIDEAN:
    return "Local Euclidean Simplification (Imai Iri)";
  case Algorithm::BUILD_DS_EUCLIDEAN:
    return "Datastructure Euclidean";

	default:
		__builtin_unreachable();
  }
}

void PerformanceLogger::begin_data_set(Algorithm algorithm, Dimension dimension, PointCount point_count, std::string name) {
	this->data_sets.push_back({
		.header = name,
		.algorithm = algorithm,
		.dimension = dimension,
		.point_count = point_count,
		.times = {},
		.simplification_sizes = {},
		.case_names = {},
	});
}

void PerformanceLogger::emit() {
  if (this->data_sets.size() == 0) {
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


	bool start = true;

	out << "[\n";
	for (auto const &data_set : this->data_sets) {
		if (!start) {
			out << ", ";
		} else {
			start = false;
			out << "  ";
		}

		auto minmax = std::ranges::minmax_element(data_set.times);
		auto avg = std::accumulate(data_set.times.begin(), data_set.times.end(), Time(0.0)) / data_set.times.size();

		out << "{\n";
		out << "    \"title\": \"" << data_set.header << "\",\n";
		out << "    \"algorithm\": \"" << algorithm_name(data_set.algorithm) << "\",\n";
		out << "    \"min\": " << minmax.min->count() << ",\n";
		out << "    \"max\": " << minmax.max->count() << ",\n";
		out << "    \"avg\": " << avg.count() << ",\n";
		out << "    \"dimension\": " << data_set.dimension << ",\n";
		out << "    \"point_count\": " << data_set.point_count << ",\n";

		out << "    \"data\": [\n";
		for (unsigned int i = 0; i < data_set.case_names.size(); i++) {
			out << "      { \"file\": \"" << data_set.case_names[i]
					<< "\", \"time\": " << data_set.times[i].count()
					<< ", \"simplification_size\": " << data_set.simplification_sizes[i]
					<< " }";
			if (i < data_set.case_names.size() - 1) {
				out << ",";
			}

			out << "\n";
		}
		out << "    ]\n";

		out << "  }\n";

	}

	out << "]";
  out.close();

	this->data_sets.clear();
}

} // namespace Log
