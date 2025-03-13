#include "datastructures.h"
#include "distance.h"
#include "log.h"
#include "simplification.h"
#include <boost/program_options.hpp>
#include <boost/program_options/detail/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

namespace po = boost::program_options;

template <Simplification::Simplification _simplification_algorithm(
              DataStructures::Polyline &, float),
          Log::Algorithm _algorithm>
static inline void _flag_action_simplify(std::string &poly_line_file_name,
                                         float epsilon) {
  Log::PerformanceLogger log(_algorithm, poly_line_file_name);
  if (std::filesystem::is_directory(poly_line_file_name)) {
    for (auto const &entry :
         std::filesystem::directory_iterator(poly_line_file_name)) {
      auto polyline =
          DataStructures::Polyline::from_file(std::filesystem::path(entry));

      auto start = std::chrono::high_resolution_clock::now();
      auto simplification_vertices =
          _simplification_algorithm(*polyline, epsilon);
      auto end = std::chrono::high_resolution_clock::now();
      log.add_data(*polyline, end - start);
    }

  } else {
    auto polyline = DataStructures::Polyline::from_file(
        std::filesystem::path(poly_line_file_name));

    auto start = std::chrono::high_resolution_clock::now();
    auto simplification_vertices =
        _simplification_algorithm(*polyline, epsilon);
    auto end = std::chrono::high_resolution_clock::now();
    log.add_data(*polyline, end - start);
  }
  log.emit();
}

static inline void handle_command_line_arguments(int argc, char *argv[]) {
  po::positional_options_description positional_options;
  positional_options.add("file", 1);
  positional_options.add("epsilon", 1);

  po::options_description description("Allowed options");
  std::string poly_line_file_name;
  float epsilon;

  auto options = description.add_options();

  options("help,h", "Show help message");

  options("file",
          po::value<std::string>(&poly_line_file_name)
              ->value_name("file")
              ->required(),
          "File from which to read the polyline");

  options("epsilon",
          po::value<float>(&epsilon)->value_name("epsilon")->required(),
          "Maximal Distance the simplification is allowed to have to the "
          "original polyline.");

  options("se",
          "Uses the Van Kreveld et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Euclidean distance.");

  options("sei",
          "Uses the Van Kreveld et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Euclidean distance. Does not "
          "explicitly compute zeros but only implicitly compare them.");

  options("sm",
          "Uses the Van Kreveld et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Manhattan distance.");

  options("sc",
          "Uses the Van Kreveld et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Chebyshev distance.");

  po::variables_map map;
  po::store(po::command_line_parser(argc, argv)
                .options(description)
                .positional(positional_options)
                .run(),
            map);

  if (map.count("help")) {
    std::cout << description;
    exit(0);
  }

  po::notify(map);

  if (map.count("se")) {
    _flag_action_simplify<Simplification::simplification_naive_euclidean,
                          Log::Algorithm::SIMPLIFICATION_SIMPLE_EUCLIDEAN>(
        poly_line_file_name, epsilon);
  } else if (map.count("sei")) {
    _flag_action_simplify<
        Simplification::simplification_naive_euclidean_implicit,
        Log::Algorithm::SIMPLIFICATION_SIMPLE_IMPLICIT_EUCLIDEAN>(
        poly_line_file_name, epsilon);
  } else if (map.count("sm")) {
    _flag_action_simplify<Simplification::simplification_naive_manhattan,
                          Log::Algorithm::SIMPLIFICATION_SIMPLE_MANHATTAN>(
        poly_line_file_name, epsilon);
  } else if (map.count("sc")) {
    _flag_action_simplify<Simplification::simplification_naive_chebyshev,
                          Log::Algorithm::SIMPLIFICATION_SIMPLE_CHEBYSHEV>(
        poly_line_file_name, epsilon);
  }

  // TESTING
  // auto polyline = DataStructures::Polyline::from_file(
  //     std::filesystem::path(poly_line_file_name));
  // auto &p = *polyline;
  // auto res = DataStructures::alt_godau_euclidean_implicit(p, 2, 2, 0, 1, 4,
  //                                                         epsilon * epsilon);
  // std::cout << p << std::endl;
  // std::cout << res << std::endl;
}

int main(int argc, char *argv[]) {
  handle_command_line_arguments(argc, argv);
  return 0;
}
