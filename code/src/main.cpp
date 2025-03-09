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
#include <utility>

namespace po = boost::program_options;

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
    Log::PerformanceLogger log(Log::Algorithm::SIMPLIFICATION_SIMPLE_EUCLIDEAN,
                               poly_line_file_name);
    if (std::filesystem::is_directory(poly_line_file_name)) {
      for (auto const &entry :
           std::filesystem::directory_iterator(poly_line_file_name)) {
        std::cout << "Test File " << entry << " in directory." << std::endl;
        auto polyline =
            DataStructures::Polyline::from_file(std::filesystem::path(entry));

        auto start = std::chrono::high_resolution_clock::now();
        auto simplification_vertices =
            Simplification::simplification_naive_euclidean(*polyline, epsilon);
        auto end = std::chrono::high_resolution_clock::now();
        log.add_data(*polyline, end - start);
      }

    } else {
      std::cout << "Test File " << poly_line_file_name << std::endl;
      auto polyline = DataStructures::Polyline::from_file(
          std::filesystem::path(poly_line_file_name));

      auto start = std::chrono::high_resolution_clock::now();
      auto simplification_vertices =
          Simplification::simplification_naive_euclidean(*polyline, epsilon);
      auto end = std::chrono::high_resolution_clock::now();
      log.add_data(*polyline, end - start);
    }
    log.emit();
  }
}

int main(int argc, char *argv[]) {
  handle_command_line_arguments(argc, argv);
  return 0;
}
