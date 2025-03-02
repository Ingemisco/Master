#include "datastructures.h"
#include "distance.h"
#include "simplification.h"
#include <boost/program_options.hpp>
#include <boost/program_options/detail/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
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
          "a distance of at most epsilon.");

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

  auto polyline = DataStructures::Polyline::from_file(
      std::filesystem::path(poly_line_file_name));
  std::cout << *polyline << std::endl;

  if (map.count("se")) {
    auto simplification_vertices =
        Simplification::simplification_naive_euclidean(*polyline, epsilon);
    std::cout << "Simplification uses: ";
    for (size_t i : *simplification_vertices) {
      std::cout << i << ", ";
    }
    std::cout << std::endl;
  }
}

int main(int argc, char *argv[]) {
  handle_command_line_arguments(argc, argv);
  return 0;
}
