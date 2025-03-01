#include "datastructures.h"
#include "distance.h"
#include "simplification.h"
#include <boost/program_options.hpp>
#include <boost/program_options/detail/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <filesystem>
#include <iostream>

namespace po = boost::program_options;

static inline void handle_command_line_arguments(int argc, char *argv[]) {
  po::options_description description("Allowed options");
  std::string poly_line_file_name;

  auto options = description.add_options();

  options("help,h", "Show help message");
  options("read,r", po::value<std::string>(&poly_line_file_name),
          "Read a polyline from a file");

  po::variables_map map;
  po::store(po::parse_command_line(argc, argv, description), map);

  po::notify(map);

  if (map.count("help")) {
    std::cout << description;
    exit(0);
  } else if (map.count("read")) {
    auto polyline = DataStructures::Polyline::from_file(
        std::filesystem::path(poly_line_file_name));
    std::cout << *polyline << std::endl;

    // auto res = DataStructures::solve_manhattan(polyline->get_point(0),
    //                                            polyline->get_point(1),
    //                                            polyline->get_point(2), 5);
    // std::cout << res.first << ", " << res.last << std::endl;

    std::cout << std::endl << "Testing algo " << std::endl;
    auto v = Simplification::simplification_naive_euclidean(*polyline, 0.5);

    std::cout << "Simplification uses: ";
    for (size_t i : *v) {
      std::cout << i << ", ";
    }
    std::cout << std::endl;

    // float res = DataStructures::alt_godau_euclidean(
    //     DataStructures::PolylineRange(*polyline, 2, 4, 0),
    //     DataStructures::LineSegment(polyline->get_point(2),
    //                                 polyline->get_point(4)),
    //     1);
    //
    // std::cout << "result " << res << std::endl;

    // std::cout << DataStructures::alt_godau_manhattan(
    //                  DataStructures::PolylineRange(*polyline, 0, 3, 0.1),
    //                  DataStructures::LineSegment(polyline->get_point(0),
    //                                              polyline->get_point(3)),
    //                  5)
    //           << std::endl;
  }
}

int main(int argc, char *argv[]) {
  handle_command_line_arguments(argc, argv);
  return 0;
}
