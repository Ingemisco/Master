#include "datastructures.h"
#include "distance.h"
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

  description.add_options()           //
      ("help,h", "Show help message") //
      ("read,r", po::value<std::string>(&poly_line_file_name),
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
  }
}

int main(int argc, char *argv[]) {
  handle_command_line_arguments(argc, argv);
  return 0;
}
