#include <boost/program_options.hpp>
#include <boost/program_options/detail/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <iostream>

namespace po = boost::program_options;

static void handle_command_line_arguments(int argc, char *argv[]) {
  po::options_description description("Allowed options");
  description.add_options()           //
      ("help,h", "Show help message") //
      ;

  po::variables_map map;
  po::store(po::parse_command_line(argc, argv, description), map);

  if (map.count("help")) {
    std::cout << description;
    exit(0);
  }
}

int main(int argc, char *argv[]) {
  handle_command_line_arguments(argc, argv);
  return 0;
}
