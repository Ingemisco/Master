#include <boost/program_options.hpp>
#include <boost/program_options/detail/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <iostream>

namespace po = boost::program_options;

static inline void handle_command_line_arguments(int argc, char *argv[]) {
  po::options_description description("Allowed options");
  std::string poly_line_file_name;
  int count;
  float min_length;
  float max_length;
  int max_angle;

  auto options = description.add_options();
  options("help,h", "Show help message");

  options("file,f",
          po::value<std::string>(&poly_line_file_name)
              ->default_value("polyline")
              ->value_name("filename"),
          "Write polyline into fie with this name. Default file name is "
          "'polyline'.");

  options("count,n",
          po::value<int>(&count)->default_value(1)->value_name("amount"),
          "Amount of files to be generated. Default is 1 "
          "and cannot be smaller "
          "than that. If more than one file is to be "
          "generated, all files will "
          "gain a number as suffix.");

  options(
      "length,l",
      po::value<std::vector<float>>()
          ->value_name("max_length / min_length max_length")
          ->multitoken(),
      "specifies the minimum and maximum length of line segments in the "
      "polyline. If only one value is given, this sets the max length and min "
      "will be 0. More than two values are not allowed.");

  options("fm", "Use the Manhattan distance for the lenth of the line segments "
                "instead of the Euclidean distance.");

  options("fc", "Use the Chebyshev distance for the lenth of the line segments "
                "instead of the Euclidean distance.");

  options("im",
          "Use the Manhattan distance for the lenth of the line segments "
          "instead of the Euclidean distance and forces all coordinates to "
          "be integers.");

  options("ic",
          "Use the Chebyshev for the lenth of the line segments "
          "instead of the Euclidean distance and forces all coordinates to "
          "be integers.");

  options(
      "a,angle",
      po::value<int>(&max_angle)->default_value(180)->value_name("max_angle"),
      "Bounds the maximum angle that two consecutive line segments in the "
      "polyline are allowed to have. Specified as an integer in degree from 1 "
      "to 180.");

  po::variables_map map;
  po::store(po::parse_command_line(argc, argv, description), map);

  po::notify(map);

  if (map.count("help")) {
    std::cout << description;
    exit(0);
  }
}

int main(int argc, char *argv[]) {
  handle_command_line_arguments(argc, argv);
  return 0;
}
