#include "generators.h"
#include <boost/program_options.hpp>
#include <boost/program_options/detail/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <filesystem>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

namespace po = boost::program_options;

static inline void handle_command_line_arguments(int argc, char *argv[]) {
  po::positional_options_description positional_options;
  positional_options.add("point_count", 1);
  positional_options.add("dimension", 1);

  po::options_description description("Allowed options");

  std::string poly_line_file_name;
  int count;
  float min_length = 0;
  float max_length = 10;
  int max_angle;
	std::string suite_name = "";

  auto options = description.add_options();

  options("help,h", "Show help message");

  options("point_count", po::value<int>()->required(),
          "Point count of polyline");
  options("dimension", po::value<int>()->required(), "Dimension of polyline");

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

  options(
      "angle,a",
      po::value<int>(&max_angle)->default_value(180)->value_name("max_angle"),
      "Bounds the maximum angle that two consecutive line segments in the "
      "polyline are allowed to have. Specified as an integer in degree from 1 "
      "to 180.");


  options(
      "suite,s",
      po::value<std::string>(&suite_name)->default_value("")->value_name("suite_name"),
			"creates a test suite with the given directory name. Uses the count flag for the amount of "
			"data sets per case and the point_count flag for the max amount of points (goes in steps of 10)"
			". Most flags are ignored."
      );

  options("integral,i", "Rounds every coordinate to the nearest integer. May "
                        "slightly violate angles and length constraints.");

  po::variables_map map;
  po::store(po::command_line_parser(argc, argv)
                .options(description)
                .positional(positional_options)
                .run(),
            map);

  if (map.count("help")) {
    std::cout << "Usage: ./datagen pointcount dimension [flags]" << std::endl;
    std::cout << description << std::endl;
    exit(0);
  }

  po::notify(map);

  if (map.count("length")) {
    auto const &values = map["length"].as<std::vector<float>>();
    switch (values.size()) {
    case 1:
      max_length = values[0];
      break;
    case 2:
      min_length = values[0];
      max_length = values[1];
      break;
    default:
      std::cerr
          << "The length flag needs exactly one or two numbers as arguments"
          << std::endl;
      exit(1);
    }
  }

  if (min_length < 0 || max_length <= 0) {
    std::cerr << "The lengths specified are not allowed to be negative and the "
                 "maximum length must not be zero.";
    exit(1);
  } else if (count <= 0) {
    std::cerr << "Amount of files to be generated must be at least one.";
    exit(1);
  } else if (max_angle <= 0 || max_angle > 180) {
    std::cerr << "Angle must be in between (0, 180].";
    exit(1);
  }

  int point_count = map["point_count"].as<int>();
  int dimension = map["dimension"].as<int>();

  std::random_device rd;
  std::mt19937 gen(rd());

	if (suite_name != "") {
		DataGeneration::make_test_suite(suite_name, 10, 10, point_count, count, gen);
		return;
	}

  auto polyline = DataGeneration::make_polyline(
      point_count, dimension, min_length, max_length,

      max_angle / 180.0 * std::numbers::pi, gen);

  if (map.count("integral") > 0) {
    DataGeneration::make_integral(*polyline.get());
  }

  if (count == 1 && !std::filesystem::exists(poly_line_file_name)) {
    DataGeneration::write_to_file(*polyline.get(),
                                  std::filesystem::path(poly_line_file_name));
  } else {
    unsigned int j = 0;
    for (int i = 0; i < count; i++) {
      while (std::filesystem::exists(poly_line_file_name + std::to_string(j))) {
        j++;
      }

      auto polyline = DataGeneration::make_polyline(
          point_count, dimension, min_length, max_length,

          max_angle / 180.0 * std::numbers::pi, gen);
      DataGeneration::write_to_file(
          *polyline.get(),
          std::filesystem::path(poly_line_file_name + std::to_string(j)));
      if (i < count - 1) {
        polyline = DataGeneration::make_polyline(
            point_count, dimension, min_length, max_length,

            max_angle / 180.0 * std::numbers::pi, gen);

        if (map.count("integral") > 0) {
          DataGeneration::make_integral(*polyline.get());
        }
      }
    }
  }
}

int main(int argc, char *argv[]) {
  handle_command_line_arguments(argc, argv);
  return 0;
}
