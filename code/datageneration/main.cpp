#include "generators.h"
#include "log.h"
#include "simplification.h"
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
#include <fstream>

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

static void test() {
	std::random_device rd;
	// uint32_t initial_seed = 4141310102; // = rd();
	// uint32_t initial_seed = 373139026;
	// uint32_t initial_seed = 2017988661;
	//uint32_t initial_seed = 2425774786;
	uint32_t initial_seed = rd();
	std::cout << "Initial seed: " << initial_seed << std::endl;
	// 
	// Save seed to file for later reproduction
	std::ofstream seed_file("data/custom/initial_seed.txt");
	seed_file << initial_seed;
	seed_file.close();

	std::mt19937 rng(initial_seed);
	std::uniform_real_distribution<float> epsilon_dist(0.01f, 100.0f);
	std::uniform_real_distribution<float> length_dist(0.1f, 100.0f);
	std::uniform_real_distribution<float> angle_dist(0.0f, 3.14159265f);

	AlgorithmConfiguration config = {
		.output_visualization = false,
		.logger = std::nullopt,
	};

	int iterations = 100000;
	int max_diff = -1;
	float best_epsilon;
	std::unique_ptr<DataStructures::Polyline> best_polyline = nullptr;

	int total_diff_sum = 0;

	for (int i = 0; i < iterations; ++i) {
		if (i % 1000 == 0) {
			std::cout << "Iteration " << i << std::endl;
		}
		float min_len = length_dist(rng);
		float max_len = length_dist(rng);
		if (min_len > max_len) std::swap(min_len, max_len);
		float angle = angle_dist(rng);

		auto poly_ = DataGeneration::make_polyline(100, 2, min_len, max_len, angle, rng);
		auto &poly = *poly_;

		float epsilon = epsilon_dist(rng);

		//if (i < 7433) {
		//	continue;
		//}
		// std::cout << "epsilon: " << epsilon << std::endl;
		// std::cout << "i: " << i << std::endl;

		//DataGeneration::write_to_file(poly, "data/custom/last_poly");
		auto result1 = Simplification::simplification_advanced_euclidean_explicit(poly, epsilon, config);
		auto result2 = Simplification::simplification_global_imai_iri_euclidean(poly, epsilon, config);
		//auto result2 = Simplification::simplification_imai_iri_euclidean(poly, epsilon, config);

		int size1 = static_cast<int>(result1->size());
		int size2 = static_cast<int>(result2->size());

		int diff = size2 - size1;
		if (diff < 0) {
			std::cout << "Something went wrong on iteration "<< i << " with a negative difference of "<< diff << std::endl;
			exit(0);
		}
		if (diff > max_diff) {
			std::cout << "Found with diff " << diff << " for epsilon " << epsilon << std::endl;
			max_diff = diff;
			best_polyline = std::move(poly_);
			best_epsilon = epsilon;
		}
		total_diff_sum += diff;
	}

	std::cout << "total diff to count: " << total_diff_sum << std::endl;

	if (total_diff_sum > 0) {
		std::filesystem::create_directories("data/custom");
		DataGeneration::write_to_file(*best_polyline, "data/custom/high_diff_poly" + std::to_string(best_epsilon));
		std::cout << "Best polyline written with difference " << max_diff << "\n";
	} else {
		std::cout << "No diff polyline found.\n";
	}
}

int main(int argc, char *argv[]) {
	test();
  //handle_command_line_arguments(argc, argv);
  return 0;
}
