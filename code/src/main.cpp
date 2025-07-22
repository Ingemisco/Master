#include "datastructures.h"
#include "distance.h"
#include "log.h"
#include "queries.h"
#include "simplification.h"
#include <boost/program_options.hpp>
#include <boost/program_options/detail/parsers.hpp>
#include <boost/program_options/errors.hpp>
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

template <Simplification::Simplification _simplification_algorithm(DataStructures::Polyline &, float),
          Log::Algorithm _algorithm>
static inline void _flag_action_simplify(po::variables_map &map, char const *flag) {

	const auto &args = map[flag].as<std::vector<std::string>>();
	if (args.size() != 2) {
		throw po::error("Flag requires exactly two inputs: file and epsilon!");
	}
	std::string const &poly_line_file_name = args[0];
	float const epsilon = std::stof(args[1]);
  Log::PerformanceLogger log(_algorithm, poly_line_file_name);
  if (std::filesystem::is_directory(poly_line_file_name)) {
    Log::measurement_directory = poly_line_file_name;
    for (auto const &entry : std::filesystem::directory_iterator(poly_line_file_name)) {
      auto polyline = DataStructures::Polyline::from_file(std::filesystem::path(entry));

      auto start = std::chrono::high_resolution_clock::now();
      auto simplification_vertices = _simplification_algorithm(*polyline, epsilon);
      auto end = std::chrono::high_resolution_clock::now();
      log.add_data(*polyline, simplification_vertices, end - start, entry.path().filename().string());
    }

  } else {
    auto polypath = std::filesystem::path(poly_line_file_name);
    auto polyline = DataStructures::Polyline::from_file(polypath);

    auto start = std::chrono::high_resolution_clock::now();
    auto simplification_vertices = _simplification_algorithm(*polyline, epsilon);
    auto end = std::chrono::high_resolution_clock::now();
    log.add_data(*polyline, simplification_vertices, end - start, polypath.filename().string());
		std::cout << "size: " << simplification_vertices->size() << std::endl;
  }
  log.emit();
}

static inline void handle_command_line_arguments(int argc, char *argv[]) {
  po::positional_options_description positional_options;
  positional_options.add("file", 1);

  po::options_description description("Allowed options");
  std::string poly_line_file_name;

  auto options = description.add_options();

  options("help,h", "Show help message");

  options("se",
					po::value<std::vector<std::string>>()->multitoken()->value_name("filename epsilon"),
          "Uses the Van Kreveld et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Euclidean distance.");

  options("sei",
					po::value<std::vector<std::string>>()->multitoken()->value_name("filename epsilon"),
          "Uses the Van Kreveld et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Euclidean distance. Does not "
          "explicitly compute zeros but only implicitly compare them.");

  options("ses",
					po::value<std::vector<std::string>>()->multitoken()->value_name("filename epsilon"),
          "Uses the Van Kreveld et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Euclidean distance. "
					"Uses Semiepxlicit computations. The second input is not epsilon but"
					"epsilon squared!");

  options("sm",
					po::value<std::vector<std::string>>()->multitoken()->value_name("filename epsilon"),
          "Uses the Van Kreveld et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Manhattan distance.");

  options("sc",
					po::value<std::vector<std::string>>()->multitoken()->value_name("filename epsilon"),
          "Uses the Van Kreveld et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Chebyshev distance.");

  options("ae",
					po::value<std::vector<std::string>>()->multitoken()->value_name("filename epsilon"),
          "Uses the Bringmann et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Euclidean distance.");

  options("aes",
					po::value<std::vector<std::string>>()->multitoken()->value_name("filename epsilon"),
          "Uses the Bringmann et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Euclidean distance. " 
					"Uses semiexplicit computations.");

  options("am",
					po::value<std::vector<std::string>>()->multitoken()->value_name("filename epsilon"),
          "Uses the Bringmann et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Manhattan distance.");

  options("ac",
					po::value<std::vector<std::string>>()->multitoken()->value_name("filename epsilon"),
          "Uses the Bringmann et al. algorithm to simplify the polyline with "
          "a distance of at most epsilon using Chebyshev distance.");

  options("bes",
					po::value<std::string>(&poly_line_file_name)->value_name("filename"),
          "Builds the datastructure to allow fast simplification queries"
					"for Euclidean distance. Simple version (n^4 space consumption).");

  options("qe", 
					po::value<std::vector<std::string>>()->multitoken()->value_name("filename epsilon"),
					"Queries for a simplification given a file to the saved datastructure and an epsilon.");

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

	if(map.count("bes")) {
		auto path = std::filesystem::path(poly_line_file_name);
		auto polyline = DataStructures::Polyline::from_file(path);
		auto ds = build_querier_simple(*polyline);
		auto file_name = path.filename().string();
		std::filesystem::path save_path = std::filesystem::path("datastructures") / file_name;
		ds->save_datastructure_to_file(save_path);
	} else if(map.count("qe")) {
		const auto &args = map["qe"].as<std::vector<std::string>>();
		if (args.size() != 2) {
			throw po::error("Flag requires exactly two inputs: file and epsilon!");
		}
		// std::string const &file_name = args[0];
		// float const epsilon = std::stof(args[1]);
		// auto querier = SimplificationQuerier::from_file(file_name);
		// querier->print();


	
	} else if (map.count("se")) {
    _flag_action_simplify<Simplification::simplification_naive_euclidean,
                          Log::Algorithm::SIMPLIFICATION_SIMPLE_EUCLIDEAN>( map, "se");
  } else if (map.count("sei")) {
    _flag_action_simplify<
        Simplification::simplification_naive_euclidean_implicit,
        Log::Algorithm::SIMPLIFICATION_SIMPLE_IMPLICIT_EUCLIDEAN>( map, "sei");
  } else if (map.count("ses")) {
    _flag_action_simplify<
        Simplification::simplification_naive_euclidean_semiexplicit,
        Log::Algorithm::SIMPLIFICATION_SIMPLE_SEMIEXPLICIT_EUCLIDEAN>( map, "ses");
  } else if (map.count("sm")) {
    _flag_action_simplify<Simplification::simplification_naive_manhattan,
                          Log::Algorithm::SIMPLIFICATION_SIMPLE_MANHATTAN>( map, "sm");
  } else if (map.count("sc")) {
    _flag_action_simplify<Simplification::simplification_naive_chebyshev,
                          Log::Algorithm::SIMPLIFICATION_SIMPLE_CHEBYSHEV>( map, "sc");
	}	else if (map.count("ae")) {
    _flag_action_simplify<Simplification::simplification_advanced_euclidean_explicit,
                          Log::Algorithm::SIMPLIFICATION_ADVANCED_EUCLIDEAN>( map, "ae");
	}	else if (map.count("am")) {
    _flag_action_simplify<Simplification::simplification_advanced_manhattan_explicit,
                          Log::Algorithm::SIMPLIFICATION_ADVANCED_MANHATTAN>( map, "am");
	}	else if (map.count("ac")) {
    _flag_action_simplify<Simplification::simplification_advanced_chebyshev_explicit,
                          Log::Algorithm::SIMPLIFICATION_ADVANCED_CHEBYSHEV>( map, "ac");
	}	else if (map.count("aes")) {
    _flag_action_simplify<Simplification::simplification_advanced_euclidean_semiexplicit,
                          Log::Algorithm::SIMPLIFICATION_ADVANCED_SEMIEXPLICIT_EUCLIDEAN>( map, "aes");
  }
}

int main(int argc, char *argv[]) {
  handle_command_line_arguments(argc, argv);
  return 0;
}
