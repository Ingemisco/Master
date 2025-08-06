#include "generators.h"
#include "datastructures.h"
#include "distance.h"
#include "global.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>

namespace DataGeneration {

// angle in radian in [0, pi]
std::unique_ptr<DataStructures::Polyline>
make_polyline(PointCount point_count, Dimension dimension, float min_length,
              float max_length, float angle, std::mt19937 &gen) {
  auto res_polyline =
      std::make_unique<DataStructures::Polyline>(point_count, dimension);
  auto &polyline = *res_polyline.get();

  std::uniform_real_distribution<float> full_angle_dist(0,
                                                        2 * std::numbers::pi);
  std::uniform_real_distribution<float> restr_angle_dist(-angle, angle);
  std::uniform_real_distribution<float> unif_dist(0, 1);

  float const min_d =
      DataStructures::integer_exponentiation(min_length, dimension);
  float const range_d =
      DataStructures::integer_exponentiation(max_length, dimension) - min_d;

  size_t const angle_count = dimension - 1;
  float *angles = new float[angle_count];
  for (unsigned int i = 0; i < angle_count; i++) {
    angles[i] = full_angle_dist(gen);
  }

  // first point fixed at zero, all others random
  for (unsigned int i = 0; i < dimension; i++) {
    polyline[0, i] = 0;
  }

  for (unsigned int i = 1; i < point_count; i++) {
    float value = std::pow(min_d + range_d * unif_dist(gen), 1.0 / dimension);
    for (unsigned int j = 0; j < dimension - 1; j++) {
      polyline[i, j] = polyline[i - 1, j] + value * std::cos(angles[j]);
      value *= std::sin(angles[j]);
    }
    polyline[i, dimension - 1] = polyline[i - 1, dimension - 1] + value;

    // restriction is only in first angle
    angles[0] += restr_angle_dist(gen);
    for (unsigned int j = 1; j < angle_count; j++) {
      angles[j] += full_angle_dist(gen);
    }
  }

  delete[] angles;
  return res_polyline;
}

void make_integral(DataStructures::Polyline &polyline) {
  for (unsigned int i = 0; i < polyline.point_count; i++) {
    for (unsigned int j = 0; j < polyline.dimension; j++) {
      polyline[i, j] = std::round(polyline[i, j]);
    }
  }
}

void write_to_file(DataStructures::Polyline const &polyline,
                   std::filesystem::path path) {
  std::ofstream output;
  output.open(path);
  if (!output.is_open()) {
    std::cerr << "Could not write to file " << path << std::endl;
    exit(1);
  }

  output << polyline.point_count << " " << polyline.dimension << "\n";
  for (unsigned int i = 0; i < polyline.point_count; i++) {
    for (unsigned int j = 0; j < polyline.dimension; j++) {
      output << polyline[i, j] << " ";
    }
    output << "\n";
  }
  std::cout << "Created file " << path << std::endl;
}

} // namespace DataGeneration
