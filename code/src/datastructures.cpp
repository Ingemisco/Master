#include "datastructures.h"
#include <cstddef>
#include <fstream>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "config.h"

using DataStructures::Polyline;

namespace DataStructures {

Polyline::Polyline(PointCount point_count, Dimension dimension)
    : data(new float[point_count * dimension]), point_count(point_count),
      dimension(dimension) {
#if DEBUG
  if (dimension == 0) {
    throw std::runtime_error(
        "Dimension of 0 is not allowed when constructing a Polyline!");
  } else if (point_count == 0) {
    throw std::runtime_error("A polyline must consist of at least 1 point!");
  }
#endif
}

Polyline::~Polyline() { delete[] this->data; }

float &Polyline::operator[](PointIndex point, Coordinate coordinate) {
#if DEBUG
  if (point >= this->point_count || coordinate >= this->dimension) {
    throw std::runtime_error("Accessed coordinate is not in bounds. Accessed " +
                             std::to_string(point) + ", " +
                             std::to_string(coordinate) +
                             " but amount of points and dimension is " +
                             std::to_string(this->point_count) + "," +
                             std::to_string(this->dimension) + ".");
  }
#endif
  return this->data[point * dimension + coordinate];
}

float Polyline::operator[](PointIndex point, Coordinate coordinate) const {
#if DEBUG
  if (point >= this->point_count || coordinate >= this->dimension) {
    throw std::runtime_error("Accessed coordinate is not in bounds. Accessed " +
                             std::to_string(point) + ", " +
                             std::to_string(coordinate) +
                             " but amount of points and dimension is " +
                             std::to_string(this->point_count) + "," +
                             std::to_string(this->dimension) + ".");
  }
#endif
  return this->data[point * dimension + coordinate];
}

std::unique_ptr<Polyline> Polyline::from_file(std::filesystem::path path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("Failed to open file: " + path.string());
  }

  // allow comments starting with '#' and empty lines at the beginning of the
  // file
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line[0] != '#') {
      break;
    }
  }

  PointCount point_count;
  Dimension dimension;
  std::istringstream first_line_stream(line);

  if (!(first_line_stream >> point_count >> dimension)) {
    throw std::runtime_error("Fist line: " + path.string());
  }

  auto polyline = std::make_unique<Polyline>(point_count, dimension);
  for (PointIndex point = 0; point < point_count; point++) {
    for (Coordinate coordinate = 0; coordinate < dimension; coordinate++) {
      if (!(input >> (*polyline)[point, coordinate])) {
        throw std::runtime_error("Invalid file format: Not enough floats.");
      }
    }
  }

  return polyline;
}

std::ostream &operator<<(std::ostream &os, Polyline &polyline) {
  os << "[";
  for (PointIndex point = 0; point < polyline.point_count; point++) {
    os << ",\n  (";
		for (Coordinate i = 0; i < polyline.dimension; i++) {
			os << polyline[point, i] << ", ";
		}
		os << ")";
  }
  os << "]";

  return os;
}

} // namespace DataStructures
