#include "datastructures.h"
#include <cstddef>
#include <fstream>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "config.h"

namespace DataStructures {

Point::Point(float *coordinates, size_t dimension)
    : coordinates(coordinates), dimension(dimension) {}

float &Point::operator[](size_t index) {
#if DEBUG
  if (index > this->dimension) {
    throw std::runtime_error("Index accessed is not in bounds. Index is " +
                             std::to_string(index) + " but dimension is only " +
                             std::to_string(this->dimension) + ".");
  }
#endif
  return this->coordinates[index];
}

float Point::operator[](size_t index) const {
#if DEBUG
  if (index > this->dimension) {
    throw std::runtime_error("Index accessed is not in bounds. Index is " +
                             std::to_string(index) + " but dimension is only " +
                             std::to_string(this->dimension) + ".");
  }
#endif
  return this->coordinates[index];
}

std::ostream &operator<<(std::ostream &os, Point const &point) {
  os << "(" << point[0];
  for (unsigned int coordinate = 1; coordinate < point.dimension;
       coordinate++) {
    os << ", " << point[coordinate];
  }
  os << ")";

  return os;
}

Polyline::Polyline(size_t point_count, size_t dimension)
    : data(new float[point_count * dimension]), point_count(point_count),
      dimension(dimension) {}

Polyline::~Polyline() { delete[] this->data; }

float &Polyline::operator[](size_t point, size_t coordinate) {
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

float Polyline::operator[](size_t point, size_t coordinate) const {
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

Point Polyline::get_point(size_t index) {
#if DEBUG
  if (index >= this->point_count) {
    throw std::runtime_error(
        "Point accessed is not in bounds. Accessed " + std::to_string(index) +
        " but amount of points is " + std::to_string(this->point_count) + ".");
  }
#endif
  return Point(&this->data[index * this->dimension], this->dimension);
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

  size_t point_count;
  size_t dimension;
  std::istringstream first_line_stream(line);

  if (!(first_line_stream >> point_count >> dimension)) {
    throw std::runtime_error("Fist line: " + path.string());
  }

  auto polyline = std::make_unique<Polyline>(point_count, dimension);
  for (unsigned int point = 0; point < point_count; point++) {
    for (unsigned int coordinate = 0; coordinate < dimension; coordinate++) {
      if (!(input >> (*polyline)[point, coordinate])) {
        throw std::runtime_error("Invalid file format: Not enough floats.");
      }
    }
  }

  return polyline;
}

std::ostream &operator<<(std::ostream &os, Polyline &polyline) {
  os << "[" << polyline.get_point(0);
  for (unsigned int point = 1; point < polyline.point_count; point++) {
    os << ";\n" << polyline.get_point(point);
  }
  os << "]";

  return os;
}
} // namespace DataStructures
