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

LineSegment::LineSegment(Point const &start, Point const &end)
    : start(start), end(end) {
#if DEBUG
  assert_compatible_points(start, end);
#endif
}

Polyline::Polyline(size_t point_count, size_t dimension)
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

template <typename T>
static inline Point _get_point(T *data, size_t index, size_t dimension,
                               size_t point_count) {
#if DEBUG
  if (index >= point_count) {
    throw std::runtime_error(
        "Point accessed is not in bounds. Accessed " + std::to_string(index) +
        " but amount of points is " + std::to_string(point_count) + ".");
  }
#endif
  return Point(&data[index * dimension], dimension);
}

Point Polyline::get_point(size_t index) {
  return _get_point(this->data, index, this->dimension, this->point_count);
}

Point const Polyline::get_point(size_t index) const {
  return _get_point(this->data, index, this->dimension, this->point_count);
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
    os << ",\n  " << polyline.get_point(point);
  }
  os << "]";

  return os;
}

PolylineRange::PolylineRange(Polyline &polyline, size_t start_point,
                             size_t end_point, float start_offset)
    : polyline(polyline), start_point_index(start_point),
      end_point_index(end_point), start_point_offset(start_offset) {}

void assert_compatible_points(Point const &point1, Point const &point2) {
  if (point1.dimension != point2.dimension) {
    throw std::runtime_error(
        "Dimensions of points do not match. First point has dimension " +
        std::to_string(point1.dimension) + " but second point has dimension " +
        std::to_string(point2.dimension) + ".");
  }
}

void assert_compatible_polylines(Polyline const &polyline1,
                                 Polyline const &polyline2) {
  if (polyline1.dimension != polyline2.dimension) {
    throw std::runtime_error(
        "Dimensions of polylines do not match. First polyline has dimension " +
        std::to_string(polyline1.dimension) +
        " but second polyline has dimension " +
        std::to_string(polyline2.dimension) + ".");
  }
}
} // namespace DataStructures
