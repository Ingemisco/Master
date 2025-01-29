#include "datastructures.h"
#include <cstddef>
#include <stdexcept>
#include <string>

#include "config.h"

using namespace DataStructures;

Point::Point(size_t dimension) { this->coordinates.resize(dimension); }

float &Point::operator[](size_t index) { return this->coordinates[index]; }

float Point::operator[](size_t index) const { return this->coordinates[index]; }

size_t Point::dimension() const { return this->coordinates.size(); }

Polyline::Polyline(size_t min_size) { this->points.reserve(min_size); }

void Polyline::add_point(Point &point) { this->points.push_back(point); }

Point &Polyline::operator[](size_t index) { return this->points[index]; }

Point interpolate(Point const &point1, Point const &point2, float ratio) {
  size_t const dimension = point1.dimension();
#if DEBUG
  if (ratio > 1 || ratio < 0) {
    throw std::runtime_error("Value of ratio is not between 0 and 1.");
  } else if (dimension != point2.dimension()) {
    std::string message =
        "Dimensions of points do not match. First point has dimension " +
        std::to_string(dimension) + " but second point has dimension " +
        std::to_string(point2.dimension()) + ".";
    throw std::runtime_error(message);
  }
#endif

  Point p(dimension);

  for (unsigned int i = 0; i < dimension; i++) {
    p[i] = (1 - ratio) * point1[i] + ratio * point2[i];
  }

  return p;
}
