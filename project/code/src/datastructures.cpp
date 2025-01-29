#include "datastructures.h"
#include <cstddef>

Point::Point(size_t dimension) { this->coordinates_.resize(dimension); }

float &Point::operator[](size_t index) { return this->coordinates_[index]; }

Polyline::Polyline(size_t min_size) { this->points_.reserve(min_size); }

void Polyline::add_point(Point &point) { this->points_.push_back(point); }

Point &Polyline::operator[](size_t index) { return this->points_[index]; }
