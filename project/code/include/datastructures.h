#ifndef INCLUDE_INCLUDE_DATASTRUCTURES_H_
#define INCLUDE_INCLUDE_DATASTRUCTURES_H_

#include <cstddef>
#include <vector>

namespace DataStructures {

class Point {
private:
  std::vector<float> coordinates;

public:
  Point(size_t dimension);

  inline float &operator[](size_t);
  inline float operator[](size_t) const;
  inline size_t dimension() const;
};

class Polyline {
private:
  std::vector<Point> points;

public:
  Polyline(size_t = 0);

  inline void add_point(Point &);
  inline Point &operator[](size_t);
};

Point interpolate(Point const &, Point const &, float);

} // namespace DataStructures

#endif // INCLUDE_INCLUDE_DATASTRUCTURES_H_
