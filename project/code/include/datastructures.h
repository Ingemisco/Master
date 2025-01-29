#ifndef INCLUDE_INCLUDE_DATASTRUCTURES_H_
#define INCLUDE_INCLUDE_DATASTRUCTURES_H_

#include <cstddef>
#include <vector>

class Point {
private:
  std::vector<float> coordinates_;

public:
  Point(size_t dimension);

  float &operator[](size_t);
};

class Polyline {
private:
  std::vector<Point> points_;

public:
  Polyline(size_t = 0);

  void add_point(Point &);
  Point &operator[](size_t);
};

#endif // INCLUDE_INCLUDE_DATASTRUCTURES_H_
