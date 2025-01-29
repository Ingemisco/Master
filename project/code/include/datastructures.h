#ifndef INCLUDE_INCLUDE_DATASTRUCTURES_H_
#define INCLUDE_INCLUDE_DATASTRUCTURES_H_

#include <cstddef>
#include <vector>

template <typename T> class Point {
private:
  std::vector<T> coordinates_;

public:
  Point(size_t dimension);
  T &operator[](size_t);
};

template <typename T> class Polyline {
private:
public:
  Polyline();
};

#endif // INCLUDE_INCLUDE_DATASTRUCTURES_H_
