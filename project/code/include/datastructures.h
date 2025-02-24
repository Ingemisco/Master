#ifndef INCLUDE_INCLUDE_DATASTRUCTURES_H_
#define INCLUDE_INCLUDE_DATASTRUCTURES_H_

#include <cstddef>
#include <filesystem>

namespace DataStructures {

// does not own the data. Only lets you view a point in a polyline
class Point {
private:
  float *const coordinates;

public:
  size_t const dimension;

  Point(float *, size_t);

  float &operator[](size_t);
  float operator[](size_t) const;
};

std::ostream &operator<<(std::ostream &, Point const &);

class Polyline {
private:
  // matrix containing all points
  float *const data;

public:
  size_t const point_count;
  size_t const dimension;

  Polyline(size_t, size_t);
  ~Polyline();

  Polyline(Polyline const &) = delete;
  Polyline(Polyline &&) = delete;
  Polyline operator=(Polyline const &) = delete;
  Polyline operator=(Polyline &&) = delete;

  float &operator[](size_t, size_t);
  float operator[](size_t, size_t) const;

  Point get_point(size_t);
  Point const get_point(size_t) const;

  static std::unique_ptr<Polyline> from_file(std::filesystem::path);

  friend std::ostream &operator<<(std::ostream &, Polyline &);
};

void assert_compatible_points(Point const &, Point const &);
void assert_compatible_polylines(Polyline const &, Polyline const &);

} // namespace DataStructures

#endif // INCLUDE_INCLUDE_DATASTRUCTURES_H_
