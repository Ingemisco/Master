#ifndef INCLUDE_INCLUDE_DATASTRUCTURES_H_
#define INCLUDE_INCLUDE_DATASTRUCTURES_H_

#include <cstddef>
#include <filesystem>

namespace DataStructures {

// does not own the data. Only lets you view a point in a polyline
struct Point final {
  float *const coordinates;
  size_t const dimension;

  Point(float *, size_t);

  float &operator[](size_t);
  float operator[](size_t) const;
};

std::ostream &operator<<(std::ostream &, Point const &);

// only used as logical grouping
struct LineSegment final {
  Point const &start;
  Point const &end;

  LineSegment(Point const &, Point const &);
};

struct Polyline final {
  // matrix containing all points
  float *const data;
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

// only used as logical grouping for Alt and Godau Algorithm and Polyline
// Simplification algorithms
struct PolylineRange final {
  Polyline &polyline;
  size_t const start_point_index;
  size_t const end_point_index;
  float const start_point_offset;

  PolylineRange(Polyline &, size_t, size_t, float);
};

void assert_compatible_points(Point const &, Point const &);
void assert_compatible_polylines(Polyline const &, Polyline const &);
} // namespace DataStructures

#endif // INCLUDE_INCLUDE_DATASTRUCTURES_H_
