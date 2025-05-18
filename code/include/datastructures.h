#ifndef INCLUDE_INCLUDE_DATASTRUCTURES_H_
#define INCLUDE_INCLUDE_DATASTRUCTURES_H_

#include <cstddef>
#include <filesystem>

namespace DataStructures {

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

  static std::unique_ptr<Polyline> from_file(std::filesystem::path);

  friend std::ostream &operator<<(std::ostream &, Polyline &);
};

} // namespace DataStructures

#endif // INCLUDE_INCLUDE_DATASTRUCTURES_H_
