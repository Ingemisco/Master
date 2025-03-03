// Uses the algorithm from Van Kreveld et al. to solve the polyline
// simplification problem

#include "datastructures.h"
#include "distance.h"
#include "simplification.h"
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

namespace Simplification {

// dynamic programming table used in the algorithm
struct DPTable final {
  float *const reachable_points;
  size_t *const point_reference_i;
  size_t *const point_reference_j;

  DPTable(size_t point_count)
      : reachable_points(
            new float[point_count * point_count * (point_count - 1)]),
        point_reference_i(
            new size_t[point_count * point_count * (point_count - 1)]),
        point_reference_j(
            new size_t[point_count * point_count * (point_count - 1)]) {
    // set first layer to be unreachable, needs to be updated accordingly. All
    // other layers will be tested individually
    std::fill_n(reachable_points, point_count * (point_count - 1),
                DataStructures::UNREACHABLE);
  }

  ~DPTable() {
    delete[] reachable_points;
    delete[] point_reference_i;
    delete[] point_reference_j;
  }
};

// dim1 and dim2 should be precomputed constants with dim2 = point_count - 1,
// dim1 =  dim2 * point_count, k and i are between 0 (inclusive) and point_count
// (exclusive), j is between 0 and point_count -1
static inline float &dp_first_reachable(DPTable &table, size_t k, size_t i,
                                        size_t j, size_t dim1, size_t dim2) {
  return table.reachable_points[k * dim1 + i * dim2 + j];
}

static inline size_t &dp_point_ref_i(DPTable &table, size_t k, size_t i,
                                     size_t j, size_t dim1, size_t dim2) {
  return table.point_reference_i[k * dim1 + i * dim2 + j];
}

static inline size_t &dp_point_ref_j(DPTable &table, size_t k, size_t i,
                                     size_t j, size_t dim1, size_t dim2) {
  return table.point_reference_j[k * dim1 + i * dim2 + j];
}

Simplification
simplification_naive_euclidean(DataStructures::Polyline &polyline,
                               float epsilon) {
  size_t const point_count = polyline.point_count;

  DPTable table(point_count);
  size_t const dim2 = point_count - 1;
  size_t const dim1 = dim2 * point_count;

  // initialization (k = 0)
  DataStructures::Point const origin = polyline.get_point(0);
  float const epsilon2 = epsilon * epsilon;

  dp_first_reachable(table, 0, 0, 0, dim1, dim2) = 0;
  unsigned int j = 1;
  while (j < point_count - 1 && DataStructures::unnormalized_euclidean_distance(
                                    origin, polyline.get_point(j)) < epsilon2) {
    dp_first_reachable(table, 0, 0, j, dim1, dim2) = 0;
    j++;
  }

  for (unsigned int k = 1; true; k++) {
    for (unsigned int i = 0; i < point_count; i++) {
      for (unsigned int j = 0; j < point_count - 1; j++) {
        auto const range = DataStructures::solve_euclidean(
            polyline.get_point(j), polyline.get_point(j + 1),
            polyline.get_point(i), epsilon);
        if (range.first == DataStructures::UNREACHABLE) {
          dp_first_reachable(table, k, i, j, dim1, dim2) =
              DataStructures::UNREACHABLE;
          continue;
        } else if (range.first ==
                   dp_first_reachable(table, k - 1, i, j, dim1, dim2)) {
          // already found best, cannot be reference (would be referenced
          // earlier because of minimality) so the i, j it references dont need
          // to be set
          dp_first_reachable(table, k, i, j, dim1, dim2) = range.first;
          continue;
        }
        // compute table[k,i,j]
        float first_reachable = 2; // valid values always between 0 and 1
        size_t ref_i = -1;
        size_t ref_j = -1;
        for (unsigned int i_ = 0; i_ < i; i_++) {
          for (unsigned int j_ = 0; j_ <= j; j_++) {
            float const val =
                dp_first_reachable(table, k - 1, i_, j_, dim1, dim2);
            if (val == DataStructures::UNREACHABLE) {
              continue;
            }
            float const reachable = DataStructures::alt_godau_euclidean(
                DataStructures::PolylineRange(polyline, j_, j + 1, val),
                DataStructures::LineSegment(polyline.get_point(i_),
                                            polyline.get_point(i)),
                epsilon);
            if (reachable == DataStructures::UNREACHABLE) {
              continue;
            } else if (reachable < first_reachable) {
              first_reachable = reachable;
              ref_i = i_;
              ref_j = j_;
            }
            if (val == range.first) {
              // skip further iterations
              j_ = j + 1;
              i_ = i;
            }
          }
        }

        if (first_reachable != 2) {
          dp_point_ref_i(table, k, i, j, dim1, dim2) = ref_i;
          dp_point_ref_j(table, k, i, j, dim1, dim2) = ref_j;
          dp_first_reachable(table, k, i, j, dim1, dim2) = first_reachable;
        } else {
          dp_first_reachable(table, k, i, j, dim1, dim2) =
              DataStructures::UNREACHABLE;
        }
      }
    }

    if (dp_first_reachable(table, k, point_count - 1, point_count - 2, dim1,
                           dim2) != DataStructures::UNREACHABLE) {
      Simplification result = std::make_unique<std::vector<size_t>>(k + 1);
      auto &p = *result;
      size_t i_ = point_count - 1;
      size_t j_ = point_count - 2;
      while (i_ > 0) {
        p[k] = i_;
        size_t const i_new = dp_point_ref_i(table, k, i_, j_, dim1, dim2);
        size_t const j_new = dp_point_ref_j(table, k, i_, j_, dim1, dim2);
        k--;
        i_ = i_new;
        j_ = j_new;
      }
      p[0] = 0;
      return result;
    }
  }
}

} // namespace Simplification
