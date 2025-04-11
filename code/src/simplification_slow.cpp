// Uses the algorithm from Van Kreveld et al. to solve the polyline
// simplification problem

#include "config.h"
#include "datastructures.h"
#include "distance.h"
#include "simplification.h"
#include "visualizer.h"
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
                DataStructures::EXPLICIT_UNREACHABLE);
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

template <DataStructures::Distance _distance>
static inline void
_simplification_initialization(DataStructures::Polyline &polyline,
                               float epsilon, size_t point_count, size_t dim1,
                               size_t dim2, DPTable &table
#if DEBUG
                               ,
                               VisualizationLog::VisualizationLogger &log
#endif

) {
  DataStructures::Point const origin = polyline.get_point(0);

  dp_first_reachable(table, 0, 0, 0, dim1, dim2) = 0;
  unsigned int j = 1;
  while (j < point_count - 1 &&
         _distance(origin, polyline.get_point(j)) < epsilon) {
    dp_first_reachable(table, 0, 0, j, dim1, dim2) = 0;
    j++;
  }

#if DEBUG
  log.leave_range(j);
#endif
}

template <DataStructures::Solver _solver, DataStructures::AltGodau _alt_godau>
static inline Simplification
_simplification_main(DataStructures::Polyline &polyline, size_t point_count,
                     float epsilon, DPTable &table, size_t dim1, size_t dim2
#if DEBUG
                     ,
                     VisualizationLog::VisualizationLogger &log
#endif

) {
  for (unsigned int k = 1; true; k++) {
#pragma omp parallel for collapse(2) schedule(dynamic, 32)
    for (unsigned int i = k; i < point_count; i++) {
      for (unsigned int j = 0; j < point_count - 1; j++) {
        auto const range =
            _solver(polyline.get_point(j), polyline.get_point(j + 1),
                    polyline.get_point(i), epsilon);
        if (range.first == DataStructures::EXPLICIT_UNREACHABLE) {
          dp_first_reachable(table, k, i, j, dim1, dim2) =
              DataStructures::EXPLICIT_UNREACHABLE;
          // optimization 2: reachability
          continue;
        }

        if (range.first == dp_first_reachable(table, k - 1, i, j, dim1, dim2)) {
          // optimization 1: global minimality
          // already found best, cannot be reference (would be referenced
          // earlier because of minimality) so the i, j it references dont need
          // to be set
          dp_first_reachable(table, k, i, j, dim1, dim2) = range.first;
          continue;
        }

        // compute table[k,i,j]
        float first_reachable = 2; // valid values always between 0 and 1
        size_t ref_i = DataStructures::IMPLICIT_UNREACHABLE;
        size_t ref_j = DataStructures::IMPLICIT_UNREACHABLE;
        for (unsigned int i_ = k - 1; i_ < i; i_++) {
          for (unsigned int j_ = 0; j_ <= j; j_++) {
            float const val =
                dp_first_reachable(table, k - 1, i_, j_, dim1, dim2);
            if (val == DataStructures::EXPLICIT_UNREACHABLE) {
              continue;
            }
            float const reachable = _alt_godau(
                DataStructures::PolylineRange(polyline, j_, j + 1, val),
                DataStructures::LineSegment(polyline.get_point(i_),
                                            polyline.get_point(i)),
                epsilon);
            if (reachable == DataStructures::EXPLICIT_UNREACHABLE) {
              continue;
            } else if (reachable < first_reachable) {
              first_reachable = reachable;
              ref_i = i_;
              ref_j = j_;
              if (reachable == range.first) {
                // optimization 3: local minimality
                // skip further iterations
                j_ = j + 1;
                i_ = i;
              }
            }
          }
        }

        if (first_reachable != 2) {
          dp_point_ref_i(table, k, i, j, dim1, dim2) = ref_i;
          dp_point_ref_j(table, k, i, j, dim1, dim2) = ref_j;
          dp_first_reachable(table, k, i, j, dim1, dim2) = first_reachable;
        } else {
          dp_first_reachable(table, k, i, j, dim1, dim2) =
              DataStructures::EXPLICIT_UNREACHABLE;
        }
      }
    }

    if (dp_first_reachable(table, k, point_count - 1, point_count - 2, dim1,
                           dim2) != DataStructures::EXPLICIT_UNREACHABLE) {
#if DEBUG
      // print whole table (relevant entries)
      for (unsigned int a = 1; a <= k; a++) {
        for (unsigned int b = a; b < point_count; b++) {
          for (unsigned int c = 0; c < point_count - 1; c++) {
            float v = dp_first_reachable(table, a, b, c, dim1, dim2);
            if (v != DataStructures::UNREACHABLE &&
                (a == 0 ||
                 dp_first_reachable(table, a - 1, b, c, dim1, dim2) != v)) {
              size_t i_ = dp_point_ref_i(table, a, b, c, dim1, dim2);
              size_t j_ = dp_point_ref_j(table, a, b, c, dim1, dim2);
              log.add_use(
                  VisualizationLog::VisualizationData(a, b, c, i_, j_, v));
            }
          }
        }
      }
      log.emit();
#endif

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

Simplification
simplification_naive_euclidean(DataStructures::Polyline &polyline,
                               float epsilon) {
  size_t const point_count = polyline.point_count;

  DPTable table(point_count);
  size_t const dim2 = point_count - 1;
  size_t const dim1 = dim2 * point_count;
  float const epsilon2 = epsilon * epsilon;

#if DEBUG
  VisualizationLog::VisualizationLogger log(
      polyline, epsilon, VisualizationLog::Distance::EUCLIDEAN);
  _simplification_initialization<
      DataStructures::unnormalized_euclidean_distance>(
      polyline, epsilon2, point_count, dim1, dim2, table, log);
  return _simplification_main<DataStructures::solve_euclidean,
                              DataStructures::alt_godau_euclidean>(
      polyline, point_count, epsilon, table, dim1, dim2, log);
#else
  _simplification_initialization<
      DataStructures::unnormalized_euclidean_distance>(
      polyline, epsilon2, point_count, dim1, dim2, table);
  return _simplification_main<DataStructures::solve_euclidean,
                              DataStructures::alt_godau_euclidean>(
      polyline, point_count, epsilon, table, dim1, dim2);
#endif
}

Simplification
simplification_naive_manhattan(DataStructures::Polyline &polyline,
                               float epsilon) {
  size_t const point_count = polyline.point_count;

  DPTable table(point_count);
  size_t const dim2 = point_count - 1;
  size_t const dim1 = dim2 * point_count;

#if DEBUG
  VisualizationLog::VisualizationLogger log(
      polyline, epsilon, VisualizationLog::Distance::MANHATTAN);
  _simplification_initialization<DataStructures::manhattan_distance>(
      polyline, epsilon, point_count, dim1, dim2, table, log);
  return _simplification_main<DataStructures::solve_manhattan,
                              DataStructures::alt_godau_manhattan>(
      polyline, point_count, epsilon, table, dim1, dim2, log);
#else
  _simplification_initialization<DataStructures::manhattan_distance>(
      polyline, epsilon, point_count, dim1, dim2, table);
  return _simplification_main<DataStructures::solve_manhattan,
                              DataStructures::alt_godau_manhattan>(
      polyline, point_count, epsilon, table, dim1, dim2);
#endif
}

Simplification
simplification_naive_chebyshev(DataStructures::Polyline &polyline,
                               float epsilon) {
  size_t const point_count = polyline.point_count;

  DPTable table(point_count);
  size_t const dim2 = point_count - 1;
  size_t const dim1 = dim2 * point_count;

#if DEBUG
  VisualizationLog::VisualizationLogger log(
      polyline, epsilon, VisualizationLog::Distance::CHEBYSHEV);
  _simplification_initialization<DataStructures::chebyshev_distance>(
      polyline, epsilon, point_count, dim1, dim2, table, log);
  return _simplification_main<DataStructures::solve_chebyshev,
                              DataStructures::alt_godau_chebyshev>(
      polyline, point_count, epsilon, table, dim1, dim2, log);
#else
  _simplification_initialization<DataStructures::chebyshev_distance>(
      polyline, epsilon, point_count, dim1, dim2, table);
  return _simplification_main<DataStructures::solve_chebyshev,
                              DataStructures::alt_godau_chebyshev>(
      polyline, point_count, epsilon, table, dim1, dim2);
#endif
}

// dynamic programming table used in the algorithm
struct DPImplicitTable final {
  size_t *const restriction;
  size_t *const point_reference_i;
  size_t *const point_reference_j;

  DPImplicitTable(size_t point_count)
      : restriction(new size_t[point_count * point_count * (point_count - 1)]),
        point_reference_i(
            new size_t[point_count * point_count * (point_count - 1)]),
        point_reference_j(
            new size_t[point_count * point_count * (point_count - 1)]) { }

  ~DPImplicitTable() {
    delete[] restriction;
    delete[] point_reference_i;
    delete[] point_reference_j;
  }
};

static inline size_t &dp_restriction(DPImplicitTable &table, size_t k, size_t i,
                                     size_t j, size_t dim1, size_t dim2) {
  return table.restriction[k * dim1 + i * dim2 + j];
}

static inline size_t &dp_point_ref_i(DPImplicitTable &table, size_t k, size_t i,
                                     size_t j, size_t dim1, size_t dim2) {
  return table.point_reference_i[k * dim1 + i * dim2 + j];
}

static inline size_t &dp_point_ref_j(DPImplicitTable &table, size_t k, size_t i,
                                     size_t j, size_t dim1, size_t dim2) {
  return table.point_reference_j[k * dim1 + i * dim2 + j];
}

Simplification
simplification_naive_euclidean_implicit(DataStructures::Polyline &polyline,
                                        float epsilon) {
  size_t const point_count = polyline.point_count;
  DPImplicitTable table(point_count);
  size_t const dim2 = point_count - 1;
  size_t const dim1 = dim2 * point_count;
  float const epsilon2 = epsilon * epsilon;

  // initialization
  auto const origin = polyline.get_point(0);
  dp_restriction(table, 0, 0, 0, dim1, dim2) = 0;
  unsigned int j = 1;
  for (; j < point_count - 1 && DataStructures::unnormalized_euclidean_distance(
                                origin, polyline.get_point(j)) < epsilon2;
       j++) {
    dp_restriction(table, 0, 0, j, dim1, dim2) = 0;
  }

#if DEBUG
  VisualizationLog::VisualizationLogger log(
      polyline, epsilon, VisualizationLog::Distance::EUCLIDEAN_IMPLICIT);
  log.leave_range(j);
#endif
  for (unsigned int i = 0; i < point_count; i++) {
		for (; j < point_count - 1; j++) {
			if (DataStructures::is_line_reachable_euclidean(DataStructures::LineSegment(polyline.get_point(j), polyline.get_point(j + 1)), polyline.get_point(i), epsilon2)) {
				dp_restriction(table, 0, i, j, dim1, dim2) = DataStructures::IMPLICIT_UNREACHABLE;
			} else {
				dp_restriction(table, 0, i, j, dim1, dim2) = DataStructures::IMPLICIT_NEVER_REACHABLE;
			}
		}
		j = 0;
	}

	std::cout << "Test" << std::endl;

  for (unsigned int k = 1; true; k++) {
#pragma omp parallel for collapse(2) schedule(dynamic, 32)
    for (size_t i = k; i < point_count; i++) {
      for (size_t j = 0; j < point_count - 1; j++) {
        if (dp_restriction(table, k - 1, i, j, dim1, dim2) == i) {
          // optimization 1: global minimality
          dp_restriction(table, k, i, j, dim1, dim2) = i;
          continue;
        } else if (dp_restriction(table, k - 1, i, j, dim1, dim2) == DataStructures::IMPLICIT_NEVER_REACHABLE) {
          dp_restriction(table, k, i, j, dim1, dim2) = DataStructures::IMPLICIT_NEVER_REACHABLE;
          continue;
        }

        size_t new_restriction = DataStructures::IMPLICIT_UNREACHABLE;
        size_t ref_i = DataStructures::IMPLICIT_UNREACHABLE;
        size_t ref_j = DataStructures::IMPLICIT_UNREACHABLE;
        for (size_t i_ = k - 1; i_ < i; i_++) {
          for (size_t j_ = 0; j_ <= j; j_++) {
            size_t const val = dp_restriction(table, k - 1, i_, j_, dim1, dim2);

            if (val >= DataStructures::IMPLICIT_NEVER_REACHABLE) {
              continue;
            }

            auto const computed_restriction =
                DataStructures::alt_godau_euclidean_implicit(
                    polyline, j_, j, i_, i, val, epsilon2);

            if (computed_restriction == DataStructures::IMPLICIT_UNREACHABLE) {
              continue;
            } else if (computed_restriction == i) {
              // optimization 3: local minimality
              new_restriction = computed_restriction;
              ref_i = i_;
              ref_j = j_;
              j_ = j + 1;
              i_ = i;
            } else if (new_restriction == DataStructures::IMPLICIT_UNREACHABLE ||
                       !DataStructures::solve_implicit_euclidean(
                           DataStructures::LineSegment(
                               polyline.get_point(j),
                               polyline.get_point(j + 1)),
                           polyline.get_point(new_restriction),
                           polyline.get_point(computed_restriction),
                           epsilon2)) {
              new_restriction = computed_restriction;
              ref_i = i_;
              ref_j = j_;
            }
          }
        }

        dp_point_ref_i(table, k, i, j, dim1, dim2) = ref_i;
        dp_point_ref_j(table, k, i, j, dim1, dim2) = ref_j;
        dp_restriction(table, k, i, j, dim1, dim2) = new_restriction;
      }
    }

    if (dp_restriction(table, k, point_count - 1, point_count - 2, dim1,
                       dim2) != DataStructures::IMPLICIT_UNREACHABLE) {
#if DEBUG
      // print whole table (relevant entries)
      for (unsigned int a = 1; a <= k; a++) {
        for (unsigned int b = a; b < point_count; b++) {
          for (unsigned int c = 0; c < point_count - 1; c++) {
            size_t v = dp_restriction(table, a, b, c, dim1, dim2);
            if (v < DataStructures::IMPLICIT_NEVER_REACHABLE &&
                (a == 0 ||
                 dp_restriction(table, a - 1, b, c, dim1, dim2) != v)) {

              size_t i_ = dp_point_ref_i(table, a, b, c, dim1, dim2);
              size_t j_ = dp_point_ref_j(table, a, b, c, dim1, dim2);
              log.add_use(
                  VisualizationLog::VisualizationData(a, b, c, i_, j_, v));
            }
          }
        }
      }
      log.emit();
#endif

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
