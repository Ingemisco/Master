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
#include <utility>
#include <vector>

namespace Simplification {

// dynamic programming table used in the algorithm
template <typename F, typename T, F UNREACHABLE_POINT>
struct DPTable final {
  F *const reachable_points;
  size_t *const point_reference_i;
  size_t *const point_reference_j;

	size_t const dim2;
	size_t const dim1;

  DPTable(size_t point_count)
      : reachable_points( new F[point_count * point_count * (point_count - 1)]), 
				point_reference_i( new size_t[point_count * point_count * (point_count - 1)]),
        point_reference_j( new size_t[point_count * point_count * (point_count - 1)]),
				dim2(point_count - 1), dim1(dim2 * point_count) {
    // set first layer to be unreachable, needs to be updated accordingly. All
    // other layers will be tested individually
    std::fill_n(reachable_points, point_count * (point_count - 1), UNREACHABLE_POINT);
  }

  ~DPTable() {
    delete[] reachable_points;
    delete[] point_reference_i;
    delete[] point_reference_j;
  }

	inline F &first_reachable(size_t k, size_t i, size_t j) {
		return reachable_points[k * dim1 + i * dim2 + j];
	}

	inline size_t &dp_point_ref_i(size_t k, size_t i, size_t j) {
		return point_reference_i[k * dim1 + i * dim2 + j];
	}

	inline size_t &dp_point_ref_j(size_t k, size_t i, size_t j) {
		return point_reference_j[k * dim1 + i * dim2 + j];
	}
};

typedef DPTable<float, float, DataStructures::EXPLICIT_UNREACHABLE> DPExplicit;
typedef DPTable<DataStructures::FRValue, DataStructures::LRValue, DataStructures::SEMIEXPLICIT_UNREACHABLE> DPSemiExplicit;

template <typename F, typename L, std::pair<F, L> empty_interval, F start, DataStructures::Distance _distance>
static inline void _simplification_initialization(DataStructures::Polyline &polyline, float epsilon, size_t point_count, DPTable<F, L, empty_interval.first> &table
#if DEBUG
                               , VisualizationLog::VisualizationLogger &log
#endif
) {
  table.first_reachable(0, 0, 0) = start;
  unsigned int j = 1;
  for (; j < point_count - 1 && _distance(polyline, 0, j) < epsilon; j++) {
    table.first_reachable(0, 0, j) = start;
  }

#if DEBUG
  log.leave_range(j);
#endif
}

template <typename F, typename L, std::pair<F, L> empty_interval, 
	std::pair<F, L> _solver(DataStructures::Polyline const &, size_t, size_t, size_t, float),
	F _alt_godau(DataStructures::Polyline const &, size_t, size_t, F, size_t, size_t, float)>
static inline Simplification
_simplification_main(DataStructures::Polyline &polyline, size_t point_count, float epsilon, DPTable<F, L, empty_interval.first> &table
#if DEBUG
                     , VisualizationLog::VisualizationLogger &log
#endif
) {
  for (unsigned int k = 1; true; k++) {
#pragma omp parallel for collapse(2) schedule(dynamic, 32)
    for (unsigned int i = k; i < point_count; i++) {
      for (unsigned int j = 0; j < point_count - 1; j++) {
        auto const range = _solver(polyline, j, j + 1, i, epsilon);
        if (range == empty_interval) {
          table.first_reachable(k, i, j) = empty_interval.first;
          // optimization 2: reachability
          continue;
        }

        if (range.first == table.first_reachable(k - 1, i, j)) {
          // optimization 1: global minimality
          // already found best, cannot be reference (would be referenced
          // earlier because of minimality) so the i, j it references dont need
          // to be set
          table.first_reachable(k, i, j) = range.first;
          continue;
        }

        // compute table[k,i,j]
        F first_reachable = empty_interval.first;
        size_t ref_i = DataStructures::IMPLICIT_UNREACHABLE;
        size_t ref_j = DataStructures::IMPLICIT_UNREACHABLE;
        for (unsigned int i_ = k - 1; i_ < i; i_++) {
          for (unsigned int j_ = 0; j_ <= j; j_++) {
            F const val = table.first_reachable(k - 1, i_, j_);
            if (val == empty_interval.first) {
              continue;
            }
            F const reachable = _alt_godau(polyline, j_, j + 1, val, i_, i, epsilon);
            if (reachable == empty_interval.first) {
              continue;
            } else if (first_reachable == empty_interval.first || reachable < first_reachable) {
              first_reachable = reachable;
              ref_i = i_;
              ref_j = j_;
              if (reachable == range.first) {
                // optimization 3: local minimality
								goto local_minimality_skip;
              }
            }
          }
        }
local_minimality_skip:
        if (first_reachable != empty_interval.first) {
          table.dp_point_ref_i(k, i, j) = ref_i;
          table.dp_point_ref_j(k, i, j) = ref_j;
          table.first_reachable(k, i, j) = first_reachable;
        } else {
          table.first_reachable(k, i, j) = empty_interval.first;
        }
      }
    }

    if (table.first_reachable(k, point_count - 1, point_count - 2) != empty_interval.first) {
#if DEBUG
      // print whole table (relevant entries)
      for (unsigned int a = 1; a <= k; a++) {
        for (unsigned int b = a; b < point_count; b++) {
          for (unsigned int c = 0; c < point_count - 1; c++) {
            float v = table.first_reachable(a, b, c);
            if (v != DataStructures::EXPLICIT_UNREACHABLE &&
                (a == 0 || table.first_reachable(a - 1, b, c) != v)) {
              size_t i_ = table.dp_point_ref_i(a, b, c);
              size_t j_ = table.dp_point_ref_j(a, b, c);
              log.add_use( VisualizationLog::VisualizationData(a, b, c, i_, j_, v));
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
        size_t const i_new = table.dp_point_ref_i(k, i_, j_);
        size_t const j_new = table.dp_point_ref_j(k, i_, j_);
        k--;
        i_ = i_new;
        j_ = j_new;
      }
      p[0] = 0;
      return result;
    }
  }
}

Simplification simplification_naive_euclidean(DataStructures::Polyline &polyline, float epsilon) {
  size_t const point_count = polyline.point_count;

  DPExplicit table(point_count);
  float const epsilon2 = epsilon * epsilon;

#if DEBUG
  VisualizationLog::VisualizationLogger log(
      polyline, epsilon, VisualizationLog::Distance::EUCLIDEAN);
  _simplification_initialization<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT,
      DataStructures::unnormalized_euclidean_distance>(
      polyline, epsilon2, point_count, table, log);
  return _simplification_main<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT, DataStructures::solve_euclidean,
                              DataStructures::alt_godau_euclidean>(
      polyline, point_count, epsilon, table, log);
#else
  _simplification_initialization<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT, 0.0f,
      DataStructures::unnormalized_euclidean_distance>(
      polyline, epsilon2, point_count, table);
  return _simplification_main<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT, DataStructures::solve_euclidean,
                              DataStructures::alt_godau_euclidean>(
      polyline, point_count, epsilon, table);
#endif
}

Simplification simplification_naive_manhattan(DataStructures::Polyline &polyline, float epsilon) {
  size_t const point_count = polyline.point_count;
  DPExplicit table(point_count);

#if DEBUG
  VisualizationLog::VisualizationLogger log(
      polyline, epsilon, VisualizationLog::Distance::MANHATTAN);
  _simplification_initialization<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT,
DataStructures::manhattan_distance>(
      polyline, epsilon, point_count, table, log);
  return _simplification_main<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT, DataStructures::solve_manhattan,
                              DataStructures::alt_godau_manhattan>(
      polyline, point_count, epsilon, table, log);
#else
  _simplification_initialization<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT, 0.0f,
DataStructures::manhattan_distance>(
      polyline, epsilon, point_count, table);
  return _simplification_main<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT, DataStructures::solve_manhattan,
                              DataStructures::alt_godau_manhattan>(
      polyline, point_count, epsilon, table);
#endif
}

Simplification simplification_naive_chebyshev(DataStructures::Polyline &polyline, float epsilon) {
  size_t const point_count = polyline.point_count;
  DPExplicit table(point_count);

#if DEBUG
  VisualizationLog::VisualizationLogger log(
      polyline, epsilon, VisualizationLog::Distance::CHEBYSHEV);
  _simplification_initialization<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT,
DataStructures::chebyshev_distance>(polyline, epsilon, point_count, table, log);
  return _simplification_main<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT, DataStructures::solve_chebyshev,
                              DataStructures::alt_godau_chebyshev>(polyline, point_count, epsilon, table, log);
#else
  _simplification_initialization<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT, 0.0f,
DataStructures::chebyshev_distance>(polyline, epsilon, point_count, table);
  return _simplification_main<float, float, DataStructures::EMPTY_INTERVAL_EXPLICIT, DataStructures::solve_chebyshev,
                              DataStructures::alt_godau_chebyshev>(polyline, point_count, epsilon, table);
#endif
}

// dynamic programming table used in the algorithm
struct DPImplicitTable final {
  size_t *const restriction;
  size_t *const point_reference_i;
  size_t *const point_reference_j;

	size_t const dim2;
	size_t const dim1;

  DPImplicitTable(size_t point_count)
      : restriction(new size_t[point_count * point_count * (point_count - 1)]),
        point_reference_i(new size_t[point_count * point_count * (point_count - 1)]),
        point_reference_j(new size_t[point_count * point_count * (point_count - 1)]),
				dim2(point_count - 1), dim1(dim2 * point_count) { }

  ~DPImplicitTable() {
    delete[] restriction;
    delete[] point_reference_i;
    delete[] point_reference_j;
  }

	inline size_t &dp_restriction(size_t k, size_t i, size_t j) {
		return restriction[k * dim1 + i * dim2 + j];
	}

	inline size_t &dp_point_ref_i(size_t k, size_t i, size_t j) {
		return point_reference_i[k * dim1 + i * dim2 + j];
	}

	inline size_t &dp_point_ref_j(size_t k, size_t i, size_t j) {
		return point_reference_j[k * dim1 + i * dim2 + j];
	}
};

Simplification simplification_naive_euclidean_implicit(DataStructures::Polyline &polyline, float epsilon) {
  size_t const point_count = polyline.point_count;
  DPImplicitTable table(point_count);
  float const epsilon2 = epsilon * epsilon;

  // initialization
  table.dp_restriction(0, 0, 0) = 0;
  unsigned int j = 1;
  for (; j < point_count - 1 && DataStructures::unnormalized_euclidean_distance(polyline, 0, j) < epsilon2; j++) {
    table.dp_restriction(0, 0, j) = 0;
  }

#if DEBUG
  VisualizationLog::VisualizationLogger log(
      polyline, epsilon, VisualizationLog::Distance::EUCLIDEAN_IMPLICIT);
  log.leave_range(j);
#endif
  for (unsigned int i = 0; i < point_count; i++) {
		for (; j < point_count - 1; j++) {
			if (DataStructures::is_line_reachable_euclidean(polyline, j, j + 1, i, epsilon2)) {
				table.dp_restriction(0, i, j) = DataStructures::IMPLICIT_UNREACHABLE;
			} else {
				table.dp_restriction(0, i, j) = DataStructures::IMPLICIT_NEVER_REACHABLE;
			}
		}
		j = 0;
	}

  for (unsigned int k = 1; true; k++) {
#pragma omp parallel for collapse(2) schedule(dynamic, 32)
    for (size_t i = k; i < point_count; i++) {
      for (size_t j = 0; j < point_count - 1; j++) {
        if (table.dp_restriction(k - 1, i, j) == i) {
          // optimization 1: global minimality
          table.dp_restriction(k, i, j) = i;
          continue;
        } else if (table.dp_restriction(k - 1, i, j) == DataStructures::IMPLICIT_NEVER_REACHABLE) {
          table.dp_restriction(k, i, j) = DataStructures::IMPLICIT_NEVER_REACHABLE;
          continue;
        }

        size_t new_restriction = DataStructures::IMPLICIT_UNREACHABLE;
        size_t ref_i = DataStructures::IMPLICIT_UNREACHABLE;
        size_t ref_j = DataStructures::IMPLICIT_UNREACHABLE;
        for (size_t i_ = k - 1; i_ < i; i_++) {
          for (size_t j_ = 0; j_ <= j; j_++) {
            size_t const val = table.dp_restriction(k - 1, i_, j_);

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
							goto local_minimality_skip;
            } else if (new_restriction == DataStructures::IMPLICIT_UNREACHABLE ||
								 !DataStructures::solve_implicit_euclidean(polyline, j, j + 1, new_restriction, computed_restriction, epsilon2)) {
              new_restriction = computed_restriction;
              ref_i = i_;
              ref_j = j_;
            }
          }
        }
local_minimality_skip:
        table.dp_point_ref_i(k, i, j) = ref_i;
        table.dp_point_ref_j(k, i, j) = ref_j;
        table.dp_restriction(k, i, j) = new_restriction;
      }
    }

    if (table.dp_restriction(k, point_count - 1, point_count - 2) != DataStructures::IMPLICIT_UNREACHABLE) {
#if DEBUG
      // print whole table (relevant entries)
      for (unsigned int a = 1; a <= k; a++) {
        for (unsigned int b = a; b < point_count; b++) {
          for (unsigned int c = 0; c < point_count - 1; c++) {
            size_t v = table.dp_restriction(a, b, c);
            if (v < DataStructures::IMPLICIT_NEVER_REACHABLE &&
                (a == 0 || table.dp_restriction(a - 1, b, c) != v)) {
              size_t i_ = table.dp_point_ref_i(a, b, c);
              size_t j_ = table.dp_point_ref_j(a, b, c);
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
        size_t const i_new = table.dp_point_ref_i(k, i_, j_);
        size_t const j_new = table.dp_point_ref_j(k, i_, j_);
        k--;
        i_ = i_new;
        j_ = j_new;
      }
      p[0] = 0;

      return result;
    }
  }
}








Simplification simplification_naive_euclidean_semiexplicit(DataStructures::Polyline &polyline, float epsilon) {
  size_t const point_count = polyline.point_count;

  DPSemiExplicit table(point_count);
	std::cout << epsilon << std::endl;

#if DEBUG
  VisualizationLog::VisualizationLogger log(polyline, epsilon, VisualizationLog::Distance::EUCLIDEAN);
  _simplification_initialization<DataStructures::FRValue, DataStructures::LRValue, DataStructures::EMPTY_INTERVAL_SEMIEXPLICIT, DataStructures::unnormalized_euclidean_distance>(polyline, epsilon, point_count, table, log);
  return _simplification_main<FRValue, LRValue, DataStructures::EMPTY_INTERVAL_SEMIEXPLICIT, DataStructures::solve_euclidean, DataStructures::alt_godau_euclidean>( polyline, point_count, epsilon, table, log);
#else
  _simplification_initialization<DataStructures::FRValue, DataStructures::LRValue, DataStructures::EMPTY_INTERVAL_SEMIEXPLICIT, DataStructures::FRValue(0,0), DataStructures::unnormalized_euclidean_distance>(polyline, epsilon, point_count, table);
  return _simplification_main<DataStructures::FRValue, DataStructures::LRValue, DataStructures::EMPTY_INTERVAL_SEMIEXPLICIT, DataStructures::solve_euclidean_se, DataStructures::alt_godau_euclidean_semiexplicit>(polyline, point_count, epsilon, table);
#endif
}








} // namespace Simplification
