#include "config.h"
#include "datastructures.h"
#include "distance.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace DataStructures {

static inline void _solver_sanity_check(Point const &point1,
                                        Point const &point2,
                                        Point const &point3) {
  assert_compatible_points(point1, point2);
  assert_compatible_points(point1, point3);
  if (manhattan_distance(point1, point2) == 0) {
    throw std::runtime_error(
        "The two points given as a line segment are the same point.");
  }
}

typedef struct {
  float offset;
  float slope;
} CoordCandidate;

// returns true if no further computation would change result
static inline bool make_solution(float &first_sol, float &last_sol,
                                 float new_sol_first, float new_sol_last) {
  if (first_sol == UNREACHABLE) {
    first_sol = new_sol_first;
    if (last_sol == UNREACHABLE) {
      last_sol = first_sol;
    } else {
      return true;
    }
  } else {
    last_sol = new_sol_last;
  }

  return false;
}

// if true returned, solution already determined
static inline bool update_solution(float &first_sol, float &last_sol, float a,
                                   float b, float range_start, float range_end,
                                   float epsilon) {
  if (b == 0) {
    return a == epsilon &&
           make_solution(first_sol, last_sol, range_start, range_end);
  }

  float const solution = (epsilon - a) / b;
  return solution >= range_start && solution <= range_end &&
         make_solution(first_sol, last_sol, solution, solution);
}

ReachabilityData solve_manhattan(Point const &point1, Point const &point2,
                                 Point const &point3, float epsilon) {
#if DEBUG
  _solver_sanity_check(point1, point2, point3);
#endif

  bool const first_reachable = manhattan_distance(point1, point3) < epsilon;
  bool const last_reachable = manhattan_distance(point2, point3) < epsilon;
  if (first_reachable && last_reachable) {
    return {0, 1};
  }

  float first_sol = first_reachable ? 0 : UNREACHABLE;
  float last_sol = last_reachable ? 1 : UNREACHABLE;

  CoordCandidate *candidates = new CoordCandidate[point1.dimension];
  unsigned int index = 0;
  float global_offset = 0;
  float global_slope = 0;
  for (unsigned int i = 0; i < point1.dimension; i++) {
    float slope = point2[i] - point1[i];
    float offset = point1[i] - point3[i];
    if (slope < 0) {
      slope *= -1;
      offset *= -1;
    } else if (slope == 0) {
      epsilon -= std::abs(offset);
      continue;
    }
    float const zero = -offset / slope;
    // if the zero is outside of [0,1] the value will never change in the
    // algorithm so we don't need to sort the values to find the intervals
    if (0 >= zero) {
      global_offset += offset;
      global_slope += slope;
      continue;
    } else if (1 <= zero) {
      global_offset -= offset;
      global_slope -= slope;
      continue;
    }
    candidates[index].slope = slope;
    candidates[index].offset = offset;

    // we start at zero, so all values (which are in  [0,1]) are positive
    global_offset -= offset;
    global_slope -= slope;

    index++;
  }

  // sort according to zero of the line segments
  std::sort(candidates, candidates + index,
            [](CoordCandidate a, CoordCandidate b) {
              return -a.offset / a.slope < -b.offset / b.slope;
            });

  float interval_start = 0;
  for (unsigned int i = 0; i < index; i++) {
    if (update_solution(first_sol, last_sol, global_offset, global_slope,
                        interval_start,
                        -candidates[i].offset / candidates[i].slope, epsilon)) {
      delete[] candidates;
      return {first_sol, last_sol};
    }

    global_offset += 2 * candidates[i].offset;
    global_slope += 2 * candidates[i].slope;
  }

  update_solution(first_sol, last_sol, global_offset, global_slope,
                  interval_start, 1, epsilon);

  delete[] candidates;
  return {first_sol, last_sol};
}

typedef struct {
  int prev;
  int next;
  float offset;
  float slope;
} ChebyshevList;

typedef struct chebyshev_queue_elem {
  float intersection;
  int prev_above;
  int prev_below;

  bool operator>(struct chebyshev_queue_elem const &other) const {
    return this->intersection > other.intersection;
  }
} ChebyshevQueueElement;

#define LINE_REMOVED -1

ReachabilityData solve_chebyshev(Point const &point1, Point const &point2,
                                 Point const &point3, float epsilon) {
#if DEBUG
  _solver_sanity_check(point1, point2, point3);
#endif
  bool const first_reachable = chebyshev_distance(point1, point3) < epsilon;
  bool const last_reachable = chebyshev_distance(point2, point3) < epsilon;
  if (first_reachable && last_reachable) {
    return {0, 1};
  }

  float first_sol = first_reachable ? 0 : UNREACHABLE;
  float last_sol = last_reachable ? 1 : UNREACHABLE;

  CoordCandidate *candidates = new CoordCandidate[2 * point1.dimension];
  unsigned int index = 0;
  for (unsigned int i = 0; i < point1.dimension; i++) {
    float const offset = point1[i] - point3[i];
    float const slope = point2[i] - point1[i];
    if (!(offset <= 0 && slope + offset <= 0)) {
      candidates[index++] = {offset, slope};
    }
    if (!(offset >= 0 && slope + offset >= 0)) {
      candidates[index++] = {-offset, -slope};
    }
  }

  std::sort(
      candidates, candidates + index,
      [](CoordCandidate a, CoordCandidate b) { return a.offset > b.offset; });
  ChebyshevList *list = new ChebyshevList[index];

  std::vector<ChebyshevQueueElement> inner_vector;
  inner_vector.reserve(index);
  std::priority_queue<ChebyshevQueueElement, std::vector<ChebyshevQueueElement>,
                      std::greater<ChebyshevQueueElement>>
      queue(std::greater<ChebyshevQueueElement>(), std::move(inner_vector));

  list[0].offset = candidates[0].offset;
  list[0].slope = candidates[0].slope;
  list[0].prev = 0; // head always has as prev itself
  size_t current = 0;
  for (unsigned int i = 1; i < index; i++) {
    float const a_ = list[current].offset;
    float const b_ = list[current].slope;
    float const a = candidates[i].offset;
    float const b = candidates[i].slope;
    if (a_ + b_ >= a + b) {
      continue;
    }

    list[current].next = current + 1;
    list[current + 1].prev = current;
    list[current + 1].offset = a;
    list[current + 1].slope = b;

    float const intersection = (a_ - a) / (b - b_);
    queue.emplace(intersection, current, current + 1);
    current++;
  }
  delete[] candidates;

  float last_intersection = 0;
  while (!queue.empty()) {
    auto data = queue.top();
    queue.pop();
    int const i = data.prev_above;
    int const j = data.prev_below;
    float const intersection = data.intersection;
    if (list[i].prev == LINE_REMOVED || list[j].next == LINE_REMOVED) {
      continue;
    } else if (list[i].prev == i) { // above is head so need to set new head and
                                    // test for solution
      list[j].prev = j;             // make new head
      list[i].next = LINE_REMOVED;
      if (update_solution(first_sol, last_sol, list[i].offset, list[i].slope,
                          last_intersection, intersection, epsilon)) {
        delete[] list;
        return {first_sol, last_sol};
      }
      continue;
    }

    int const before = list[i].prev;
    list[before].next = j;
    list[j].prev = before;
    list[i].prev = LINE_REMOVED;
    if (list[before].prev != before) {
      float const a = list[j].offset;
      float const b = list[j].slope;
      float const a_ = list[before].offset;
      float const b_ = list[before].slope;
      float const intersection = (a_ - a) / (b - b_);
      queue.emplace(intersection, before, j);
    }
  }

  update_solution(first_sol, last_sol, list[0].offset, list[0].slope,
                  last_intersection, 1, epsilon);

  delete[] list;
  return {first_sol, last_sol};
}

ReachabilityData solve_euclidean(Point const &point1, Point const &point2,
                                 Point const &point3, float epsilon) {
#if DEBUG
  _solver_sanity_check(point1, point2, point3);
#endif
  float dot_product = 0;
  for (unsigned int i = 0; i < point1.dimension; i++) {
    dot_product += (point2[i] - point1[i]) * (point1[i] - point3[i]);
  }

  float const a2 =
      DataStructures::unnormalized_euclidean_distance(point1, point2);
  float const a1 = 2 * dot_product;
  float const a0 =
      DataStructures::unnormalized_euclidean_distance(point1, point3) -
      epsilon * epsilon;

  float const discriminant = a1 * a1 - 4 * a2 * a0;
  if (discriminant < 0) {
    return {UNREACHABLE, UNREACHABLE};
  }

  float const root = std::sqrt(discriminant);
  float const t0 = (-a1 - root) / (2 * a2);
  float const t1 = (-a1 + root) / (2 * a2);
  if (t0 > 1 || t1 < 0) {
    return {UNREACHABLE, UNREACHABLE};
  }
  return {t0 < 0 ? 0 : t0, t1 > 1 ? 1 : t1};
}
} // namespace DataStructures
