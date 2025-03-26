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
    interval_start = -candidates[i].offset / candidates[i].slope;
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
  size_t head = 0;
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
      head = j;
      last_intersection = intersection;
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
    last_intersection = intersection;
  }

  update_solution(first_sol, last_sol, list[head].offset, list[head].slope,
                  last_intersection, 1, epsilon);

  delete[] list;
  return {first_sol, last_sol};
}

// computes the product <v - u | u - w>
static inline float scalar_product(Point const &u, Point const &v,
                                   Point const &w) {
  float dot_product = 0;
  for (unsigned int i = 0; i < u.dimension; i++) {
    dot_product += (v[i] - u[i]) * (u[i] - w[i]);
  }
  return dot_product;
}

ReachabilityData solve_euclidean(Point const &point1, Point const &point2,
                                 Point const &point3, float epsilon) {
#if DEBUG
  _solver_sanity_check(point1, point2, point3);
#endif

  float const a2 =
      DataStructures::unnormalized_euclidean_distance(point1, point2);
  float const a1 = 2 * scalar_product(point1, point2, point3);
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

// decides if point 1 (true) or point 2 (false) comes before on the line
// segment. Assumes both have solutions if both same second comes first
bool solve_implicit_euclidean(LineSegment line, Point const &point1,
                              Point const &point2, float epsilon2) {
#if DEBUG
  assert_compatible_points(line.start, point1);
  assert_compatible_points(line.start, point2);
#endif
  float const point1_dist = unnormalized_euclidean_distance(point1, line.start);
  float const point2_dist = unnormalized_euclidean_distance(point2, line.start);

  float const a0_1 = point1_dist - epsilon2;
  float const a1_1 = 2 * scalar_product(line.start, line.end, point1);

  float const a0_2 = point2_dist - epsilon2;
  float const a1_2 = 2 * scalar_product(line.start, line.end, point2);

  float const a2 = unnormalized_euclidean_distance(line.start, line.end);

  float const discriminant_1 = a1_1 * a1_1 - 4 * a0_1 * a2;
  float const discriminant_2 = a1_2 * a1_2 - 4 * a0_2 * a2;

  float const x = a1_2 - a1_1;
  float const y = discriminant_1 + discriminant_2 - x * x;
  float const y2 = y * y;
  float const discr_prod = 4 * discriminant_1 * discriminant_2;

  return a0_2 > 0 && ((x <= 0 && discriminant_1 >= discriminant_2) ||
                      (x >= 0 && discriminant_1 >= discriminant_2 &&
                       (y >= 0 && discr_prod <= y2)) ||
                      (x <= 0 && discriminant_1 <= discriminant_2 &&
                       (y <= 0 || discr_prod > y2)));
}

static inline bool _is_in_01(float a1, float a2, float discriminant) {
  float const z = 2 * a2 + a1;
  return a1 <= 0 && discriminant <= a1 * a1 &&
         (z >= 0 || z * z <= discriminant);
}

// similar to the other implicit euclidean function. Returns 0 (unreachable), 1
// or
// 2. for the solutions [t01, t11], [t02, t12] on the line segment returns the
// index of the bigger of t01, t02 with same additional checks as above but also
// checks that t01 <= t12 and returns 0 if this is not the case
size_t solve_implicit_euclidean_in(LineSegment line, Point const &restriction,
                                   Point const &point, float epsilon2) {
#if DEBUG
  assert_compatible_points(line.start, restriction);
  assert_compatible_points(line.start, point);
#endif
  float const restriction_dist =
      unnormalized_euclidean_distance(restriction, line.start);

  float const point_dist = unnormalized_euclidean_distance(point, line.start);
  float const a0_1 = restriction_dist - epsilon2;
  float const a1_1 = 2 * scalar_product(line.start, line.end, restriction);

  float const a0_2 = point_dist - epsilon2;
  float const a1_2 = 2 * scalar_product(line.start, line.end, point);

  float const a2 = unnormalized_euclidean_distance(line.start, line.end);

  float const discriminant_1 = a1_1 * a1_1 - 4 * a0_1 * a2;
  float const discriminant_2 = a1_2 * a1_2 - 4 * a0_2 * a2;

  float const x = a1_1 - a1_2;
  float const y = discriminant_1 + discriminant_2 - x * x;
  float const y2 = y * y;
  float const discr_prod = 4 * discriminant_1 * discriminant_2;

  if (discriminant_2 < 0 ||
      (point_dist > epsilon2 && !_is_in_01(a1_2, a2, discriminant_2)) ||
      (x < 0 && y < 0 && discr_prod < y2)) {
    return 0;
  } else if (point_dist <= epsilon2 ||
             !((x >= 0 && discriminant_2 >= discriminant_1 &&
                (y <= 0 || y2 <= discr_prod)) ||
               (x < 0 && discriminant_2 < discriminant_1 &&
                (y >= 0 && y2 >= discr_prod)))) {
    return 1;
  }
  return 2;
}

bool is_line_reachable_euclidean(LineSegment line, Point const &point,
                                 float epsilon2) {
  float const dist_uw = unnormalized_euclidean_distance(point, line.start);
  float const dist_vw = unnormalized_euclidean_distance(point, line.end);
  if (dist_uw <= epsilon2 || dist_vw <= epsilon2) {
    return true;
  }
  float const dist_uv = unnormalized_euclidean_distance(line.start, line.end);
  float const prod = scalar_product(line.start, line.end, point);

  return 0 >= prod &&
         -prod <= unnormalized_euclidean_distance(line.start, line.end) &&
         dist_uw * dist_uv - prod * prod <= epsilon2 * dist_uv;
}

} // namespace DataStructures
