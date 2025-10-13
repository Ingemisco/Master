#include <cassert>
#include <cstddef>
#include <generator>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>
#include <tuple>

#include "datastructures.h"
#include "distance.h"
#include "global.h"
#include "log.h"
#include "simplification.h"

using DataStructures::Polyline;
using DataStructures::Queue;
using DataStructures::ReachabilityData;

namespace Simplification {

struct Interval final {
	PointIndex start_vertex; // start vertex of line segment on which interval starts
	PointIndex end_vertex;   // start vertex of line segment on which interval ends
	float rel_interval_start;
	float rel_interval_end;
};

// contains possible mapping of one interval to another (I times_<= J in thesis)
typedef std::pair<Interval, Interval> IntervalMapping;

struct SolutionIntervals {
	std::vector<Interval> intervals;
};

struct ShortcutMappings {
	std::vector<IntervalMapping> shortcut_mappings;
};

typedef std::tuple<size_t, size_t, size_t> RangeInterval; // first value is index of domain interval, second and third are first and last index of range interval


	struct QueueData final {
		PointIndex first_interval_index;
		PointIndex last_interval_index;
		float first_reachable;

	};

		std::ostream& operator<<(std::ostream& os, const QueueData& qd) {
				os << "("<< qd.first_interval_index << ", "
					 << qd.last_interval_index << ", "
					 << qd.first_reachable << ")";
				return os;
}


struct GlobalShortcutGraph final {
	Polyline const &polyline;
	PointCount const point_count;
	float const epsilon;
	SolutionIntervals *solution_intervals;
	// ShortcutMappings *global_shortcut_mappings;


	GlobalShortcutGraph(Polyline const &polyline, float epsilon) : polyline(polyline), point_count(polyline.point_count), epsilon(epsilon) {
		this->solution_intervals = new SolutionIntervals[polyline.point_count];
		// this->global_shortcut_mappings = new ShortcutMappings[polyline.point_count * (polyline.point_count - 1)];

		this->initialize_solution_intervals();
		// this->initialize_shortcut_mappings(polyline); 
	}

	~GlobalShortcutGraph() {
		delete[] this->solution_intervals;
		// delete[] this->global_shortcut_mappings;
	}

	// inline ShortcutMappings const &interval_mappings(PointIndex start, PointIndex end) const {
	// 	PointIndex const index = (this->point_count - 1) * start - start * (start + 1) / 2 + end - 1;
	// 	return this->global_shortcut_mappings[index];
	// }

	inline void initialize_solution_intervals() {
		for (PointCount i = 0u; i < this->point_count; i++) {
			for (auto j = 0u; j < this->point_count - 1; j++) {
				auto interval = DataStructures::solve_euclidean(polyline, j, j + 1, i, this->epsilon);
				if (interval == DataStructures::EMPTY_INTERVAL_EXPLICIT) {
					continue;
				}
				
				auto const interval_start = interval.first;
				auto const vertex_start = j;
				auto interval_end = interval.second;

				// merge adjacent intervals to one major interval.
				// The solutions here can never be empty because 0 is valid solution since 1 is valid solution in previous interval
				j++;
				for (; interval_end == 1 && j < this->point_count - 1; j++) {
					interval_end = DataStructures::solve_euclidean(polyline, j, j + 1, i, this->epsilon).second;
				}
				j--;
				if (vertex_start < i && i <= j + 1 && i != 0 && i != this->point_count - 1) {
					this->solution_intervals[i].intervals.emplace_back(vertex_start, i - 1, interval_start, 1);
					this->solution_intervals[i].intervals.emplace_back(i, j, 0, interval_end);
				} else {
					this->solution_intervals[i].intervals.emplace_back(vertex_start, j, interval_start, interval_end);
				}
			}
		}
	}

	// inline void initialize_shortcut_mappings(Polyline const &polyline) {
	// 	auto index = 0u;
	// 	Queue<QueueData> queue(this->point_count);
	// 	for(PointIndex shortcut_start = 0u; shortcut_start < this->point_count - 1; shortcut_start++) {
	// 		for(PointIndex shortcut_end = shortcut_start + 1; shortcut_end < this->point_count; shortcut_end++) {
	// 			queue.reset();
	// 			this->cell_reachability(polyline, queue, shortcut_start, shortcut_end, index);
	// 			index++;
	// 		}
	// 	}
	// }
	

	inline std::generator<RangeInterval> proceeding_interval_mappings(Queue<QueueData> &queue, PointIndex shortcut_start, PointIndex shortcut_end) const {
		auto const &start_intervals = this->solution_intervals[shortcut_start].intervals;
		auto const start_interval_number = start_intervals.size();
		auto start_interval_index = 0u;

		auto &end_intervals = this->solution_intervals[shortcut_end].intervals;
		auto const end_interval_number = end_intervals.size();
		auto end_interval_index = 0u;
		auto end_interval_index_end = 0u;

		while(start_interval_index < start_interval_number) {
			// skip all end intervals that are before *all* starting intervals and thus cannot be the end of admissible subpolylines
			for (; end_interval_index < end_interval_number && 
				end_intervals[end_interval_index].end_vertex + end_intervals[end_interval_index].rel_interval_end
				< start_intervals[start_interval_index].start_vertex + start_intervals[start_interval_index].rel_interval_start
				; end_interval_index++);

			// no end intervals in which a subpolyline can end
			if (end_interval_index >= end_interval_number) {
				co_return;
			}

			queue.push_front({
				.first_interval_index = start_interval_index,
				.last_interval_index = start_interval_index,
				.first_reachable = 0.0f,
			});

			PointIndex current_vertex = start_intervals[start_interval_index].end_vertex + 1;
			start_interval_index++;
			while (!queue.is_empty()) {
				ReachabilityData const solution = current_vertex == this->point_count? DataStructures::EMPTY_INTERVAL_EXPLICIT : DataStructures::solve_euclidean(polyline, shortcut_start, shortcut_end, current_vertex, this->epsilon);

				// possibly merge starting intervals
				unsigned int merged_interval_start_index = this->point_count + 1;
				unsigned int merged_interval_end_index = queue.peek_front().last_interval_index;
				while (!queue.is_empty() && queue.peek_front().first_reachable < solution.first) {
					merged_interval_start_index = queue.peek_front().first_interval_index;
					queue.pop_front();
				}

				// new interval may start in current line segment so needs to be considered for merging 
				if(start_interval_index < start_interval_number && start_intervals[start_interval_index].start_vertex == current_vertex - 1) {
					merged_interval_start_index = std::min(start_interval_index, merged_interval_start_index);
					merged_interval_end_index = start_interval_index; 
					start_interval_index++;
				}

				// merge intervals if there actually is something to merge 
				if(merged_interval_end_index != this->point_count + 1) {
					queue.push_front({
						.first_interval_index = merged_interval_start_index,
						.last_interval_index  = merged_interval_end_index,
						// solution interval may be empty (value -1). To not cause problems later when removing, set value to 0 in this case
						.first_reachable = std::max(solution.first, 0.0f), 
					});
				}


				while (!queue.is_empty() && queue.peek_back().first_reachable > solution.second) {
					auto const [first_index, last_index, _] = queue.peek_back();
					queue.pop_back();

					for (auto i = first_index; i <= last_index; i++) {
						for (; end_interval_index < end_interval_number && 
							end_intervals[end_interval_index].end_vertex + end_intervals[end_interval_index].rel_interval_end
							< start_intervals[i].start_vertex + start_intervals[i].rel_interval_start
							; end_interval_index++) ;
						if (end_interval_index >= end_interval_number) {
							// no more intervals to map to
							co_return;
						} else if(end_intervals[end_interval_index].start_vertex >= current_vertex) {
							continue;
						}

						for (; end_interval_index_end < end_interval_number && 
							end_intervals[end_interval_index_end].end_vertex + end_intervals[end_interval_index_end].rel_interval_end <= current_vertex ; end_interval_index_end++) ;
						end_interval_index_end--;

						// update solution interval of end to be smaller to avoid ugly mappings
						if (end_intervals[end_interval_index].start_vertex + end_intervals[end_interval_index].rel_interval_start < start_intervals[i].start_vertex + start_intervals[i].rel_interval_start) {
							end_intervals[end_interval_index].start_vertex = start_intervals[i].start_vertex;
							end_intervals[end_interval_index].rel_interval_start = start_intervals[i].rel_interval_start;
						}

						co_yield RangeInterval(i, end_interval_index, end_interval_index_end);
					}
				}
				current_vertex++;
			}
		}

	}


	#ifdef test
	inline void cell_reachability(Polyline const &polyline, Queue<QueueData> &queue, PointIndex shortcut_start, PointIndex shortcut_end, unsigned int index) {
		auto const &start_intervals = this->solution_intervals[shortcut_start].intervals;
		auto const start_interval_number = start_intervals.size();
		auto start_interval_index = 0u;

		auto const &end_intervals = this->solution_intervals[shortcut_end].intervals;
		auto const end_interval_number = end_intervals.size();
		auto end_interval_index = 0u;
		auto end_interval_index_end = 0u;

		while(start_interval_index < start_interval_number) {
			// skip all end intervals that are before *all* starting intervals and thus cannot be the end of admissible subpolylines
			for (; end_interval_index < end_interval_number && 
				end_intervals[end_interval_index].end_vertex + end_intervals[end_interval_index].rel_interval_end
				< start_intervals[start_interval_index].start_vertex + start_intervals[start_interval_index].rel_interval_start
				; end_interval_index++);

			// no end intervals in which a subpolyline can end
			if (end_interval_index >= end_interval_number) {
				return;
			}

			queue.push_front({
				.first_interval_index = start_interval_index,
				.last_interval_index = start_interval_index,
				.first_reachable = 0.0f,
			});

			PointIndex current_vertex = start_intervals[start_interval_index].end_vertex + 1;
			start_interval_index++;
			while (!queue.is_empty()) {
				ReachabilityData const solution = current_vertex == this->point_count? DataStructures::EMPTY_INTERVAL_EXPLICIT : DataStructures::solve_euclidean(polyline, shortcut_start, shortcut_end, current_vertex, this->epsilon);

				// possibly merge starting intervals
				unsigned int merged_interval_start_index = this->point_count + 1;
				unsigned int merged_interval_end_index = queue.peek_front().last_interval_index;
				while (!queue.is_empty() && queue.peek_front().first_reachable < solution.first) {
					merged_interval_start_index = queue.peek_front().first_interval_index;
					queue.pop_front();
				}

				// new interval may start in current line segment so needs to be considered for merging 
				if(start_interval_index < start_interval_number && start_intervals[start_interval_index].start_vertex == current_vertex - 1) {
					merged_interval_start_index = std::min(start_interval_index, merged_interval_start_index);
					merged_interval_end_index = start_interval_index; 
					start_interval_index++;
				}

				// merge intervals if there actually is something to merge 
				if(merged_interval_end_index != this->point_count + 1) {
					queue.push_front({
						.first_interval_index = merged_interval_start_index,
						.last_interval_index  = merged_interval_end_index,
						// solution interval may be empty (value -1). To not cause problems later when removing, set value to 0 in this case
						.first_reachable = std::max(solution.first, 0.0f), 
					});
				}


				while (!queue.is_empty() && queue.peek_back().first_reachable > solution.second) {
					auto const [first_index, last_index, _] = queue.peek_back();
					queue.pop_back();

					for (auto i = first_index; i <= last_index; i++) {
						for (; end_interval_index < end_interval_number && 
							end_intervals[end_interval_index].end_vertex + end_intervals[end_interval_index].rel_interval_end
							< start_intervals[i].start_vertex + start_intervals[i].rel_interval_start
							; end_interval_index++) ;
						if (end_interval_index >= end_interval_number) {
							// no more intervals to map to
							return;
						} else if(end_intervals[end_interval_index].start_vertex >= current_vertex) {
							continue;
						}

						for (; end_interval_index_end < end_interval_number && 
							end_intervals[end_interval_index_end].end_vertex + end_intervals[end_interval_index_end].rel_interval_end <= current_vertex ; end_interval_index_end++) ;
						end_interval_index_end--;

						// otherwise use the end interval 
						bool const use_start_interval_as_bound_for_start = start_intervals[i].start_vertex + start_intervals[i].rel_interval_start >= end_intervals[end_interval_index].start_vertex + end_intervals[end_interval_index].rel_interval_start;

						this->global_shortcut_mappings[index].shortcut_mappings.emplace_back(start_intervals[i], Interval{
							// overestimation of actual interval, but with necessary intersections correct
							.start_vertex = use_start_interval_as_bound_for_start? start_intervals[i].start_vertex: end_intervals[end_interval_index].start_vertex,
							.end_vertex = end_intervals[end_interval_index_end].end_vertex,
							.rel_interval_start = use_start_interval_as_bound_for_start? start_intervals[i].rel_interval_start: end_intervals[end_interval_index].rel_interval_start,
							.rel_interval_end = end_intervals[end_interval_index_end].rel_interval_end,
						});
					}
				}
				current_vertex++;
			}
		}
	}

#endif


  void print() {
		std::cout << "Printing the computed intervals using epsilon = " << this->epsilon << ":" << std::endl;
		for (auto i = 0u; i < this->point_count; i++) {
			std::cout << "Point " << i << ": ";
			for (auto const &interval : this->solution_intervals[i].intervals) {
				std::cout << "(" << interval.start_vertex << " + " << interval.rel_interval_start << ", " << interval.end_vertex << " + " << interval.rel_interval_end << "), ";
			}
			std::cout << std::endl;
		}

		Queue<QueueData> queue(this->point_count);
		for(PointIndex shortcut_start = 0u; shortcut_start < this->point_count - 1; shortcut_start++) {
			for(PointIndex shortcut_end = shortcut_start + 1; shortcut_end < this->point_count; shortcut_end++) {
				std::cout << "Line segment " << " (" << shortcut_start << ", " << shortcut_end << ")" << ": ";
				for (auto const &[dom, interval_start, interval_end] : this->proceeding_interval_mappings(queue, shortcut_start, shortcut_end) ) {
					std::cout << "(" << dom << "[" << interval_start << ", "<< interval_end <<  "]); ";
					//std::cout << "[(" << interval_start.start_vertex << " + " << interval_start.rel_interval_start << ", " << interval_start.end_vertex << " + " << interval_start.rel_interval_end << "), (" << interval_end.start_vertex << " + " << interval_end.rel_interval_start << ", " << interval_end.end_vertex << " + " << interval_end.rel_interval_end << ")]; ";
				}
				queue.reset();
				std::cout << std::endl;
			}
		}


		std::cout << "Again: Printing the computed intervals using epsilon = " << this->epsilon << ":" << std::endl;
		for (auto i = 0u; i < this->point_count; i++) {
			std::cout << "Point " << i << ": ";
			for (auto const &interval : this->solution_intervals[i].intervals) {
				std::cout << "(" << interval.start_vertex << " + " << interval.rel_interval_start << ", " << interval.end_vertex << " + " << interval.rel_interval_end << "), ";
			}
			std::cout << std::endl;
		}
	}
};


static inline bool is_local_shortcut(GlobalShortcutGraph const &graph, PointIndex start, PointIndex end) {
	// only need to find mapping that contains start in the start->end mappings and check if end is in mapped range. 
	// if number of intervals is small use linear search, else binary search (can also be done in constant time with more precomputations but those are not needed in global case and require cubic time)
	
	Queue<QueueData> queue(graph.point_count);
	for (auto const &[i, i_start, i_end] : graph.proceeding_interval_mappings(queue, start, end)) {
		auto const &domain      = graph.solution_intervals[start].intervals[i];
		auto const &range_start = graph.solution_intervals[end].intervals[i_start];
		auto const &range_end   = graph.solution_intervals[end].intervals[i_end];
		if (domain.end_vertex + domain.rel_interval_end >= start) {
			return domain.start_vertex + domain.rel_interval_start <= start && (range_start.start_vertex + range_start.rel_interval_start <= end && range_end.end_vertex + range_end.rel_interval_end >= end);
		}
	}

	return false;
}

Simplification simplification_local_imai_iri_from_gsg(GlobalShortcutGraph const &graph) {
	unsigned int *sizes = new unsigned int[graph.point_count];
	unsigned int *path  = new unsigned int[graph.point_count];
	sizes[0] = 0;
	for (unsigned int i = 1; i < graph.point_count; i++) {
		sizes[i] = sizes[i-1] + 1;
		path[i] = i - 1;
		for (unsigned int j = 0; j < i; j++) {
			if (sizes[j] + 1 < sizes[i] && is_local_shortcut(graph, j, i)) {
				sizes[i] = sizes[j] + 1;
				path[i] = j;
			}
		}
	}

	unsigned int k = graph.point_count - 1;
	unsigned int index = sizes[k];
	Simplification result = std::make_unique<std::vector<size_t>>(index + 1);
	while (index > 0) {
		(*result)[index] = k;
		k = path[k];
		index--;
	}
	
	(*result)[0] = 0;
	return result;
}



Simplification simplification_global_imai_iri_euclidean(Polyline const &polyline, float epsilon, AlgorithmConfiguration &config) {
	auto time_start = std::chrono::high_resolution_clock::now();

	GlobalShortcutGraph graph(polyline, epsilon);
	graph.print();

	auto local_simpl = simplification_local_imai_iri_from_gsg(graph);

	if (config.logger.has_value()) {
		auto end = std::chrono::high_resolution_clock::now();
		if (config.logger.has_value()) {
			// config.logger.value().add_data(simplification->size(), end - time_start, "");
			config.logger.value().add_data(0, end - time_start, "");
		}
	}

	return local_simpl;
}
}
