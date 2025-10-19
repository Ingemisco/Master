#if __has_include(<generator>)
#include <cassert>
#include <cstddef>
#include <generator>
#include <iostream>
#include <limits>
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

struct SolutionIntervals final {
	std::vector<Interval> intervals;
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


	GlobalShortcutGraph(Polyline const &polyline, float epsilon) : polyline(polyline), point_count(polyline.point_count), epsilon(epsilon) {
    this->solution_intervals = new SolutionIntervals[polyline.point_count];
		this->initialize_solution_intervals();
	}

	~GlobalShortcutGraph() {
		delete[] this->solution_intervals;
	}

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

			// just in case of numerical problems, add interval [i, i] for vertex i if no intervals have been found 
			if (this->solution_intervals[i].intervals.empty()) {
				this->solution_intervals[i].intervals.emplace_back(i, i, 0, 0);
			}
		}

		// None of this should be necessary but somehow it is???
		if (this->solution_intervals[0].intervals[0].start_vertex + this->solution_intervals[0].intervals[0].rel_interval_start > 0) {
			this->solution_intervals[0].intervals.insert(this->solution_intervals[0].intervals.begin(), Interval{0,0,0,0});
		}
		if (this->solution_intervals[point_count-1].intervals.back().end_vertex + this->solution_intervals[point_count-1].intervals.back().rel_interval_end < point_count - 1) {
			this->solution_intervals[point_count-1].intervals.emplace_back(point_count - 2, point_count - 2, 1, 1);
		}

	}

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

			PointIndex current_vertex = std::min(start_intervals[start_interval_index].end_vertex + 1, this->point_count - 1);
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
				if(merged_interval_start_index != this->point_count + 1) {
					queue.push_front({
						.first_interval_index = merged_interval_start_index,
						.last_interval_index  = merged_interval_end_index,
						// solution interval may be empty (value -1). To not cause problems later when removing, set value to 0 in this case
						.first_reachable = std::max(solution.first, 0.0f), 
					});
				}

				while (!queue.is_empty() && (queue.peek_back().first_reachable > solution.second || current_vertex >= this->point_count) ) {
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
						} else if(end_intervals[end_interval_index].start_vertex + end_intervals[end_interval_index].rel_interval_start > current_vertex) {
							continue;
						}

						for (; end_interval_index_end < end_interval_number && 
							end_intervals[end_interval_index_end].start_vertex + end_intervals[end_interval_index_end].rel_interval_start < current_vertex ; end_interval_index_end++) ;

							// I have absolutely no idea how this value can be zero, but it can.
						if (end_interval_index_end != 0) {
							end_interval_index_end--;
						}

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

  void print() {
		std::cout << "Printing the computed intervals using epsilon = " << this->epsilon << ":" << std::endl;
		for (auto i = 0u; i < this->point_count; i++) {
			std::cout << "Point " << i << ": ";
			for (auto const &interval : this->solution_intervals[i].intervals) {
				std::cout << "(" << interval.start_vertex << " + " << interval.rel_interval_start << ", " << interval.end_vertex << " + " << interval.rel_interval_end << "), ";
			}
			std::cout << std::endl;
		}

		Queue<QueueData> queue(this->point_count + 5);
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
	
	Queue<QueueData> queue(graph.point_count + 5);
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
		for (unsigned int j = 0; j < i - 1; j++) {
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

	struct SimplifificationData final {
		size_t size;
		size_t parent_vertex;
		size_t parent_interval;
	};


void printSimplificationDebugInfo(SimplifificationData* simplification_data, 
                                 size_t* offset_array, 
                                 size_t point_count,
                                 GlobalShortcutGraph& graph) {
    std::cout << "\n=== SIMPLIFICATION DEBUG INFO ===" << std::endl;
    
    // Print offset array
    std::cout << "\nOFFSET ARRAY (size: " << point_count + 1 << "):" << std::endl;
    std::cout << "Index | Offset" << std::endl;
    std::cout << "------|-------" << std::endl;
    for (size_t i = 0; i <= point_count; i++) {
        std::cout << i << "     | " << offset_array[i] << std::endl;
    }
    
    std::cout << "\nTotal simplification_data entries: " << offset_array[point_count] << std::endl;
    
    // Print simplification data for each vertex
    for (size_t vertex = 0; vertex < point_count; vertex++) {
        size_t start_idx = offset_array[vertex];
        size_t end_idx = offset_array[vertex + 1];
        size_t interval_count = end_idx - start_idx;
        
        std::cout << "\n--- Vertex " << vertex << " ---" << std::endl;
        std::cout << "Intervals: " << interval_count << " (indices " << start_idx << " to " << (end_idx - 1) << ")" << std::endl;
        
        if (interval_count == 0) {
            std::cout << "  No intervals for this vertex" << std::endl;
            continue;
        }
        
        // Print actual intervals from graph
        auto& intervals = graph.solution_intervals[vertex].intervals;
        std::cout << "Actual intervals in graph:" << std::endl;
        for (size_t interval_idx = 0; interval_idx < intervals.size(); interval_idx++) {
            auto& interval = intervals[interval_idx];
            std::cout << "  [" << interval_idx << "] (" << interval.start_vertex 
                      << " + " << interval.rel_interval_start << ", " 
                      << interval.end_vertex << " + " << interval.rel_interval_end << ")" << std::endl;
        }
        
        // Print simplification data for each interval of this vertex
        std::cout << "Simplification data:" << std::endl;
        for (size_t local_idx = 0; local_idx < interval_count; local_idx++) {
            size_t global_idx = start_idx + local_idx;
            auto& data = simplification_data[global_idx];
            
            std::cout << "  Interval " << local_idx << " (global " << global_idx << "): ";
            
            if (data.size == std::numeric_limits<size_t>::max()) {
                std::cout << "UNREACHABLE";
            } else {
                std::cout << "size=" << data.size 
                          << ", parent_vertex=" << data.parent_vertex 
                          << ", parent_interval=" << data.parent_interval;
                
                // Validate parent indices
                if (data.parent_vertex < point_count) {
                    size_t parent_start = offset_array[data.parent_vertex];
                    size_t parent_end = offset_array[data.parent_vertex + 1];
                    if (data.parent_interval < (parent_end - parent_start)) {
                        std::cout << " [VALID PARENT]";
                    } else {
                        std::cout << " [INVALID parent_interval: " << data.parent_interval 
                                  << " >= " << (parent_end - parent_start) << "]";
                    }
                } else {
                    std::cout << " [INVALID parent_vertex: " << data.parent_vertex 
                              << " >= " << point_count << "]";
                }
            }
            std::cout << std::endl;
        }
    }
    
    // Print final path reconstruction info
    std::cout << "\n=== FINAL PATH RECONSTRUCTION ===" << std::endl;
    size_t final_vertex = point_count - 1;
    size_t final_offset = offset_array[final_vertex];
    size_t final_interval_count = offset_array[final_vertex + 1] - final_offset;
    
    std::cout << "Final vertex: " << final_vertex << std::endl;
    std::cout << "Final interval count: " << final_interval_count << std::endl;
    
    if (final_interval_count > 0) {
        size_t last_global_idx = offset_array[point_count] - 1;
        auto& final_data = simplification_data[last_global_idx];
        
        std::cout << "Using simplification_data[" << last_global_idx << "]:" << std::endl;
        std::cout << "  Size: " << final_data.size << std::endl;
        std::cout << "  Parent vertex: " << final_data.parent_vertex << std::endl;
        std::cout << "  Parent interval: " << final_data.parent_interval << std::endl;
        
        if (final_data.size == std::numeric_limits<size_t>::max()) {
            std::cout << "  WARNING: Final vertex is UNREACHABLE!" << std::endl;
        }
    } else {
        std::cout << "ERROR: Final vertex has no intervals!" << std::endl;
    }
    
    std::cout << "==================================\n" << std::endl;
}












Simplification simplification_global_imai_iri_euclidean(Polyline const &polyline, float epsilon, AlgorithmConfiguration &config) {
	auto time_start = std::chrono::high_resolution_clock::now();
	GlobalShortcutGraph graph(polyline, epsilon);
	//graph.print();

	//auto local_simpl = simplification_local_imai_iri_from_gsg(graph);

	// offset_array[i] is offset in simplification_data array at which data for the vertex i can be found. 
	size_t *offset_array = new size_t[polyline.point_count + 1];
	offset_array[0] = 0;
	for (auto i = 0u; i < polyline.point_count; i++) {
		offset_array[i + 1] = offset_array[i] + graph.solution_intervals[i].intervals.size(); 
	}

	SimplifificationData *simplification_data = new SimplifificationData[offset_array[polyline.point_count]];
	for (auto i = 1u; i < offset_array[polyline.point_count]; i++) {
		simplification_data[i].size = std::numeric_limits<size_t>::max();
	}
	simplification_data[0].size = 1; // parentdata for start is never used
	Queue<QueueData> queue(graph.point_count + 5);
	
	unsigned int last_interval_that_contains_vertex = 0u;
	for (auto i = 1u; i < polyline.point_count; i++) {
		for (auto j = 0u; j < i; j++) {
			for (auto const &[dom, i_start, i_end] : graph.proceeding_interval_mappings(queue, j, i)) {
				for (auto k = i_start; k <= i_end; k++) {
					if (simplification_data[offset_array[j] + dom].size < std::numeric_limits<size_t>::max() && 
						simplification_data[offset_array[i] + k].size > simplification_data[offset_array[j] + dom].size + 1) {
						simplification_data[offset_array[i] + k].size = simplification_data[offset_array[j] + dom].size + 1;
						simplification_data[offset_array[i] + k].parent_vertex = j;
						simplification_data[offset_array[i] + k].parent_interval = dom;
					}
				}
			}
			queue.reset();
		}


		unsigned int interval_that_contains_vertex; // completely invalid. Interval must be found
		for (unsigned int interval_index = 0u; interval_index < graph.solution_intervals[i].intervals.size(); interval_index++) {
			if (graph.solution_intervals[i].intervals[interval_index].end_vertex + graph.solution_intervals[i].intervals[interval_index].rel_interval_end >= i) {
				interval_that_contains_vertex = interval_index;
				goto proceed;
			}
		}
		// should not be reachable
		std::cout << "ERROR did not find interval for index i = " << i << std::endl; graph.print();
		exit(0);
	proceed:
		if (simplification_data[offset_array[i] + interval_that_contains_vertex].size == std::numeric_limits<size_t>::max()) {
			simplification_data[offset_array[i] + interval_that_contains_vertex].size = simplification_data[offset_array[i-1] + last_interval_that_contains_vertex].size + 1;
			simplification_data[offset_array[i] + interval_that_contains_vertex].parent_vertex = i - 1;
			simplification_data[offset_array[i] + interval_that_contains_vertex].parent_interval = last_interval_that_contains_vertex;
		}
		last_interval_that_contains_vertex = interval_that_contains_vertex;
	}

	//printSimplificationDebugInfo(simplification_data, offset_array, polyline.point_count, graph);

	auto vertex = simplification_data[offset_array[polyline.point_count] - 1].parent_vertex;
	auto interval = simplification_data[offset_array[polyline.point_count] - 1].parent_interval;
	auto const size = simplification_data[offset_array[polyline.point_count] - 1].size;
	auto i = size - 1;

	Simplification result = std::make_unique<std::vector<size_t>>(size);
	(*result)[i] = polyline.point_count - 1;
	i--;
	while (i > 0) {
		(*result)[i] = vertex;
		auto const new_vertex = simplification_data[offset_array[vertex] + interval].parent_vertex;
		interval = simplification_data[offset_array[vertex] + interval].parent_interval;
		vertex = new_vertex;
		i--;
	}
	
	(*result)[0] = 0;

	if (config.logger.has_value()) {
		auto end = std::chrono::high_resolution_clock::now();
		if (config.logger.has_value()) {
			config.logger.value().add_data(result->size(), end - time_start, "");
		}
	}

	delete[] offset_array;
	delete[] simplification_data;

	return result;
}
}

#endif
