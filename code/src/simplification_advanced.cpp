#include "datastructures.h"
#include "distance.h"
#include "simplification.h"
#include "log.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <memory>
#include <ostream>
#include <tuple>
#include <utility>
#include <vector>


using DataStructures::Polyline, DataStructures::Queue;

struct CRDataSemiExplicit final {
	size_t k;
	float a;
	float d;
};

// Important: intervals must have size at least as large as entry_costs, even though one entry less is required.
// The last entry must not contain -1 as first component, Just set it to [0, 1] for simplicity
//
// first component of exit_costs is value, second is the index of entry_costs it came from (only if the exit_cost is not INDEX_UNREACHABLE)
template <typename F, typename L, std::pair<F, L> EMPTY_INTERVAL>
static void cell_reachability_explicit(size_t const *entry_costs, std::pair<F, L> const *intervals,
																			 std::pair<size_t, size_t> *exit_costs, Queue<std::tuple<F, size_t, size_t>> &queue) {
	size_t const n = queue.capacity;
	auto last_interval = EMPTY_INTERVAL;
	for (size_t j = 0; j < n; j++) {
		if (last_interval == EMPTY_INTERVAL) {
			queue.reset();
			do {
				exit_costs[j].first = DataStructures::INDEX_UNREACHABLE;
				j++;
			} while (intervals[j-1] == EMPTY_INTERVAL);
			if (j == n) return;
			last_interval = intervals[j-1];
		}

		size_t const lambda_j = entry_costs[j - 1];
		size_t k_left = lambda_j;
		size_t ref = j - 1;
		while(!queue.is_empty()) {
			auto const [t, k, ref_index] = queue.peek_front();
			if (k < lambda_j && last_interval.first < t) break;
			if (k < k_left) {
				k_left = k;
				ref = ref_index;
			}
			queue.pop_front();
		}

		queue.push_front(std::tuple<F, size_t, size_t>(last_interval.first, k_left, ref));

		while(!queue.is_empty()) {
			F const t = std::get<0>(queue.peek_back());
			if(t <= last_interval.second) break;
			queue.pop_back();
		}

		auto const [start, val, i] = queue.peek_back(); 
		exit_costs[j].first = val;
		exit_costs[j].second = i;
		last_interval = intervals[j];
	}
}



template <typename T>
struct _DP1Data final {
	T first_reachable_point;
	size_t i_;
};

struct KappaData final {
	size_t smallest_k;
	size_t i_;
	size_t j_;
};






void printKappaData(KappaData* kappa, KappaData* kappa2, size_t n) {
    std::cout << "\n=== KAPPA DATA DUMP ===" << std::endl;
    std::cout << "Array dimensions: " << n << " x " << (n-1) << " = " << (n * (n-1)) << " total elements" << std::endl;
    
    size_t total_kappa_unreachable = 0;
    size_t total_kappa2_unreachable = 0;
    
    // Print kappa array
    std::cout << "\n--- KAPPA ARRAY ---" << std::endl;
    for (size_t i = 0; i < n; i++) {
        bool has_valid_data = false;
        std::vector<size_t> valid_indices;
        
        for (size_t j = 0; j < n-1; j++) {
            size_t index = i * (n-1) + j;
            if (kappa[index].smallest_k != DataStructures::INDEX_UNREACHABLE) {
                has_valid_data = true;
                valid_indices.push_back(j);
            }
        }
        
        if (has_valid_data) {
            std::cout << "Vertex " << i << " - Valid intervals: ";
            for (size_t j : valid_indices) {
                size_t index = i * (n-1) + j;
                std::cout << "j=" << j << "(k=" << kappa[index].smallest_k 
                          << ", i_=" << kappa[index].i_ 
                          << ", j_=" << kappa[index].j_ << ") ";
            }
            std::cout << std::endl;
        } else {
            total_kappa_unreachable++;
        }
    }
    
    // Print kappa2 array  
    std::cout << "\n--- KAPPA2 ARRAY ---" << std::endl;
    for (size_t i = 0; i < n; i++) {
        bool has_valid_data = false;
        std::vector<size_t> valid_indices;
        
        for (size_t j = 0; j < n-1; j++) {
            size_t index = i * (n-1) + j;
            if (kappa2[index].smallest_k != DataStructures::INDEX_UNREACHABLE) {
                has_valid_data = true;
                valid_indices.push_back(j);
            }
        }
        
        if (has_valid_data) {
            std::cout << "Vertex " << i << " - Valid intervals: ";
            for (size_t j : valid_indices) {
                size_t index = i * (n-1) + j;
                std::cout << "j=" << j << "(k=" << kappa2[index].smallest_k 
                          << ", i_=" << kappa2[index].i_ 
                          << ", j_=" << kappa2[index].j_ << ") ";
            }
            std::cout << std::endl;
        } else {
            total_kappa2_unreachable++;
        }
    }
    
    // Print summary statistics
    std::cout << "\n=== SUMMARY ===" << std::endl;
    std::cout << "Kappa: " << total_kappa_unreachable << "/" << n << " vertices completely unreachable" << std::endl;
    std::cout << "Kappa2: " << total_kappa2_unreachable << "/" << n << " vertices completely unreachable" << std::endl;
    
    // Print the critical last elements
    std::cout << "\n=== CRITICAL INDICES ===" << std::endl;
    size_t last_index_kappa = (n-1) * (n-1);
    size_t last_index_kappa_alt = (n-1) * n - 1;
    
    std::cout << "Last index used: " << last_index_kappa << std::endl;
    std::cout << "Alternative last index: " << last_index_kappa_alt << std::endl;
    std::cout << "Array size: " << (n * (n-1)) << std::endl;
    
    if (last_index_kappa < n * (n-1)) {
        std::cout << "kappa[" << last_index_kappa << "]: smallest_k=" << kappa[last_index_kappa].smallest_k
                  << ", i_=" << kappa[last_index_kappa].i_ 
                  << ", j_=" << kappa[last_index_kappa].j_ << std::endl;
    } else {
        std::cout << "ERROR: last_index_kappa out of bounds!" << std::endl;
    }
    
    if (last_index_kappa_alt < n * (n-1)) {
        std::cout << "kappa[" << last_index_kappa_alt << "]: smallest_k=" << kappa[last_index_kappa_alt].smallest_k
                  << ", i_=" << kappa[last_index_kappa_alt].i_ 
                  << ", j_=" << kappa[last_index_kappa_alt].j_ << std::endl;
    } else {
        std::cout << "ERROR: last_index_kappa_alt out of bounds!" << std::endl;
    }
    
    // Find the actual maximum k value
    size_t max_k = 0;
    size_t max_k_index = 0;
    for (size_t i = 0; i < n * (n-1); i++) {
        if (kappa[i].smallest_k != DataStructures::INDEX_UNREACHABLE && 
            kappa[i].smallest_k > max_k) {
            max_k = kappa[i].smallest_k;
            max_k_index = i;
        }
    }
    std::cout << "Maximum k value: " << max_k << " at index " << max_k_index << std::endl;
}












// EMPTY_INTERVAL must satisfy that its first component is the value to be used for invalid entries, so it should not be a valid entry
// NON_EMPTY_INTERVAL must be a non empty interval whose first point is something equivalent to 0 (i.e. start of the line segment)
template <typename F, typename L, std::pair<F const, L const> eq_solver(Polyline const &, size_t, size_t, size_t, float),
	std::pair<F, L> NON_EMPTY_INTERVAL, std::pair<F, L> EMPTY_INTERVAL>
struct Simplifier final {
	
	using Interval  = std::pair<F, L>;
	using CellQueue = Queue<std::tuple<F, size_t, size_t>>;
	using ExitCost  = std::pair<size_t, size_t>;


	Polyline const &polyline;
	_DP1Data<F> *dp1;
	KappaData *kappa;
	KappaData *kappa2;


	Interval *reachability_data;
	Interval *intervals;
	size_t *entry_costs;
	ExitCost *exit_costs;
	CellQueue queue;
	float const epsilon;
	size_t const n;


	// NOTE: first coordinate is now i, second j, third k
	Simplifier(Polyline const &polyline, float epsilon) 
		: polyline(polyline), queue(polyline.point_count), epsilon(epsilon), n(polyline.point_count) {
		
		this->reachability_data = new Interval[n * (n - 1)];
		this->entry_costs   = new size_t[n];
		this->exit_costs    = new ExitCost[n];
		this->intervals     = new Interval[n];
		this->dp1           = new _DP1Data<F>[n * (n - 1) * n];
		this->kappa         = new KappaData[n * (n - 1)];
		this->kappa2        = new KappaData[n * (n - 1)];
	}

	~Simplifier() {
		delete[] this->reachability_data;
		delete[] this->entry_costs;
		delete[] this->exit_costs;
		delete[] this->intervals;
		delete[] this->dp1;
		delete[] this->kappa2;
		delete[] this->kappa;
	}

	inline void initialize() {
		this->intervals[n - 1] = NON_EMPTY_INTERVAL;

		unsigned int index = 0;
		for (unsigned int i = 0; i < n; i++) {
			for (unsigned int j = 0; j < n - 1; j++) {
				this->reachability_data[index++] = eq_solver(polyline, j, j+1, i, epsilon);
			}
		}

		size_t j = 1;
		for(; j < n - 1 && reachability_data[j].first == NON_EMPTY_INTERVAL.first; j++);
		for (size_t j_ = 0; j_ < j; j_++) {
			for (size_t k = 0; k < n; k++) {
				dp1[k + n * j_].first_reachable_point = NON_EMPTY_INTERVAL.first;
				dp1[k + n * j_].i_ = 0;
			}
		}
		for (size_t j_ = j; j_ < n - 1; j_++) {
			for (size_t k = 0; k < n; k++) {
				dp1[k + n * j_].first_reachable_point = EMPTY_INTERVAL.first;
			}
		}
		for (size_t i = 1; i < n; i++) {
			for (size_t j_ = 0; j_ < n - 1; j_++) {
				dp1[n * (j_ + (n - 1) * i) ].first_reachable_point = EMPTY_INTERVAL.first;
			}
		}
		for (size_t k = 0; k < j; k++) {
			kappa[k].smallest_k = 0;
			kappa[k].i_ = 0;
			kappa[k].j_ = 0;
		}
		for (size_t k = j; k < n-1; k++) {
			kappa[k].smallest_k = DataStructures::INDEX_UNREACHABLE;
		}
		for (size_t t = 0; t < n * (n-1); t++) {
			kappa2[t].smallest_k = DataStructures::INDEX_UNREACHABLE;
		}
	}

	inline F dp(size_t k, size_t i, size_t j) {
		if (k == 0 && 0 < i) {
			return EMPTY_INTERVAL.first;
		}
		size_t const pos = (n - 1) * i + j;
		auto const range = reachability_data[pos];
		if (kappa2[pos].smallest_k <= k) {
			return range.first;
		} 
		F const t = dp1[pos * n +  k].first_reachable_point;
		if (t != EMPTY_INTERVAL.first && t <= range.second) {
			return std::max(t, range.first);
		}
		return EMPTY_INTERVAL.first;
	}

	inline void kappa2_subroutine(size_t i) {
		for (size_t i_ = 0; i_ < i; i_++) {
			for (size_t j = 0; j < n - 1; j++) {
				this->intervals[j] = eq_solver(polyline, i_, i, j+1, epsilon);

				auto const k = kappa[j + (n-1) * i_].smallest_k;
				this->entry_costs[j] = DataStructures::INDEX_UNREACHABLE;
				if (k != DataStructures::INDEX_UNREACHABLE) {
					this->entry_costs[j] = k + 1;
				}
			}

			cell_reachability_explicit<F, L, EMPTY_INTERVAL>(const_cast<size_t const *>(entry_costs), const_cast<Interval const *>(intervals), exit_costs, queue);
			

			size_t temp2 = i * (n - 1);
			for (size_t j = 0; j < n - 1; j++) {
				auto const exit = this->exit_costs[j];
				if (reachability_data[j + (n-1)*i].first == EMPTY_INTERVAL.first) {
					this->kappa2[temp2].smallest_k = DataStructures::INDEX_UNREACHABLE;
				} else if (exit.first < this->kappa2[temp2].smallest_k) {
					this->kappa2[temp2].smallest_k = exit.first;
					this->kappa2[temp2].i_ = i_;
					this->kappa2[temp2].j_ = exit.second;
				}
				temp2++;
			}
		}

		// from i-1 to i is always valid with j' = i-1, j = i and t' = t = 0. This is sometimes not found for small epsilon so artificially add it if necessary (also ensures that a simplification is always found)
		// In a just world this would not be necessary...
		if (i < n - 1 && this->kappa2[i * (n-1) + i].smallest_k - 1 > this->kappa[(i-1) * (n-1) + (i - 1)].smallest_k) {
			this->kappa2[i * (n-1) + (i)].smallest_k = this->kappa[(i-1) * (n-1) + (i - 1)].smallest_k + 1;
			this->kappa2[i * (n-1) + (i)].i_ = i - 1;
			this->kappa2[i * (n-1) + (i)].j_ = i - 1;
		}

		// ... and because the world hates me, this is also necessary...
		if (this->kappa2[i * (n-1) + i - 1].smallest_k - 1 > this->kappa[(i-1) * (n-1) + (i - 1)].smallest_k) {
			this->kappa2[i * (n-1) + (i - 1)].smallest_k = this->kappa[(i-1) * (n-1) + (i - 1)].smallest_k + 1;
			this->kappa2[i * (n-1) + (i - 1)].i_ = i - 1;
			this->kappa2[i * (n-1) + (i - 1)].j_ = i - 1;
		}
	}

  inline Simplification::Simplification simplify() {
		this->initialize();
		for (size_t i = 1; i < n; i++) {
			kappa2_subroutine(i);

			for (size_t j = 0; j < n - 1; j++) {
				for (size_t k = 1; k < n; k++) {
					size_t const pos = k + n * (j + (n - 1) * i);
					auto const v1 = dp1[k + n * (j + (n - 1) * (i-1))];
					auto const v2 = dp(k - 1, i - 1, j);
					dp1[pos].first_reachable_point = EMPTY_INTERVAL.first;
					dp1[pos].i_ = DataStructures::INDEX_UNREACHABLE;
					if (v2 == EMPTY_INTERVAL.first) {
						dp1[pos].first_reachable_point = v1.first_reachable_point;
						dp1[pos].i_ = v1.i_;
					} else if(v1.first_reachable_point == EMPTY_INTERVAL.first || v2 < v1.first_reachable_point) {
						dp1[pos].first_reachable_point = v2;
						dp1[pos].i_ = i - 1;
					} else {
						dp1[pos].first_reachable_point = v1.first_reachable_point;
						dp1[pos].i_ = v1.i_;
					}
				}

				size_t kappa1 = DataStructures::INDEX_UNREACHABLE;
				size_t kappa1_ref = DataStructures::INDEX_UNREACHABLE;
				auto const sij = reachability_data[j + (n-1) * i].second;
				for (size_t k = 1; k < n; k++) {
					auto const val = dp1[k + n * (j + (n - 1) * i)].first_reachable_point;
					if (val != EMPTY_INTERVAL.first && val <= sij) {
						kappa1 = k;
						kappa1_ref = dp1[k + n * (j + (n - 1) * i)].i_;
						break;
					}
				}

				size_t const pos = (n-1) * i + j;
				if (kappa1 < kappa2[pos].smallest_k) {
					kappa[pos].smallest_k = kappa1;
					kappa[pos].i_ = kappa1_ref;
					kappa[pos].j_ = j;
				} else {
					kappa[pos] = kappa2[pos];
				}
			}
		}


		// ignore the warning because it is guaranteed by polyline construction that the size is >= 2 (because n > 1)
		[[assume(n >= 2)]];
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
		size_t k = kappa[(n-1)*n - 1].smallest_k;
#pragma GCC diagnostic pop 

		Simplification::Simplification result = std::make_unique<std::vector<size_t>>(1 + k);
		auto &p = *result;
		size_t i_ = n - 1;
		size_t j_ = n - 2;

		// printKappaData(kappa, kappa2, n);
		
		while (0 < i_) {
			p[k] = i_;
			auto const temp = kappa2[i_ * (n-1) + j_];
			if (temp.smallest_k <= k) {
				i_ = temp.i_;
				j_ = temp.j_;
			} else {
				i_ = dp1[k + n * (j_ + (n-1) * i_)].i_;
			}
			k--;
		}


		p[0] = 0;
		return result;
	}
};



Simplification::Simplification Simplification::simplification_advanced_euclidean_explicit(DataStructures::Polyline const &polyline, float epsilon, AlgorithmConfiguration &config) {
	auto time_start = std::chrono::high_resolution_clock::now();
	Simplifier<float, float, DataStructures::solve_euclidean, std::pair<float, float>(0, 1), DataStructures::EMPTY_INTERVAL_EXPLICIT> simplifier(polyline, epsilon);

	auto simplification = simplifier.simplify();

	// assumes test case was started before calling this function
	if (config.logger.has_value()) {
		auto end = std::chrono::high_resolution_clock::now();
		if (config.logger.has_value()) {
			config.logger.value().add_data(simplification->size(), end - time_start, "");
		}
	}

	return simplification;
}



Simplification::Simplification Simplification::simplification_advanced_manhattan_explicit(DataStructures::Polyline const &polyline, float epsilon, AlgorithmConfiguration &config) {
	auto time_start = std::chrono::high_resolution_clock::now();
	Simplifier<float, float, DataStructures::solve_manhattan, std::pair<float, float>(0, 1), DataStructures::EMPTY_INTERVAL_EXPLICIT> simplifier(polyline, epsilon);
	auto simplification = simplifier.simplify();

	// assumes test case was started before calling this function
	if (config.logger.has_value()) {
		auto end = std::chrono::high_resolution_clock::now();
		if (config.logger.has_value()) {
			config.logger.value().add_data(simplification->size(), end - time_start, "");
		}
	}

	return simplification;
}

Simplification::Simplification Simplification::simplification_advanced_chebyshev_explicit(DataStructures::Polyline const &polyline, float epsilon, AlgorithmConfiguration &config) {
	auto time_start = std::chrono::high_resolution_clock::now();
	Simplifier<float, float, DataStructures::solve_chebyshev, std::pair<float, float>(0, 1), DataStructures::EMPTY_INTERVAL_EXPLICIT> simplifier(polyline, epsilon);
	auto simplification = simplifier.simplify();

	// assumes test case was started before calling this function
	if (config.logger.has_value()) {
		auto end = std::chrono::high_resolution_clock::now();
		if (config.logger.has_value()) {
			config.logger.value().add_data(simplification->size(), end - time_start, "");
		}
	}

	return simplification;
}



Simplification::Simplification Simplification::simplification_advanced_euclidean_semiexplicit(DataStructures::Polyline const &polyline, float epsilon, AlgorithmConfiguration &config) {
	auto time_start = std::chrono::high_resolution_clock::now();

	Simplifier<DataStructures::FRValue, DataStructures::LRValue, DataStructures::solve_euclidean_se, DataStructures::NONEMPTY_INTERVAL_SEMIEXPLICIT, DataStructures::EMPTY_INTERVAL_SEMIEXPLICIT> simplifier(polyline, epsilon);
	auto simplification = simplifier.simplify();

	// assumes test case was started before calling this function
	if (config.logger.has_value()) {
		auto end = std::chrono::high_resolution_clock::now();
		if (config.logger.has_value()) {
			config.logger.value().add_data(simplification->size(), end - time_start, "");
		}
	}

	return simplification;
}





















