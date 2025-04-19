#include "datastructures.h"
#include "distance.h"
#include "simplification.h"
#include <algorithm>
#include <cstring>
#include <memory>
#include <ostream>
#include <tuple>

using DataStructures::Point, DataStructures::Polyline;


template<typename T>
class Queue final {
private:
	T *array;
	size_t start; // position **before** first element, i.e., index where push_front would insert 
	size_t end; // position last element 

public:
	size_t const capacity;

	Queue(size_t capacity) : start(capacity - 1), end(capacity - 1), capacity(capacity) {
    this->array = new T[2 * capacity - 1];
	}

	~Queue() {
		delete[] this->array;
	}

	Queue(const Queue&) = delete;
	Queue(Queue &&) = delete;
	Queue& operator=(const Queue &) = delete;
	Queue& operator=(Queue &&) = delete;

	void reset() {
		this->start = this->capacity - 1;
		this->end = this->start;
	}

	void push_front(T value) {
		this->array[this->start] = std::move(value);
    this->start--;
	}

	void push_back(T value) {
    this->end++;
		this->array[this->end] = std::move(value);
	}

	void pop_front() {
		this->start++;
	}

	void pop_back() {
		this->end--;
	}

	T peek_front() const {
		return this->array[this->start + 1];
	}

	T peek_back() const {
		return this->array[this->end];
	}

	bool is_empty() const {
		return this->start == this->end;
	}

	void print() {
		for(unsigned int i = start + 1; i <= end; i++) {
			std::cout << this->array[i] << ", ";
		}
		std::cout << std::endl;
	}
};

struct CRDataSemiExplicit final {
	size_t k;
	float a;
	float d;
};

// Important: intervals must have size at least as large as entry_costs, even though one entry less is required.
// The last entry must not contain -1 as first component, Just set it to [0, 1] for simplicity
//
// first component of exit_costs is value, second is the index of entry_costs it came from (only if the exit_cost is not INDEX_UNREACHABLE)
template <typename T, std::pair<T, T> EMPTY_INTERVAL>
static void cell_reachability_explicit(size_t const *entry_costs, std::pair<T, T> const *intervals,
																			 std::pair<size_t, size_t > *exit_costs, Queue<std::tuple<T, size_t, size_t>> &queue) {
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
		size_t ref = j;
		while(!queue.is_empty()) {
			auto const [t, k, ref_index] = queue.peek_front();
			if (k < lambda_j && t > last_interval.first) break;
			if (k < k_left) {
				k_left = k;
				ref = ref_index;
			}
			queue.pop_front();
		}
		queue.push_front(std::tuple<T, size_t, size_t>(last_interval.first, k_left, ref));

		while(!queue.is_empty()) {
			T const t = std::get<0>(queue.peek_back());
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

// EMPTY_INTERVAL must satisfy that its first component is the value to be used for invalid entries, so it should not be a valid entry
template <typename T, std::pair<T, T> eq_solver(Point const &, Point const &, Point const &, float),
	std::pair<T, T> NON_EMPTY_INTERVAL, std::pair<T, T> EMPTY_INTERVAL>
struct Simplifier final {
	Polyline const &polyline;
	_DP1Data<T> *dp1;
	KappaData *kappa;
	KappaData *kappa2;


	std::pair<T, T> *reachability_data;
	std::pair<T, T> *intervals;
	size_t *entry_costs;
	std::pair<size_t, size_t> *exit_costs;
	Queue<std::tuple<T, size_t, size_t>> queue;
	float const epsilon;
	size_t const n;


	// NOTE: first coordinate is now i, second j, third k
	Simplifier(Polyline const &polyline, float epsilon) 
		: polyline(polyline), queue(polyline.point_count), epsilon(epsilon), n(polyline.point_count) {
		this->reachability_data = new std::pair<T, T>[n * (n - 1)];
		this->entry_costs = new size_t[n];
		this->exit_costs = new std::pair<size_t, size_t>[n];
		this->intervals = new std::pair<T, T>[n];
		this->dp1 = new _DP1Data<T>[n * (n - 1) * n];
		this->kappa = new KappaData[n * (n - 1)];
		this->kappa2 = new KappaData[n * (n - 1)];
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
				this->reachability_data[index++] = eq_solver(polyline.get_point(j), polyline.get_point(j+1), polyline.get_point(i), epsilon);
			}
		}

		size_t j = 1;
		for(; j < n - 1 && reachability_data[j].first == 0; j++);
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

	inline T dp(size_t k, size_t i, size_t j) {
		if (k == 0 && i > 0) {
			return EMPTY_INTERVAL.first;
		}
		size_t const pos = (n - 1) * i + j;
		auto const range = reachability_data[pos];
		if (k >= kappa2[pos].smallest_k) {
			return range.first;
		} 
		float const t = dp1[pos * n +  k].first_reachable_point;
		if (t != EMPTY_INTERVAL.first && t <= range.second) {
			return std::max(t, range.first);
		}
		return EMPTY_INTERVAL.first;
	}

	inline void kappa2_subroutine(size_t i) {
		// size_t temp = 0;
		for (size_t i_ = 0; i_ < i; i_++) {
			for (size_t j = 0; j < n - 1; j++) {
				this->intervals[j] = eq_solver(polyline.get_point(i_), polyline.get_point(i), polyline.get_point(j+1), epsilon);

				auto const k = kappa[j + (n-1) * i_].smallest_k;
				this->entry_costs[j] = DataStructures::INDEX_UNREACHABLE;
				if (k != DataStructures::INDEX_UNREACHABLE) {
					this->entry_costs[j] = k + 1;
				}
			}

			// auto const k = kappa[n - 2 + (n-1) * i_].smallest_k;
			// this->entry_costs[n-1] = DataStructures::INDEX_UNREACHABLE;
			// if (k != DataStructures::INDEX_UNREACHABLE) {
			// 	this->entry_costs[n-1] = k + 1;
			// }

			cell_reachability_explicit<T, EMPTY_INTERVAL>(const_cast<size_t const *>(entry_costs), 
																								 const_cast<std::pair<T, T> const *>(intervals), exit_costs, queue);
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
	}

  inline void simplify() {
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
					} else if(v1.first_reachable_point == EMPTY_INTERVAL.first || v1.first_reachable_point > v2) {
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

		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n - 1; j++) {
				std::cout << i << "," << j << ": " << kappa[j + (n-1)*i].smallest_k << "\n";
				std::cout << i << "," << j << ": " << kappa2[j + (n-1)*i].smallest_k << "\n";
			
			}
		}

		std::cout << "simplification size is " << kappa[(n-1) * n - 1].smallest_k << std::endl;
	}
};



Simplification::Simplification 
Simplification::simplification_advanced_euclidean_explicit(DataStructures::Polyline &polyline, float epsilon) {
	Simplifier<float, DataStructures::solve_euclidean, std::pair<float, float>(0, 1), DataStructures::EMPTY_INTERVAL_EXPLICIT> 
		simplifier(polyline, epsilon);
	simplifier.simplify();
  return std::make_unique<std::vector<size_t>>();
}
