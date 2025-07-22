#ifndef INCLUDE_INCLUDE_DATASTRUCTURES_H_
#define INCLUDE_INCLUDE_DATASTRUCTURES_H_

#include <cstddef>
#include <filesystem>
#include <iostream>

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


template<typename T>
class Queue final {
private:
	T *array;
	size_t start; // position **before** first element, i.e., index where push_front would insert 
	size_t end; // position last element 

public:
	size_t const capacity;

	Queue(size_t capacity) : start(capacity - 1), end(capacity - 1), capacity(capacity) {
		// NOTE: queue starts in the middle and has to both sides enough space to avoid modulo computations/wrap around
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

} // namespace DataStructures

#endif // INCLUDE_INCLUDE_DATASTRUCTURES_H_
