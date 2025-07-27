#ifndef INCLUDE_GLOBAL_H
#define INCLUDE_GLOBAL_H

#include <ostream>

typedef size_t LineStart;
typedef size_t LineEnd;
typedef size_t Dimension;
typedef size_t PointCount;
typedef size_t PointIndex;
typedef size_t Coordinate;

template<typename Tuple, std::size_t... Is>
inline void print_tuple(std::ostream &os, const Tuple &t, std::index_sequence<Is...>) {
    ((os << (Is == 0 ? "" : ", ") << std::get<Is>(t)), ...);
}

// Generalized operator<< for any tuple
template<typename... Args>
inline std::ostream &operator<<(std::ostream &os, const std::tuple<Args...> &t) {
    os << "(";
    print_tuple(os, t, std::index_sequence_for<Args...>{});
    os << ")";
    return os;
}

template<typename A, typename B>
inline std::ostream &operator<<(std::ostream &os, const std::pair<A, B> &t) {
	os << "(" << t.first << ", " << t.second <<")";
	return os;
}






#endif // INCLUDE_GLOBAL_H
