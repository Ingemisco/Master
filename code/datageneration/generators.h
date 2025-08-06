#ifndef INCLUDE_DATAGENERATION_GENERATORS_H_
#define INCLUDE_DATAGENERATION_GENERATORS_H_

#include "datastructures.h"
#include "global.h"
#include <filesystem>
#include <memory>
#include <random>

namespace DataGeneration {

std::unique_ptr<DataStructures::Polyline>
make_polyline(PointCount, Dimension, float, float, float, std::mt19937 &);
void make_integral(DataStructures::Polyline &);
void write_to_file(DataStructures::Polyline const &, std::filesystem::path);

bool make_test_suite(std::string &, size_t, size_t, size_t, size_t, std::mt19937 &);

} // namespace DataGeneration

#endif // INCLUDE_DATAGENERATION_GENERATORS_H_
