#ifndef INCLUDE_DATAGENERATION_GENERATORS_H_
#define INCLUDE_DATAGENERATION_GENERATORS_H_

#include "datastructures.h"
#include <cmath>
#include <filesystem>
#include <memory>
#include <random>

namespace DataGeneration {

std::unique_ptr<DataStructures::Polyline>
make_polyline(size_t, size_t, float, float, float, std::mt19937 &);
void make_integral(DataStructures::Polyline &);
void write_to_file(DataStructures::Polyline &, std::filesystem::path);
} // namespace DataGeneration

#endif // INCLUDE_DATAGENERATION_GENERATORS_H_
