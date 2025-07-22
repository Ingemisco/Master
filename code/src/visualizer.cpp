#include "visualizer.h"
#include "datastructures.h"
#include "distance.h"
#include <cstdlib>
#include <fstream>
namespace VisualizationLog {

VisualizationLogger::VisualizationLogger(
    DataStructures::Polyline const &polyline, float epsilon, Distance distance)
    : polyline(polyline), epsilon(epsilon), distance(distance) {}

// creates a range of leaves of the dynamic program for the entries (0, 0, 0) to
// (0, 0, last_reachable)
void VisualizationLogger::leave_range(size_t last_reachable) {
	switch (this->distance) {
		case Distance::EUCLIDEAN_IMPLICIT:
		for (unsigned int j = 0; j < last_reachable; j++) {
			this->data.push_back(VisualizationData(0, 0, j, 0, 0, 0lu));
		}
		break;

		case Distance::EUCLIDEAN_SEMIEXPLICIT:
		for (unsigned int j = 0; j < last_reachable; j++) {
			this->data.push_back(VisualizationData(0, 0, j, 0, 0, DataStructures::FRValue(0,0)));
		}
		break;

		default: // explicit 
		for (unsigned int j = 0; j < last_reachable; j++) {
			this->data.push_back(VisualizationData(0, 0, j, 0, 0, 0.0f));
		}
		break;
	}
}

// sets the dynamic program table of (k, i, j) to use the value of (k - 1, i_,
// j_), sets first reachable on (j,j+1)  to t
void VisualizationLogger::add_use(VisualizationData d) {
  this->data.push_back(d);
}

void VisualizationLogger::emit() {
  if (!std::filesystem::exists("visualization")) {
    std::filesystem::create_directory("visualization");
  }
  if (!std::filesystem::is_directory("visualization")) {
    throw std::runtime_error(
        "File called 'visualization' exists. Cannot create "
        "a directory of same name.");
  }
  // Get the current time
  auto now = std::chrono::system_clock::now();
  auto now_time_t = std::chrono::system_clock::to_time_t(now);
  auto now_tm = *std::localtime(&now_time_t);

  // Format the timestamp
  std::ostringstream oss;
  oss << std::put_time(&now_tm, "%Y-%m-%d-%H-%M-%S-") << rand();
  auto filename = "visualization/" + oss.str() + ".log";

  std::ofstream out(filename);

  switch (this->distance) {
  case Distance::EUCLIDEAN:              out << "E ";  break;
  case Distance::MANHATTAN:              out << "M ";  break;
  case Distance::CHEBYSHEV:              out << "C ";  break;
  case Distance::EUCLIDEAN_IMPLICIT:     out << "IE "; break;
  case Distance::EUCLIDEAN_SEMIEXPLICIT: out << "SE "; break;
  }

  out << this->epsilon << " " << this->polyline.point_count << " " << this->polyline.dimension << "\n";
  for (unsigned int i = 0; i < this->polyline.point_count; i++) {
    for (unsigned int j = 0; j < this->polyline.dimension; j++) {
      out << this->polyline[i, j] << " ";
    }
    out << "\n";
  }

	switch (this->distance) {
		case Distance::EUCLIDEAN_IMPLICIT:
		for (auto &d : this->data) {
			out << d.k << " " << d.i << " " << d.j << " " << d.i_ << " " << d.j_ << " " << d.r << "\n";
		}
		break;

		case Distance::EUCLIDEAN_SEMIEXPLICIT:
		for (auto &d : this->data) {
			out << d.k << " " << d.i << " " << d.j << " " << d.i_ << " " << d.j_ << " " << d.f << "\n";
		}
		break;

		default: // explicit
		for (auto &d : this->data) {
			out << d.k << " " << d.i << " " << d.j << " " << d.i_ << " " << d.j_ << " " << d.t << "\n";
		}
		break;
	}
  out.close();
}

VisualizationData::VisualizationData(size_t k, size_t i, size_t j, size_t i_, size_t j_, float t)
    : k(k), i(i), j(j), i_(i_), j_(j_), t(t) {}

VisualizationData::VisualizationData(size_t k, size_t i, size_t j, size_t i_, size_t j_, size_t r)
    : k(k), i(i), j(j), i_(i_), j_(j_), r(r) {}

VisualizationData::VisualizationData(size_t k, size_t i, size_t j, size_t i_, size_t j_, DataStructures::FRValue f)
    : k(k), i(i), j(j), i_(i_), j_(j_), f(f) {}

} // namespace VisualizationLog
