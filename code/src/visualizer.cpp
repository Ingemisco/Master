#include "visualizer.h"
#include "datastructures.h"
#include <fstream>
namespace VisualizationLog {

VisualizationLogger::VisualizationLogger(
    DataStructures::Polyline const &polyline, float epsilon, Distance distance)
    : polyline(polyline), epsilon(epsilon), distance(distance) {}

// creates a range of leaves of the dynamic program for the entries (0, 0, 0) to
// (0, 0, last_reachable)
void VisualizationLogger::leave_range(size_t last_reachable) {
  for (unsigned int j = 0; j < last_reachable; j++) {
    this->data.push_back(VisualizationData(0, 0, j, 0, 0, 0));
  }
}

// sets the dynamic program table of (k, i, j) to use the value of (k - 1, i_,
// j_), sets first reachable on (j,j+1)  to t
void VisualizationLogger::add_use(VisualizationData d) {
  this->data.push_back(d);
}

void VisualizationLogger::emit() {
  // Get the current time
  auto now = std::chrono::system_clock::now();
  auto now_time_t = std::chrono::system_clock::to_time_t(now);
  auto now_tm = *std::localtime(&now_time_t);

  // Format the timestamp
  std::ostringstream oss;
  oss << std::put_time(&now_tm, "%Y-%m-%d-%H-%M-%S");
  auto filename = "vislog-" + oss.str() + ".log";

  std::ofstream out(filename);

  switch (this->distance) {
  case Distance::EUCLIDEAN:
    out << "E ";
    break;
  case Distance::MANHATTAN:
    out << "M ";
    break;
  case Distance::CHEBYSHEV:
    out << "C ";
    break;
  }

  out << this->epsilon << " " << this->polyline.point_count << " "
      << this->polyline.dimension << "\n";
  for (unsigned int i = 0; i < this->polyline.point_count; i++) {
    for (unsigned int j = 0; j < this->polyline.dimension; j++) {
      out << this->polyline[i, j] << " ";
    }
    out << "\n";
  }

  for (auto &d : this->data) {
    out << d.k << " " << d.i << " " << d.j << " " << d.i_ << " " << d.j_ << " "
        << d.t << "\n";
  }
}

VisualizationData::VisualizationData(size_t k, size_t i, size_t j, size_t i_,
                                     size_t j_, float t)
    : k(k), i(i), j(j), i_(i_), j_(j_), t(t) {}

} // namespace VisualizationLog
