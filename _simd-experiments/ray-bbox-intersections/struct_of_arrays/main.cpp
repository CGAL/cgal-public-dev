
#include "../util.h"

#include "vvector3.h"
#include "vbbox.h"

#include "intersection_strategies/xsimd.h"
//#include "intersection_strategies/branchless.h"

#include <iostream>
#include <numeric>

template<typename T>
struct VQuery {
  Ray<T> ray;
  VBBox<T> vbox;
};

template<typename T>
auto pack_queries(std::vector<Query<T>> queries) {
  std::vector<VQuery<T>> vqueries;

  for (const Query<T> &query : queries)
    vqueries.push_back({query.ray, VBBox<T>(query.boxes)});

  return vqueries;
}

int main() {


  long N = 3742217;
//  long R = 100;

  // Load test data
  auto file = std::ifstream("../data/remeshing_intersections_3742217.txt");
  if (!file.is_open()) return EXIT_FAILURE;
  auto queries = load_queries<double>(file, N);
  std::cout << "Loaded "
            << std::accumulate(queries.begin(), queries.end(), 0,
                               [](const auto &a, const Query<double> &b) { return a + b.boxes.size(); })
            << " scenarios "
            << "divided into " << queries.size() << " queries." << std::endl;

  // Convert test data to vector format
  auto vqueries = pack_queries(queries);

}