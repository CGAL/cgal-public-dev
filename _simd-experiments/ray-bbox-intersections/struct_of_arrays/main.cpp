
#include "../util.h"

#include "vvector3.h"
#include "vbbox.h"

#include <fstream>

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

}