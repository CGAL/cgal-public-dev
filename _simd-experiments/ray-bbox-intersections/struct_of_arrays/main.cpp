
#include "vvector3.h"
#include "vbbox.h"
#include "vray.h"

auto load(std::ifstream &file, long N) {

}

double time(const std::function<void(void)> &f) {

  auto start = std::chrono::_V2::high_resolution_clock::now();
  {
    f();
  }
  auto end = std::chrono::_V2::high_resolution_clock::now();

  return (end - start).count();
}

int main() {

}