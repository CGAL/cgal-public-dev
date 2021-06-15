#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <vector>

#include "xvector3.h"
#include "xray.h"
#include "xbbox.h"


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

