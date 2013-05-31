#include <iostream>
#include <sstream>
#include <fenv.h>
#include <csignal>

class Floating_point_exception {
  const int signal;
public:
  Floating_point_exception (int signal)
    : signal (signal) {
    std::clog << "Constructed " << *this << std::endl;
  }
  const char * what() const {
    switch (signal) {
    case (FE_DIVBYZERO):
      return "Floating point exception FE_DIVBYZERO";
    case (FE_INVALID):
      return "Floating point exception FE_INVALID";
    case (FE_OVERFLOW):
      return "Floating point exception FE_OVERFLOW";
    default:
      return "Floating point exception (unknown signal)";
    }
  }
  friend std::ostream & operator<< (std::ostream &out,
                                    const Floating_point_exception &e) {
    return out << e.what();
  }
};

void throw_fpe (int signal) { throw Floating_point_exception (signal); }

class Floating_point_exception_handler {
  fenv_t old_environment;
public:
  Floating_point_exception_handler () {
    feholdexcept (&old_environment);
    //feenableexcept (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    feenableexcept (FE_ALL_EXCEPT);
    signal (SIGFPE, throw_fpe);
    std::cerr << "Installed floating point exception handling" << std::endl;
  }
  ~Floating_point_exception_handler () {
    fesetenv (&old_environment);
    std::cerr << "Cleared floating point exception handling"
              << " (to: " << fegetexcept() << ")" << std::endl;
  }
};

double f (double x, double y) {
  Floating_point_exception_handler handler;
  // std::cerr << "fegetexcept() at start of f: " << fegetexcept() << " should be "
  //           << FE_ALL_EXCEPT << " = FE_ALL_EXCEPT" << std::endl;
  try {
    //if (y == 0.) throw Floating_point_exception (FE_OVERFLOW);

    volatile double res = (x/y);
    return res;
  } catch (...) {
    std::cerr << "Caught in f." << std::endl;
    return 100000.;
  }
}

double g (double x, double y) {
  try {
    return f (x,y);
  } catch (...) {
    std::cerr << "Caught in g." << std::endl;
    return 100000.;
  }
}

int main () {
  //Floating_point_exception_handler handler;
  using namespace std;
  double x, y;
  string X = "1.0e100";
  string Y = "0.0";
  istringstream (X) >> x;
  istringstream (Y) >> y;
  try {
    cout << "g (" << x << ", " << y << ") = " << g (x, y) << endl;
  } catch (...) {
    std::cerr << "Caught in main." << std::endl;
    cout << "g (" << x << ", " << y << ") = " << 100000. << endl;
  }

  return (int)(g (0., 1.));
}
