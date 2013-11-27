include(CGAL_library)

CGAL_define_library(
  NAME CGAL_Core
  BUILD_DEFAULT ON
  REQUIRED_DEPENDENCIES CGAL GMP MPFR
  OPTIONAL_DEPENDENCIES
  SOURCES BigFloat.cpp CoreAux.cpp CoreDefs.cpp CoreIO.cpp Expr.cpp extLong.cpp GmpIO.cpp Real.cpp)
