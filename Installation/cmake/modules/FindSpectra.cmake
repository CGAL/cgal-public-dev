# - Try to find Spectra
# Once done this will define
#
#  SPECTRA_FOUND = SPECTRA_FOUND - TRUE
#  SPECTRA_INCLUDE_DIR - include directory for Spectra

if(SPECTRA_INCLUDE_DIR)
  # in cache already
  set(SPECTRA_FOUND TRUE)
  set(SPECTRA_INCLUDE_DIRS "${SPECTRA_INCLUDE_DIR}")
else(SPECTRA_INCLUDE_DIR)

# first look in user defined locations
find_path(SPECTRA_INCLUDE_DIR
          NAMES SymGEigsSolver.h
                MatOp/SparseSymMatProd.h
                MatOp/SparseCholesky.h
          PATHS /usr/local/include/Spectra
                /usr/local/include/Spectra/include 
          ENV   SPECTRA_INC_DIR
         )

set(SPECTRA_FOUND TRUE)
set(SPECTRA_INCLUDE_DIRS "${SPECTRA_INCLUDE_DIR}")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SPECTRA_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SPECTRA  DEFAULT_MSG
                                  SPECTRA_INCLUDE_DIR)

mark_as_advanced(SPECTRA_INCLUDE_DIR)


endif(SPECTRA_INCLUDE_DIR)

