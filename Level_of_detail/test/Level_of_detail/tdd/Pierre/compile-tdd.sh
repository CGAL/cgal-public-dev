COMPNAME=$1

cd tdd-build
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug ../cgal-public-dev/Level_of_detail/test/Level_of_detail/tdd/$COMPNAME/

cd ..