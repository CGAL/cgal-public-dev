#include "gmock/gmock.h"
#include <CGAL/LOD_log.h>

using namespace testing;

class LogTest: public Test {

public:
   CGAL::LOD_log log;
};

TEST_F(LogTest, PrintsOk) {
   ASSERT_THAT(log.print(), Eq("ok")); 
}