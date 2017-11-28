#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#define PN "\r\n"
#else 
#define PS "/" 
#define PN "\n"
#endif

// Google test includes.
#include "gmock/gmock.h"

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

using namespace testing;

class MyLogTest: public Test {

public:
	using Log = CGAL::LOD::Mylog;
	Log log;
};

TEST_F(MyLogTest, PrintsOk) {
   
   ASSERT_THAT(log.state(), Eq("ok")); 
}

TEST_F(MyLogTest, AddsIntValueToTheLogStream) {

	auto value = 100;
	log.out << "test value: " << value;

	ASSERT_THAT(log.data(), Eq("test value: 100"));
}

TEST_F(MyLogTest, AppendsDoubleValueToTheLogObject) {

	auto value = 24.2;
	log.append(value);

	ASSERT_THAT(log.data(), Eq("24.2"));
}

TEST_F(MyLogTest, SavesMultipleLinesToTheFile) {

	auto intValue = 25;
	auto doubleValue = 34.5;
	auto stringValue = ": double";

	log.out << intValue    << "" + std::string(PN) + ""
			<< doubleValue << " "
			<< stringValue << std::endl;

	const auto saved = true;

	ASSERT_THAT(log.save("test"), Eq(saved));
}

TEST_F(MyLogTest, ClearsData) {
	
	log.out << "some data" << std::endl;
	log.append("some other data" + std::string(PN) + "");

	log.clear();

	ASSERT_THAT(log.data(), Eq(""));
}