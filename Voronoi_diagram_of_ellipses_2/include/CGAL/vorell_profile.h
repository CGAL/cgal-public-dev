//    (c) 2007-2009 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED

#ifndef CGAL_VORELL_PROFILE_H
#define CGAL_VORELL_PROFILE_H

#ifdef PROFILE

#include <CGAL/Timer.h>
#define PROF_DECL CGAL::Timer timer
#define PROF_START timer.start()
#define PROF(msg) std::cerr << msg << timer.time() << std::endl
#define PROF_STOP timer.stop
#define PROF_RESET timer.reset

#else

#define PROF_DECL
#define PROF_START
#define PROF(msg)
#define PROF_STOP
#define PROF_RESET

#endif

#endif
