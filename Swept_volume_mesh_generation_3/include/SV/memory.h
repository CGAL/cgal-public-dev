// Copyright (c) 2011 Andreas von Dziegielewski and Michael Hemmer (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
//
// Author(s)     : Michael Hemmer (mhsaar@googlemail.com)
//
// ================================================================================

#ifndef SV_MEMORY_H
#define SV_MEMORY_H

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include <CGAL/Memory_sizer.h>

namespace SV{

// return memory usage in % 
double inline memory_usage() {
  rlimit rlim;  
  //rusage usag;
  getrlimit(RLIMIT_RSS, &rlim);
  //getrusage(RUSAGE_SELF, &usag); 

  /// this is a bad hack, it will just work for the server at mainz and my laptop .-) 
  double max_memory = (std::min)(double(rlim.rlim_cur), double(858993459)*double(60));
  CGAL::Memory_sizer msizer; 
  //std::cout << "\n" << msizer.virtual_size () << " " << rlim.rlim_cur<< " " << max_memory << std::endl; 
  return double(msizer.virtual_size ())/ max_memory; 
  //return msizer.resident_size () /double(rlim.rlim_cur);
  //  return (double(usag.ru_maxrss) *1024 )/(double(rlim.rlim_cur ));
}

double inline memory_usage_max() {
  rlimit rlim;  
  rusage usag;
  getrlimit(RLIMIT_RSS, &rlim);
  getrusage(RUSAGE_SELF, &usag); 
  //CGAL::Memory_sizer msizer; 
  
  return (double(usag.ru_maxrss) *1024 )/(double(rlim.rlim_cur ));
  //return msizer.virtual_size ()/double(rlim.rlim_cur); 
  //return msizer.resident_size () /double(rlim.rlim_cur);
}

}//namespace SV 
#endif // SV_MEMORY_H
