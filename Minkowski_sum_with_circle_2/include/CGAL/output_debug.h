#ifndef MINKOWSKI_SUM_2_DEBUG_H
#define MINKOWSKI_SUM_2_DEBUG_H

#define SHOULD_LOG // comment this definition to see no logging
#ifdef SHOULD_LOG
#define LOG_DEBUG std::clog
#else
#define LOG_DEBUG if(false) std::clog
#endif

#define SHOULD_OUT // comment this definition to see no output
#ifdef SHOULD_OUT
#define OUT_DEBUG if(Base::verbose()) std::cout
#else
#define OUT_DEBUG if(false) std::cout
#endif

#ifdef SHOULD_OUT
#define LIM_DEBUG std::cout
#else
#define LIM_DEBUG if(false) std::cout
#endif

#define SHOULD_ERR // comment this definition to see no errors
#ifdef SHOULD_ERR
#define ERR_DEBUG std::cerr
#else
#define ERR_DEBUG if(false) std::cerr
#endif

#endif // MINKOWSKI_SUM_2_DEBUG_H
