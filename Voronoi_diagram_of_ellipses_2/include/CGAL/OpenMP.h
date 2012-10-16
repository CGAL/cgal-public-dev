//    (c) 2011-2012 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED

#ifndef OPENMP_H
#define OPENMP_H

#ifdef _OPENMP

#include <omp.h>

class OMP_lock {
    omp_lock_t mutex;
    bool init;
    bool locked;
public:
    OMP_lock(): init(false) { }
    ~OMP_lock() { if (init) omp_destroy_lock(&mutex); }

    void lock() {
//        std::cerr << "lock " << this << std::endl;
        if (!init) {
            omp_init_lock(&mutex); init = true;
        }
        if (omp_in_parallel()) {
            omp_set_lock(&mutex);
            locked = true;
        } else locked = false;
    }
    void unlock() {
//        std::cerr << "unlock " << this << std::endl;
        if (locked) omp_unset_lock(&mutex);
    }
};

class OMP_nest_lock {
    omp_nest_lock_t mutex;
    bool init;
    int ignored; // TODO: fixme, needs stack
public:
    OMP_nest_lock(): init(false) { }
    ~OMP_nest_lock() { if (init) omp_destroy_nest_lock(&mutex); }

    void lock() {
//        std::cerr << "Nlock " << this << std::endl;
        if (!init) {
            omp_init_nest_lock(&mutex); init = true;
            ignored = 0;
        }
        if (omp_in_parallel()) {
            omp_set_nest_lock(&mutex);
        } else ignored++;
    }
    void unlock() {
//        std::cerr << "Nunlock " << this << std::endl;
        if (ignored) ignored--;
        else omp_unset_nest_lock(&mutex);
    }
};

template<class T>
class OMP_guard {
    T *m;
public:
    OMP_guard(T &mutex): m(&mutex) { m->lock(); }
    ~OMP_guard() { m->unlock(); }
};

#else

struct OMP_lock {
    void lock() const { }
    void unlock() const { }
};

struct OMP_nest_lock {
    void lock() const { }
    void unlock() const { }
};

template<class T=OMP_lock>
struct OMP_guard {
    OMP_guard(const T &mutex) { }
};

#define omp_in_parallel() (0)

#endif

#endif // OPENMP_H
