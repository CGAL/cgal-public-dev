//    (c) 2008-2009 National and Kapodistrian University of Athens
//    (c) 2009-2011 INRIA Nancy
//    (c) 2011-2012 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED

#ifndef OBJECT_CACHE_H
#define OBJECT_CACHE_H

#include<map>
//#include<boost/unordered_map.hpp>

#include <CGAL/OpenMP.h>

#ifdef CGAL_OBJECT_CACHE_DISABLE

template<class K, class V>
struct Value_cache {
    typedef std::pair<bool, V> Found_value;

    Found_value read(const K k) const {
        return std::make_pair(false, V());
    }

    void write(const K k, const V& obj) { }
    void lock() { }
    void unlock() { }
};

template<class B, class K = typename B::Key_type>
struct Object_cache: Value_cache<K, B*> {
    B* read(const B& obj) { return 0; }
    void write(const B& obj) { }
    void lock() { }
    void unlock() { }
};

#else

template<class K, class V>
class Value_cache {
protected:
    typedef K Key;
    typedef ::std::map<K, V> Cache;
//    typedef ::boost::unordered::unordered_map<K, V> Cache;

    typedef typename Cache::const_iterator Cache_iter;
    static Cache cache;
    static OMP_nest_lock omp_mutex;

public:
    typedef std::pair<bool, V> Found_value;

    Found_value read(const K k) const {
        OMP_guard<OMP_nest_lock> x_(omp_mutex);
        Cache_iter it = cache.find(k);
        if (it == cache.end()) return std::make_pair(false, V());
        return std::make_pair(true, it->second);
    }

    void write(const K k, const V& obj) {
        OMP_guard<OMP_nest_lock> x_(omp_mutex);
        cache[k] = obj;
    }

#ifdef _OPENMP
    void lock() { omp_mutex.lock(); }
    void unlock() { omp_mutex.unlock(); }
#else
    void lock() { }
    void unlock() { }
#endif
};

template<class B, class K = typename B::Key_type>
class Object_cache: Value_cache<K, B*> {
protected:
    typedef Value_cache<K, B*> Base;
    typedef typename Base::Key Key;
    typedef typename Base::Found_value Found_value;

public:    

    B* read(const B& obj) {
        OMP_guard<OMP_nest_lock> x_(Base::omp_mutex);
        Found_value fv = Base::read(obj.key());
        return fv.second;
    }

    void write(const B& obj) {
        OMP_guard<OMP_nest_lock> x_(Base::omp_mutex);
        Key k = obj.key();
        Found_value fv = Base::read(k);
        if (fv.first) delete fv.second;
        Base::cache[k] = new B(obj);
    }

#ifdef _OPENMP
    void lock() { Base::lock(); }
    void unlock() { Base::unlock(); }
#else
    void lock() { }
    void unlock() { }
#endif

};

template<class K, class V>
typename Value_cache<K,V>::Cache Value_cache<K,V>::cache;

//template<class B, class K>
//typename Object_cache<K,V>::Cache Value_cache<K,V>::cache;

template<class K, class V>
OMP_nest_lock Value_cache<K,V>::omp_mutex;

#endif

#endif // OBJECT_CACHE_H
