// Copyright (c) 1999,2001,2003  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra, Sylvain Pion
 
#ifndef CGAL_HANDLE_FOR_H
#define CGAL_HANDLE_FOR_H

#include <boost/config.hpp>
#include <CGAL/memory.h>
#include <CGAL/assertions.h>
#include <algorithm>
#include <cstddef>
#include <vector>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4345) // Avoid warning  http://msdn.microsoft.com/en-us/library/wewb47ee(VS.80).aspx
#endif
namespace CGAL {

// The recycle parameter has to go last for compatibility :-(
// recycle == 0: don't recycle anything
// recycle == 1: recycle the allocated space
// recycle == 2: recycle the objects (don't call the destructor)
template <class T, class Alloc = CGAL_ALLOCATOR(T), int recycle=0>
class Handle_for
{
    // Wrapper that adds the reference counter.
    struct RefCounted {
        T t;
        unsigned int count;
    };

    typedef typename Alloc::template rebind<RefCounted>::other  Allocator;
    typedef typename Allocator::pointer                         pointer;

    static Allocator   allocator;
    pointer            ptr_;

    static std::vector<RefCounted*> pool;

public:

    typedef T element_type;
    
    typedef std::ptrdiff_t Id_type ;

    // With recycle == 2, this may return an old object with a
    // non-default value. We may want to use a different constructor
    // Handle_for(struct Whatever) for that.
    //
    // This constructor is the only one used by Gmpq (besides copy/move)
    Handle_for()
    {
	if (recycle>=1 && !pool.empty()) {
		ptr_=pool.back();
		pool.pop_back();
		if (recycle<=1) 
			new (&(ptr_->t)) element_type();
	} else {
		ptr_=allocator.allocate(1);
		new (&(ptr_->t)) element_type(); // we get the warning here 
	}
        ptr_->count = 1;
    }

    Handle_for(const element_type& t)
    {
	if (recycle>=1 && !pool.empty()) {
		ptr_=pool.back();
		pool.pop_back();
		if (recycle<=1) 
			new (&(ptr_->t)) element_type(t);
		else ptr_->t = t; // recycle == 2 requires assignable here
	} else {
		ptr_=allocator.allocate(1);
		new (&(ptr_->t)) element_type(t);
	}
        ptr_->count = 1;
    }

#ifndef CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
    Handle_for(element_type && t)
    {
	if (recycle>=1 && !pool.empty()) {
		ptr_=pool.back();
		pool.pop_back();
		if (recycle<=1) 
			new (&(ptr_->t)) element_type(std::move(t));
		else ptr_->t = std::move(t);
	} else {
		ptr_=allocator.allocate(1);
		new (&(ptr_->t)) element_type(std::move(t));
	}
        ptr_->count = 1;
    }
#endif

/* I comment this one for now, since it's preventing the automatic conversions
   to take place.  We'll see if it's a problem later.
    template < typename T1 >
    Handle_for(const T1& t1)
      : ptr_(allocator.allocate(1))
    {
        new (&(ptr_->t)) T(t1);
        ptr_->count = 1;
    }
*/

#if !defined CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES && !defined CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
    template < typename T1, typename T2, typename... Args >
    Handle_for(T1 && t1, T2 && t2, Args && ... args)
    {
	if (recycle>=1 && !pool.empty()) {
		ptr_=pool.back();
		pool.pop_back();
		if (recycle==2) allocator.destroy(ptr_);
	} else {
		ptr_=allocator.allocate(1);
	}
	new (&(ptr_->t)) element_type(std::forward<T1>(t1), std::forward<T2>(t2), std::forward<Args>(args)...);
        ptr_->count = 1;
    }
#else
    template < typename T1, typename T2 >
    Handle_for(const T1& t1, const T2& t2)
    {
	if (recycle>=1 && !pool.empty()) {
		ptr_=pool.back();
		pool.pop_back();
		if (recycle==2) allocator.destroy(ptr_);
	} else {
		ptr_=allocator.allocate(1);
	}
        new (&(ptr_->t)) element_type(t1, t2);
        ptr_->count = 1;
    }

    template < typename T1, typename T2, typename T3 >
    Handle_for(const T1& t1, const T2& t2, const T3& t3)
    {
	if (recycle>=1 && !pool.empty()) {
		ptr_=pool.back();
		pool.pop_back();
		if (recycle==2) allocator.destroy(ptr_);
	} else {
		ptr_=allocator.allocate(1);
	}
        new (&(ptr_->t)) element_type(t1, t2, t3);
        ptr_->count = 1;
    }

    template < typename T1, typename T2, typename T3, typename T4 >
    Handle_for(const T1& t1, const T2& t2, const T3& t3, const T4& t4)
    {
	if (recycle>=1 && !pool.empty()) {
		ptr_=pool.back();
		pool.pop_back();
		if (recycle==2) allocator.destroy(ptr_);
	} else {
		ptr_=allocator.allocate(1);
	}
        new (&(ptr_->t)) element_type(t1, t2, t3, t4);
        ptr_->count = 1;
    }
#endif // CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

    Handle_for(const Handle_for& h)
      : ptr_(h.ptr_)
    {
        ++(ptr_->count);
    }

    Handle_for&
    operator=(const Handle_for& h)
    {
        Handle_for tmp = h;
        swap(tmp);
        return *this;
    }

    Handle_for&
    operator=(const element_type &t)
    {
        if (is_shared())
            Handle_for(t).swap(*this);
        else
            ptr_->t = t;

        return *this;
    }

#ifndef CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
    // Note : I don't see a way to make a useful move constructor, apart
    //        from e.g. using NULL as a ptr value, but this is drastic.

    Handle_for&
    operator=(Handle_for && h)
    {
        swap(h);
        return *this;
    }

    Handle_for&
    operator=(element_type && t)
    {
        if (is_shared())
            Handle_for(std::move(t)).swap(*this);
        else
            ptr_->t = std::move(t);

        return *this;
    }
#endif

    ~Handle_for()
    {
      if (--(ptr_->count) == 0) {
          if (recycle<=1) allocator.destroy( ptr_);
          if (recycle==0) allocator.deallocate( ptr_, 1);
	  else pool.push_back( ptr_);
      }
    }

    void
    initialize_with(const element_type& t)
    {
        // kept for backward compatibility.  Use operator=(t) instead.
        *this = t;
    }

    Id_type id() const { return Ptr() - static_cast<T const*>(0); }
    
    bool identical(const Handle_for& h) const { return Ptr() == h.Ptr(); }


    // Ptr() is the "public" access to the pointer to the object.
    // The non-const version asserts that the instance is not shared.
    const element_type *
    Ptr() const
    {
       return &(ptr_->t);
    }

    /*
    // The assertion triggers in a couple of places, so I comment it for now.
    T *
    Ptr()
    {
      CGAL_assertion(!is_shared());
      return &(ptr_->t);
    }
    */

    bool
    is_shared() const
    {
	return ptr_->count > 1;
    }

    bool
    unique() const
    {
	return !is_shared();
    }

    long
    use_count() const
    {
	return ptr_->count;
    }

    void
    swap(Handle_for& h)
    {
      std::swap(ptr_, h.ptr_);
    }

protected:

    void
    copy_on_write()
    {
      if ( is_shared() ) Handle_for(ptr_->t).swap(*this);
    }

    // ptr() is the protected access to the pointer.  Both const and non-const.
    // Redundant with Ptr().
    element_type *
    ptr()
    { return &(ptr_->t); }

    const element_type *
    ptr() const
    { return &(ptr_->t); }
};


template <class T, class Allocator, int r>
typename Handle_for<T, Allocator, r>::Allocator
Handle_for<T, Allocator, r>::allocator;

template <class T, class Allocator, int r>
std::vector<typename Handle_for<T, Allocator, r>::RefCounted*>
Handle_for<T, Allocator, r>::pool;

template <class T, class Allocator, int r>
inline
void
swap(Handle_for<T, Allocator, r> &h1, Handle_for<T, Allocator, r> &h2)
{
    h1.swap(h2);
}

template <class T, class Allocator, int r>
inline
bool
identical(const Handle_for<T, Allocator, r> &h1,
          const Handle_for<T, Allocator, r> &h2)
{
    return h1.identical(h2);
}

template <class T> inline bool identical(const T &t1, const T &t2) { return &t1 == &t2; }

template <class T, class Allocator, int r>
inline
const T&
get(const Handle_for<T, Allocator, r> &h)
{
    return *(h.Ptr());
}

template <class T>
inline
const T&
get(const T &t)
{
    return t;
}

} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_HANDLE_FOR_H
