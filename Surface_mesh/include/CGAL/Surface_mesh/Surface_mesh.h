//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
// Copyright (C) 2014 GeometryFactory
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//


#ifndef CGAL_SURFACE_MESH_H
#define CGAL_SURFACE_MESH_H

#include <iterator>
#include <algorithm>
#include <utility>
#include <iostream>
#include <cstddef>
#include <vector>
#include <string>
#include <typeinfo>

#include <boost/cstdint.hpp>
#include <boost/array.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>
#include <boost/property_map/property_map.hpp>

#include <CGAL/Iterator_range.h>
#include <CGAL/circulator.h>
#include <CGAL/assertions.h>
#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include <CGAL/Surface_mesh/IO.h>
//#include <CGAL/Surface_mesh/Properties.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {
  /// \ingroup PkgSurface_mesh
  /// This class is a data structure that can be used as halfedge data structure or polyhedral
  /// surface. It is an alternative to the classes `HalfedgeDS` and `Polyhedron_3`
  /// defined in the packages  \ref PkgHDSSummary and \ref PkgPolyhedronSummary. 
  /// The main difference is that it is indexed based and not pointer based,
  /// and that the mechanism for adding information to vertices, halfedges,
  /// and faces is much simpler and done at runtime and not at compile time.
  /// When elements are removed, they are only marked as removed, and a garbage
  /// collection function must be called to really remove them. 
  /// @tparam P The type of the \em point property of a vertex. There is no requirement on `P`,
  ///         besides being default constructible and assignable. 
  ///         In typical use cases it will be a 2D or 3D point type.
  /// \cgalModels `MutableFaceGraph` and `FaceListGraph`

template <typename P>
class Surface_mesh
{

    typedef Surface_mesh<P> Self;

    template<typename>
    class Handle_iterator;
public:



/// \addtogroup PkgSurface_mesh
///
/// @{

/// @cond CGAL_DOCUMENT_INTERNALS
class Base_property_array
{
public:

    /// Default constructor
    Base_property_array(const std::string& name) : name_(name) {}

    /// Destructor.
    virtual ~Base_property_array() {}

    /// Reserve memory for n elements.
    virtual void reserve(size_t n) = 0;

    /// Resize storage to hold n elements.
    virtual void resize(size_t n) = 0;

    /// Free unused memory.
    virtual void shrink_to_fit() = 0;

    /// Extend the number of elements by one.
    virtual void push_back() = 0;

    virtual bool transfer(const Base_property_array& other) = 0;

    /// Let two elements swap their storage place.
    virtual void swap(size_t i0, size_t i1) = 0;

    /// Return a deep copy of self.
    virtual Base_property_array* clone () const = 0;

    /// Return the type_info of the property
    virtual const std::type_info& type() = 0;

    /// Return the name of the property
    const std::string& name() const { return name_; }


protected:

    std::string name_;
};

  /// @endcond


//== CLASS DEFINITION =========================================================

/// @cond CGAL_DOCUMENT_INTERNALS

template <class T>
class Property_array : public Base_property_array
{
public:

    typedef T                                       value_type;
    typedef std::vector<value_type>                 vector_type;
    typedef typename vector_type::reference         reference;
    typedef typename vector_type::const_reference   const_reference;

    Property_array(const std::string& name, T t=T()) : Base_property_array(name), value_(t) {}

public: // virtual interface of Base_property_array

    virtual void reserve(size_t n)
    {
        data_.reserve(n);
    }

    virtual void resize(size_t n)
    {
        data_.resize(n, value_);
    }

    virtual void push_back()
    {
        data_.push_back(value_);
    }

    bool transfer(const Base_property_array& other)
    {
      const Property_array<T>* pa = dynamic_cast<const Property_array*>(&other);
      if(pa != NULL){
        std::copy((*pa).data_.begin(), (*pa).data_.end(), data_.end()-(*pa).data_.size());
        return true;
      } 
      return false;
    }

    virtual void shrink_to_fit()
    {
        vector_type(data_).swap(data_);
    }

    virtual void swap(size_t i0, size_t i1)
    {
        T d(data_[i0]);
        data_[i0]=data_[i1];
        data_[i1]=d;
    }

    virtual Base_property_array* clone() const
    {
        Property_array<T>* p = new Property_array<T>(this->name_, this->value_);
        p->data_ = data_;
        return p;
    }

    virtual const std::type_info& type() { return typeid(T); }


public:

    /// Get pointer to array (does not work for T==bool)
    const T* data() const
    {
        return &data_[0];
    }

    /// Access the i'th element. No range check is performed!
    reference operator[](int _idx)
    {
        CGAL_assertion( size_t(_idx) < data_.size() );
        return data_[_idx];
    }

    /// Const access to the i'th element. No range check is performed!
    const_reference operator[](int _idx) const
    {
        CGAL_assertion( size_t(_idx) < data_.size());
        return data_[_idx];
    }



private:
    vector_type data_;
    value_type  value_;
};


#if 0
// specialization for bool properties
template <>
inline const bool*
Property_array<bool>::data() const
{
    CGAL_assertion(false);
    return NULL;
}
#endif 

  /// @endcond

//== CLASS DEFINITION =========================================================

/// @cond CGAL_DOCUMENT_INTERNALS

template<typename>
class Property_container;
/// @endcond 




//== CLASS DEFINITION =========================================================
/// @cond CGAL_DOCUMENT_INTERNALS

template <class, class>
class Property_map;

template<typename Key>
class Property_container
{
public:

    // default constructor
    Property_container() : size_(0) {}

    // destructor (deletes all property arrays)
    virtual ~Property_container() { clear(); }

    // copy constructor: performs deep copy of property arrays
    Property_container(const Property_container& _rhs) { operator=(_rhs); }

    // assignment: performs deep copy of property arrays
    Property_container& operator=(const Property_container& _rhs)
    {
        if (this != &_rhs)
        {
            clear();
            parrays_.resize(_rhs.n_properties());
            size_ = _rhs.size();
            for (unsigned int i=0; i<parrays_.size(); ++i)
                parrays_[i] = _rhs.parrays_[i]->clone();
        }
        return *this;
    }

    void transfer(const Property_container& _rhs)
    {
      for(unsigned int i=0; i<parrays_.size(); ++i){
        for (unsigned int j=0; j<_rhs.parrays_.size(); ++j){
          if(parrays_[i]->name() ==  _rhs.parrays_[j]->name()){
            parrays_[i]->transfer(* _rhs.parrays_[j]);
            break;
          }
        }
      }
    }

    // returns the current size of the property arrays
    size_t size() const { return size_; }

    // returns the number of property arrays
    size_t n_properties() const { return parrays_.size(); }

    // returns a vector of all property names
    std::vector<std::string> properties() const
    {
        std::vector<std::string> names;
        for (unsigned int i=0; i<parrays_.size(); ++i)
            names.push_back(parrays_[i]->name());
        return names;
    }

    // add a property with name \c name and default value \c t
    template <class T>
    std::pair<Property_map<Key, T>, bool>
    add(const std::string& name, const T t=T())
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
        {
            if (parrays_[i]->name() == name)
            {
              return std::make_pair(Property_map<Key, T>(dynamic_cast<Property_array<T>*>(parrays_[i])), false);
            }
        }

        // otherwise add the property
        Property_array<T>* p = new Property_array<T>(name, t);
        p->resize(size_);
        parrays_.push_back(p);
        return std::make_pair(Property_map<Key, T>(p), true);
    }


    // get a property by its name. returns invalid property if it does not exist.
    template <class T> 
    std::pair<Property_map<Key, T>,bool>
    get(const std::string& name) const
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            if (parrays_[i]->name() == name)
              return std::make_pair(Property_map<Key, T>(dynamic_cast<Property_array<T>*>(parrays_[i])), true);
        return std::make_pair(Property_map<Key, T>(), false);
    }


    // returns a property if it exists, otherwise it creates it first.
    template <class T>
    Property_map<Key, T> 
    get_or_add(const std::string& name, const T t=T())
    {
      Property_map<Key, T> p;
      bool b;
      boost::tie(p,b)= get<T>(name);
        if (!b) p = add<T>(name, t).first;
        return p;
    }


    // get the type of property by its name. returns typeid(void) if it does not exist.
    const std::type_info& 
    get_type(const std::string& name)
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            if (parrays_[i]->name() == name)
                return parrays_[i]->type();
        return typeid(void);
    }


    // delete a property
    template <class T> 
    void
    remove(Property_map<Key, T>& h)
    {
        typename std::vector<Base_property_array*>::iterator it=parrays_.begin(), end=parrays_.end();
        for (; it!=end; ++it)
        {
            if (*it == h.parray_)
            {
                delete *it;
                parrays_.erase(it);
                h.reset();
                break;
            }
        }
    }


    // delete all properties
    void clear()
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            delete parrays_[i];
        parrays_.clear();
        size_ = 0;
    }


    // reserve memory for n entries in all arrays
    void reserve(size_t n) const
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            parrays_[i]->reserve(n);
    }

    // resize all arrays to size n
    void resize(size_t n)
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            parrays_[i]->resize(n);
        size_ = n;
    }

    // free unused space in all arrays
    void shrink_to_fit() const
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            parrays_[i]->shrink_to_fit();
    }

    // add a new element to each vector
    void push_back()
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            parrays_[i]->push_back();
        ++size_;
    }

    // swap elements i0 and i1 in all arrays
    void swap(size_t i0, size_t i1) const
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            parrays_[i]->swap(i0, i1);
    }


private:
    std::vector<Base_property_array*>  parrays_;
    size_t  size_;
};

  /// @endcond

#ifndef DOXYGEN_RUNNING
/// 
///
/// `Property_map` enables to attach properties to the simplices of a 
///  surface mesh.
/// 
/// @tparam Key The key type of the property map. It must be a model of `Index`.
/// @tparam Value The value type of the property.
///
/// \cgalModels `LvaluePropertyMap`
///
template <class I, class T>
class Property_map
/// @cond CGAL_DOCUMENT_INTERNALS
  : public boost::put_get_helper< 
           typename Property_array<T>::reference,
           Property_map< I, T > >
/// @endcond
{
    typedef void (Property_map::*bool_type)() const;
    void this_type_does_not_support_comparisons() const {}
public:
    typedef I key_type;
    typedef T value_type;
    typedef boost::lvalue_property_map_tag category;

#ifndef DOXYGEN_RUNNING

    typedef typename Property_array<T>::reference reference;

    typedef typename Property_array<T>::const_reference const_reference;
#else 
    /// A reference to the value type of the property.
  typedef unspecified_type reference;

    /// A const reference to the value type of the property.
  typedef unspecified_type const_reference;
#endif

#ifndef DOXYGEN_RUNNING
    friend class Property_container<I>;

    template <typename K>  friend class Surface_mesh;
#endif

public:
/// @cond CGAL_DOCUMENT_INTERNALS
    Property_map(Property_array<T>* p=NULL) : parray_(p) {}

    void reset()
    {
        parray_ = NULL;
    }
  /// @endcond 

public:
    /// \name Accessing Properties
    //@{
#ifdef DOXYGEN_RUNNING
    /// Conversion to a Boolean. It is \c true when the property map
    /// can be used, and \c false otherwise.  
  operator bool () const;
#else
    operator bool_type() const {
        return parray_ != NULL ?
            &Property_map::this_type_does_not_support_comparisons : 0;
    }
#endif
    /// Access the property associated with the key \c i.
    reference operator[](const I& i)
    {
      CGAL_assertion(parray_ != NULL);
      return (*parray_)[i];
    }

    /// Access the property associated with the key \c i.
    reference operator[](const I& i) const
    {
      CGAL_assertion(parray_ != NULL);
      return (*parray_)[i];
    }

    bool transfer (const Property_map& other)
    {
      return parray_->transfer(*(other.parray_));
    }

    /// Allows access to the underlying storage of the property. This
    /// is useful when the key associated with the properties is
    /// unimportant and only the properties are of interest
    /// (e.g. rendering).
    ///
    /// \returns a pointer to the underlying storage of the property.
    const T* data() const
    {
      CGAL_assertion(parray_ != NULL);
      return parray_->data();
    }

    //@}
private:

    Property_array<T>& array()
    {
        CGAL_assertion(parray_ != NULL);
        return *parray_;
    }

    const Property_array<T>& array() const
    {
        CGAL_assertion(parray_ != NULL);
        return *parray_;
    }

    Property_array<T>* parray_;
};

#endif // DOXYGEN_RUNNING

///@}



    /// \name Basic Types
    ///
    ///@{

    /// The point type.
    typedef P Point;

    /// The type used to represent an index.
    typedef boost::uint32_t size_type;

    ///@}

    /// \name Basic Elements
    ///
    ///@{

#ifndef DOXYGEN_RUNNING
    /// Base class for vertex, halfedge, edge, and face index. 
    ///
    /// \attention Note that `Index` is not a model of the concept `Handle`,
    /// because it cannot be dereferenced.
    /// \sa `Vertex_index`, `Halfedge_index`, `Edge_index`, `Face_index`.
    template<typename T>
    class Index
    {
    public:
        /// Constructor. %Default construction creates an invalid index.
        /// We write -1, which is <a href="http://en.cppreference.com/w/cpp/concept/numeric_limits">
        /// <tt>std::numeric_limits<size_type>::max()</tt></a>
        /// as `size_type` is an unsigned type. 
        explicit Index(size_type _idx=-1) : idx_(_idx) {}

        /// Get the underlying index of this index
        operator size_type() const { return idx_; }

        /// reset index to be invalid (index=-1)
        void reset() { idx_=-1; }

        /// return whether the index is valid, i.e., the index is not equal to -1.
        bool is_valid() const { 
          size_type inf = -1;
          return idx_ != inf;
        }

        /// are two indices equal?
        bool operator==(const T& _rhs) const {
            return idx_ == _rhs.idx_;
        }

        /// are two indices different?
        bool operator!=(const T& _rhs) const {
            return idx_ != _rhs.idx_;
        }

        /// Comparison by index.
        bool operator<(const T& _rhs) const {
            return idx_ < _rhs.idx_;
        }

        /// increments the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// increment.
        Index& operator++() { ++idx_; return *this; }
        /// decrements the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// decrement.
        Index& operator--() { --idx_; return *this; }

        /// increments the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// increment.
        Index operator++(int) { Index tmp(*this); ++idx_; return tmp; }
        /// decrements the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// decrement.
        Index operator--(int) { Index tmp(*this); --idx_; return tmp; }
    private:
        size_type idx_;
    };

#endif

    /// This class represents a vertex.
    /// \cgalModels `Index`
    /// \sa `Halfedge_index`, `Edge_index`, `Face_index`
    class Vertex_index
#ifndef DOXYGEN_RUNNING
 : public Index<Vertex_index>
#endif
    {
    public:
        /// %Default constructor.
        Vertex_index() : Index<Vertex_index>(-1) {}

        explicit Vertex_index(size_type _idx) : Index<Vertex_index>(_idx) {}

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Vertex_index const& v)
        {
          return (os << 'v' << (size_type)v );
        }
    };

    /// This class represents a halfedge.
    /// \cgalModels `Index`
    /// \sa `Vertex_index`, `Edge_index`, `Face_index`
    class Halfedge_index
#ifndef DOXYGEN_RUNNING
      : public Index<Halfedge_index>
#endif
    {
    public:
        /// %Default constructor
        Halfedge_index() : Index<Halfedge_index>(-1) {}

        explicit Halfedge_index(size_type _idx) : Index<Halfedge_index>(_idx) {}

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Halfedge_index const& h)
        {
          return (os << 'h' << (size_type)h );
        }

    };

    /// This class represents a face
    /// \cgalModels `Index`
    /// \sa `Vertex_index`, `Halfedge_index`, `Edge_index`
    class Face_index
#ifndef DOXYGEN_RUNNING
      : public Index<Face_index>
#endif
    {
    public:
        /// %Default constructor
        Face_index() : Index<Face_index>(-1) {}

        explicit Face_index(size_type _idx) : Index<Face_index>(_idx) {}

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Face_index const& f)
        {
          return (os << 'f' << (size_type)f );
        }
    };

    /// This class represents an edge.
    /// \cgalModels `Index`
    /// \sa `Vertex_index`, `Halfedge_index`, `Face_index`
    class Edge_index
    {
    public:
        /// %Default constructor
        Edge_index() : halfedge_(-1) { }

        Edge_index(size_type idx) : halfedge_(idx * 2) { }

        /// constructs an `Edge_index` from a halfedge.
        Edge_index(Halfedge_index he) : halfedge_(he) { }
        /// @cond CGAL_DOCUMENT_INTERNALS
        /// returns the internal halfedge.
        Halfedge_index halfedge() const { return halfedge_; }

        /// returns the underlying index of this index.
        operator size_type() const { return (size_type)halfedge_ / 2; }

        /// resets index to be invalid (index=-1)
        void reset() { halfedge_.reset(); }

        /// returns whether the index is valid, i.e., the index is not equal to -1.
        bool is_valid() const { return halfedge_.is_valid(); }

        /// Are two indices equal?
        bool operator==(const Edge_index& other) const { return (size_type)(*this) == (size_type)other; }

        /// Are two indices different?
        bool operator!=(const Edge_index& other) const { return (size_type)(*this) != (size_type)other; }

        /// compares by index.
        bool operator<(const Edge_index& other) const { return (size_type)(*this) < (size_type)other;}

        /// decrements the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// decrement.
        Edge_index& operator--() { halfedge_ = Halfedge_index((size_type)halfedge_ - 2); return *this; }

        /// increments the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// increment.
        Edge_index& operator++() { halfedge_ = Halfedge_index((size_type)halfedge_ + 2); return *this; }

        /// decrements internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// decrement.
        Edge_index operator--(int) { Edge_index tmp(*this); halfedge_ = Halfedge_index((size_type)halfedge_ - 2); return tmp; }

        /// increments internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// increment.
        Edge_index operator++(int) { Edge_index tmp(*this); halfedge_ = Halfedge_index((size_type)halfedge_ + 2); return tmp; }

        /// @endcond 

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Edge_index const& e)
        {
          return (os << 'e' << (size_type)e << " on " << e.halfedge());
        }
    private:
        Halfedge_index halfedge_;
    };


 
    ///@}
private: //-------------------------------------------------- connectivity types

    /// This type stores the vertex connectivity
    /// \sa `Halfedge_connectivity`, `Face_connectivity`
    struct Vertex_connectivity
    {
        /// an incoming halfedge per vertex (it will be a border halfedge for border vertices)
        Halfedge_index  halfedge_;
    };


    /// This type stores the halfedge connectivity
    /// \sa `Vertex_connectivity`, `Face_connectivity`
    struct Halfedge_connectivity
    {
        /// face incident to halfedge
        Face_index      face_;
        /// vertex the halfedge points to
        Vertex_index    vertex_;
        /// next halfedge within a face (or along a border)
        Halfedge_index  next_halfedge_;
        /// previous halfedge within a face (or along a border)
        Halfedge_index  prev_halfedge_;
    };


    /// This type stores the face connectivity
    /// \sa `Vertex_connectivity`, `Halfedge_connectivity`
    struct Face_connectivity
    {
        /// a halfedge that is part of the face
        Halfedge_index  halfedge_;
    };

private: //------------------------------------------------------ iterator types
    template<typename Index_>
    class Index_iterator
      : public boost::iterator_facade< Index_iterator<Index_>,
                                       Index_,
                                       std::bidirectional_iterator_tag
                                       >
    {
        typedef boost::iterator_facade< Index_iterator<Index_>,
                                        Index_,
                                        std::bidirectional_iterator_tag
                                        > Facade;
    public:
        Index_iterator() : hnd_(), mesh_(NULL) {}
        Index_iterator(const Index_& h, const Surface_mesh* m)
          : hnd_(h), mesh_(m) {
            if (mesh_ && mesh_->has_garbage())
              while (mesh_->has_valid_index(hnd_) && mesh_->is_removed(hnd_)) ++hnd_;
        }
    private:
        friend class boost::iterator_core_access;
        void increment()
        {
            ++hnd_;
            CGAL_assertion(mesh_ != NULL);
            while (mesh_->has_garbage() && mesh_->has_valid_index(hnd_) && mesh_->is_removed(hnd_)) ++hnd_;
        }

        void decrement()
        {
            --hnd_;
            CGAL_assertion(mesh_ != NULL);
            while (mesh_->has_garbage() && mesh_->has_valid_index(hnd_) && mesh_->is_removed(hnd_)) --hnd_;
        }

        bool equal(const Index_iterator& other) const
        {
            return this->hnd_ == other.hnd_;
        }

        Index_& dereference() const { return const_cast<Index_&>(hnd_); }

        Index_ hnd_;
        const Surface_mesh* mesh_;
    };
public:
    /// \name Range Types
    ///
    /// Each range `R` in this section has a nested type `R::iterator`, 
    /// is convertible to `std:pair<R::iterator,R::iterator>`, so that one can use `boost::tie()`,
    /// and can be used with `BOOST_FOREACH()`, as well as with the C++11 range based for-loop.

    ///@{

#ifndef DOXYGEN_RUNNING
    typedef Index_iterator<Vertex_index> Vertex_iterator;
#endif

    /// \brief The range over all vertex indices.
    ///
    /// A model of <a href="http://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a> with value type `Vertex_index`.
    /// \sa `vertices()`
    /// \sa `Halfedge_range`, `Edge_range`, `Face_range`
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type Vertex_range;
#else
    typedef Iterator_range<Vertex_iterator> Vertex_range;
#endif

#ifndef DOXYGEN_RUNNING
    typedef Index_iterator<Halfedge_index> Halfedge_iterator;
#endif

    /// \brief The range over all halfedge indices.
    ///
    /// A model of <a href="http://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a> with value type `Halfedge_index`.
    /// \sa `halfedges()`
    /// \sa `Vertex_range`, `Edge_range`, `Face_range`
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type Halfedge_range;
#else
    typedef Iterator_range<Halfedge_iterator> Halfedge_range;
#endif

#ifndef DOXYGEN_RUNNING
    typedef Index_iterator<Edge_index> Edge_iterator;
#endif

    /// \brief The range over all edge indices.
    ///
    /// A model of <a href="http://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a> with value type `Edge_index`.
    /// \sa `edges()`
    /// \sa `Halfedge_range`, `Vertex_range`, `Face_range`
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type Edge_range;
#else
    typedef Iterator_range<Edge_iterator> Edge_range;
#endif


#ifndef DOXYGEN_RUNNING
    typedef Index_iterator<Face_index> Face_iterator;
#endif
    /// \brief The range over all face indices.
    ///
    /// A model of <a href="http://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a> with value type `Face_index`.
    /// \sa `faces()`
    /// \sa `Vertex_range`, `Halfedge_range`, `Edge_range`
 #ifdef DOXYGEN_RUNNING
    typedef unspecified_type Face_range;
#else
   typedef Iterator_range<Face_iterator> Face_range;
#endif

#ifndef DOXYGEN_RUNNING 

  typedef CGAL::Vertex_around_target_iterator<Surface_mesh> Vertex_around_target_iterator;
  typedef Iterator_range<Vertex_around_target_iterator> Vertex_around_target_range;

  typedef CGAL::Halfedge_around_target_iterator<Surface_mesh>  Halfedge_around_target_iterator;
  typedef Iterator_range<Halfedge_around_target_iterator> Halfedge_around_target_range;

  typedef CGAL::Face_around_target_iterator<Surface_mesh>  Face_around_target_iterator;
  typedef Iterator_range<Face_around_target_iterator> Face_around_target_range;

  typedef CGAL::Vertex_around_face_iterator<Surface_mesh>  Vertex_around_face_iterator;
  typedef Iterator_range<Vertex_around_face_iterator> Vertex_around_face_range;

  typedef CGAL::Halfedge_around_face_iterator<Surface_mesh>  Halfedge_around_face_iterator;
  typedef Iterator_range<Halfedge_around_face_iterator> Halfedge_around_face_range;

  typedef CGAL::Face_around_face_iterator<Surface_mesh>  Face_around_face_iterator;
  typedef Iterator_range<Face_around_face_iterator> Face_around_face_range;
#endif

    /// @cond CGAL_BEGIN_END
    /// Start iterator for vertices.
    Vertex_iterator vertices_begin() const
    {
        return Vertex_iterator(Vertex_index(0), this);
    }

    /// End iterator for vertices.
    Vertex_iterator vertices_end() const
    {
        return Vertex_iterator(Vertex_index(num_vertices()), this);
    }
    /// @endcond


    /// returns the iterator range of the vertices of the mesh.
    Vertex_range vertices() const {
      return make_range(vertices_begin(), vertices_end());
    }


    /// @cond CGAL_BEGIN_END
    /// Start iterator for halfedges.
    Halfedge_iterator halfedges_begin() const
    {
        return Halfedge_iterator(Halfedge_index(0), this);
    }

    /// End iterator for halfedges.
    Halfedge_iterator halfedges_end() const
    {
        return Halfedge_iterator(Halfedge_index(num_halfedges()), this);
    }
    /// @endcond


    /// returns the iterator range of the halfedges of the mesh.
    Halfedge_range halfedges() const {
      return make_range(halfedges_begin(), halfedges_end());
    }


    /// @cond CGAL_BEGIN_END
    /// Start iterator for edges.
    Edge_iterator edges_begin() const
    {
        return Edge_iterator(Edge_index(0), this);
    }

    /// End iterator for edges.
    Edge_iterator edges_end() const
    {
        return Edge_iterator(Edge_index(num_edges()), this);
    }
    /// @endcond


    /// returns the iterator range of the edges of the mesh.
    Edge_range edges() const
    {
        return make_range(edges_begin(), edges_end());
    }


    /// @cond CGAL_BEGIN_END
    /// Start iterator for faces.
    Face_iterator faces_begin() const
    {
        return Face_iterator(Face_index(0), this);
    }

    /// End iterator for faces.
    Face_iterator faces_end() const
    {
        return Face_iterator(Face_index(num_faces()), this);
    }
    /// @endcond

    /// returns the iterator range of the faces of the mesh.
    Face_range faces() const {
      return make_range(faces_begin(), faces_end());
    }

#ifndef DOXYGEN_RUNNING
    /// returns the iterator range for vertices around vertex `target(h)`, starting at `source(h)`.
    Vertex_around_target_range vertices_around_target(Halfedge_index h) const
    {
      return CGAL::vertices_around_target(h,*this);
    }

    /// returns the iterator range for incoming halfedges around vertex `target(h)`, starting at `h`.
    Halfedge_around_target_range halfedges_around_target(Halfedge_index h) const
    {
      return CGAL::halfedges_around_target(h,*this);
    }

    /// returns the iterator range for faces around vertex `target(h)`, starting at `face(h)`.
    Face_around_target_range faces_around_target(Halfedge_index h) const
    {
      return CGAL::faces_around_target(h,*this);
    }

    /// returns the iterator range for vertices around face `face(h)`, starting at `target(h)`.
    Vertex_around_face_range vertices_around_face(Halfedge_index h) const
     {
       return CGAL::vertices_around_face(h,*this);
     }

    /// returns the iterator range for halfedges around face `face(h)`, starting at `h`.
    Halfedge_around_face_range halfedges_around_face(Halfedge_index h) const
    {
      return CGAL::halfedges_around_face(h,*this);
    }

    /// returns the iterator range for halfedges around face `face(h)`, starting at `h`.
    Face_around_face_range faces_around_face(Halfedge_index h) const
    {
       return CGAL::faces_around_face(h,*this);
    }

#endif

    ///@}


public: 

#ifndef DOXYGEN_RUNNING
    /// \name Circulator Types
    ///
    /// The following circulators enable to iterate through the elements around a face or vertex.
    /// As explained in the \ref SurfaceMeshOrientation "User Manual", we can speak of a  
    /// *clockwise* or *counterclockwise*
    /// traversal, by looking at the surface from the right side.  
    ///@{

    /// \brief This class circulates clockwise through all 
    /// one-ring neighbors of a vertex. 
    ///  A model of `BidirectionalCirculator` with value type `Vertex_index`.
    /// \sa `Halfedge_around_target_circulator`, `Face_around_target_circulator`

  typedef CGAL::Vertex_around_target_circulator<Surface_mesh> Vertex_around_target_circulator;



    /// \brief This class circulates clockwise through all incident faces of a vertex.
    ///  A model of `BidirectionalCirculator` with value type `Face_index`.
    /// \sa `Vertex_around_target_circulator`, `Halfedge_around_target_circulator`

  typedef CGAL::Face_around_target_circulator<Surface_mesh> Face_around_target_circulator;


    /// \brief This class circulates clockwise through all halfedges around a vertex that have this vertex as target.
    ///  A model of `BidirectionalCirculator` with value type `Halfedge_index`.
    /// \sa `Vertex_around_target_circulator`, `Halfedge_around_target_circulator`

  typedef CGAL::Halfedge_around_target_circulator<Surface_mesh> Halfedge_around_target_circulator;


    /// \brief This class circulates clockwise through all halfedges around a vertex that have this vertex as source.
    ///  A model of `BidirectionalCirculator` with value type `Halfedge_index`.
    /// \sa `Vertex_around_target_circulator`, `Halfedge_around_target_circulator`

  typedef CGAL::Halfedge_around_source_circulator<Surface_mesh> Halfedge_around_source_circulator;

    /// \brief This class circulates counterclockwise through all vertices around a face.
    ///  A model of `BidirectionalCirculator` with value type `Vertex_index`.

  typedef  CGAL::Vertex_around_face_circulator<Surface_mesh> Vertex_around_face_circulator;


    /// \brief This class circulates counterclockwise through all halfedges around a face.
    ///  A model of `BidirectionalCirculator` with value type `Halfedge_index`.

  typedef  CGAL::Halfedge_around_face_circulator<Surface_mesh> Halfedge_around_face_circulator;

   /// \brief This class circulates counterclockwise through all faces around a face.
   ///  A model of `BidirectionalCirculator` with value type `Face_index`.
   ///  Note that the face index is the same after `operator++`, if the neighboring faces share 
   ///  several halfedges.

  typedef  CGAL::Face_around_face_circulator<Surface_mesh> Face_around_face_circulator;
  /// @}
#endif

  /// @cond CGAL_DOCUMENT_INTERNALS
  // typedefs which make it easier to write the partial specialisation of boost::graph_traits

  typedef Vertex_index   vertex_index;
  typedef P                   vertex_property_type;
  typedef Halfedge_index halfedge_index;
  typedef Edge_index     edge_index;
  typedef Face_index     face_index;

  typedef Vertex_iterator     vertex_iterator;
  typedef Halfedge_iterator   halfedge_iterator;
  typedef Edge_iterator       edge_iterator;
  typedef Face_iterator      face_iterator;
  typedef CGAL::Out_edge_iterator<Self>     out_edge_iterator;

  typedef boost::undirected_tag             directed_category;
  typedef boost::disallow_parallel_edge_tag edge_parallel_category; 

  struct traversal_category : public virtual boost::bidirectional_graph_tag,
                              public virtual boost::vertex_list_graph_tag,
                              public virtual boost::edge_list_graph_tag
  {};

  typedef size_type vertices_size_type;
  typedef size_type halfedges_size_type;
  typedef size_type edges_size_type;
  typedef size_type faces_size_type;
  typedef size_type degree_size_type;

 /// @endcond
public:

    /// \name Construction, Destruction, Assignment
    ///
    ///  Copy constructors as well as assignment do also copy simplices marked as removed.
    ///@{

    /// %Default constructor.
    Surface_mesh();

    /// Copy constructor: copies `rhs` to `*this`. Performs a deep copy of all properties.
    Surface_mesh(const Surface_mesh& rhs) { *this = rhs; }

    /// assigns `rhs` to `*this`. Performs a deep copy of all properties.
    Surface_mesh& operator=(const Surface_mesh& rhs);

    /// assigns `rhs` to `*this`. Does not copy custom properties.
    Surface_mesh& assign(const Surface_mesh& rhs);

    ///@}

public:

    /// \name Adding Vertices, Edges, and Faces
    ///@{

   /// adds a new vertex, and resizes vertex properties if necessary.
    Vertex_index add_vertex()
    {
      size_type inf = -1;
      if(vertices_freelist_ != inf){
        size_type idx = vertices_freelist_;
        vertices_freelist_ = (size_type)vconn_[Vertex_index(vertices_freelist_)].halfedge_;
        --removed_vertices_;
        vremoved_[Vertex_index(idx)] = false;
        return Vertex_index(idx);
      } else {
        vprops_.push_back();
        return Vertex_index(num_vertices()-1);
      }
    }

    /// adds a new vertex, resizes vertex properties if necessary,
    /// and sets the \em point property to `p`.
    /// \note Several vertices may have the same point property.
    Vertex_index add_vertex(const Point& p) 
    {
        Vertex_index v = add_vertex();
        vpoint_[v] = p;
        return v;
    }



public:

    /// adds a new edge, and resizes edge and halfedge properties if necessary.
    Halfedge_index add_edge()
    {
      Halfedge_index h0, h1;
      size_type inf = -1;
      if(edges_freelist_ != inf){
        size_type idx = edges_freelist_;
        edges_freelist_ = (size_type)hconn_[Halfedge_index(edges_freelist_)].next_halfedge_;
        --removed_edges_;
        eremoved_[Edge_index(Halfedge_index(idx))] = false;
        return Halfedge_index(idx);
      } else {
        eprops_.push_back();
        hprops_.push_back();
        hprops_.push_back();

        return Halfedge_index(num_halfedges()-2);
      }
    }

    /// adds two opposite halfedges, and resizes edge and halfedge properties if necessary.
    /// Sets the targets of the halfedge to the given vertices, but does not modify the halfedge
    /// associated to the vertices.
    /// \note The function does not check whether there is already an edge between the vertices.
    /// \returns the halfedge with `v1` as target

    Halfedge_index add_edge(Vertex_index v0, Vertex_index v1)
    {
        CGAL_assertion(v0 != v1);
        Halfedge_index h = add_edge();

        set_target(h, v1);
        set_target(opposite(h), v0);

        return h;
    }

    /// adds a new face, and resizes face properties if necessary.
    Face_index add_face()
    {
      size_type inf = -1;
      if(faces_freelist_ != inf){
        size_type idx = faces_freelist_;
        faces_freelist_ = (size_type)fconn_[Face_index(faces_freelist_)].halfedge_;
        --removed_faces_;
        fremoved_[Face_index(idx)] = false;
        return Face_index(idx);
      } else {
        fprops_.push_back();
        return Face_index(num_faces()-1);
      }
    }

    /// if possible, adds a new face with vertices from a range with value type `Vertex_index`.
    /// The function adds halfedges between successive vertices if they are not yet indicent to halfedges,
    /// or updates the connectivity of halfedges already in place. 
    /// Resizes halfedge, edge, and face properties if necessary.
    /// \returns the face index of the added face, or `Surface_mesh::null_face()` if the face could not be added.
    template <typename Range>
    Face_index add_face(const Range& vertices);


    /// adds a new triangle connecting vertices `v0`, `v1`, `v2`.
    /// \returns the face index of the added face, or `Surface_mesh::null_face()` if the face could not be added.
    Face_index add_face(Vertex_index v0, Vertex_index v1, Vertex_index v2)
    {
        boost::array<Vertex_index, 3> 
            v = {{v0, v1, v2}};
        return add_face(v);
    }

    /// adds a new quad connecting vertices `v0`, `v1`, `v2`, `v3`.
    /// \returns the face index of the added face, or `Surface_mesh::null_face()` if the face could not be added.
    Face_index add_face(Vertex_index v0, Vertex_index v1, Vertex_index v2, Vertex_index v3)
    {
        boost::array<Vertex_index, 4> 
            v = {{v0, v1, v2, v3}};
        return add_face(v);
    }

    ///@}


 
    /// \name Low-Level Removal Functions 
    ///
    /// Although the elements are only marked as removed
    /// their connectivity and properties should not be used.
    ///
    /// \warning Functions in this group do not adjust any of
    /// connected elements and usually leave the surface mesh in an
    /// invalid state.
    /// 
    ///
    /// @{

    /// removes vertex `v` from the halfedge data structure without
    /// adjusting anything.
    void remove_vertex(Vertex_index v)
    {
        vremoved_ = add_property_map<Vertex_index, bool>("v:removed", false).first;
        vremoved_[v] = true; ++removed_vertices_; garbage_ = true;
        vconn_[v].halfedge_ = Halfedge_index(vertices_freelist_);
        vertices_freelist_ = (size_type)v;
    }

    /// removes the two halfedges corresponding to `e` from the halfedge data structure without
    /// adjusting anything.
    void remove_edge(Edge_index e)
    {
        eremoved_ = add_property_map<Edge_index, bool>("e:removed", false).first;
        eremoved_[e] = true; ++removed_edges_; garbage_ = true;
        hconn_[Halfedge_index((size_type)e << 1)].next_halfedge_ = Halfedge_index(edges_freelist_ );
        edges_freelist_ = ((size_type)e << 1);
    }

    /// removes  face `f` from the halfedge data structure without
    /// adjusting anything.

    void remove_face(Face_index f)
    {
        fremoved_ = add_property_map<Face_index, bool>("f:removed", false).first;
        fremoved_[f] = true; ++removed_faces_; garbage_ = true;
        fconn_[f].halfedge_ = Halfedge_index(faces_freelist_);
        faces_freelist_ = (size_type)f;
    }


    ///@}


    /// \name Memory Management
    ///
    /// Functions to check the number of elements, the amount of space
    /// allocated for elements, and to clear the structure.
    ///@{

  /// returns the number of vertices in the mesh.
  size_type number_of_vertices() const
  {
    return num_vertices() - number_of_removed_vertices();
  }
 
  /// returns the number of halfedges in the mesh.
  size_type number_of_halfedges() const
  {
    return num_halfedges() - number_of_removed_halfedges();
  }

  /// returns the number of edges in the mesh.
  size_type number_of_edges() const
  {
    return num_edges() - number_of_removed_edges();
  }

  /// returns the number of faces in the mesh.
  size_type number_of_faces() const
  {
    return num_faces() - number_of_removed_faces();
  }

    /// returns `true` iff the mesh is empty, i.e., has no vertices, halfedges and faces.
    bool is_empty() const
  {
    return ( num_vertices() == number_of_removed_vertices()
             && num_halfedges() == number_of_removed_halfedges()
             && num_faces() == number_of_removed_faces());
  }

    /// removes all vertices, halfedge, edges and faces. Collects garbage and clears all properties.
    void clear();

 
    /// reserves space for vertices, halfedges, edges, faces, and their currently
    /// associated properties.
    void reserve(size_type nvertices,
                 size_type nedges,
                 size_type nfaces )
    {
        vprops_.reserve(nvertices);
        hprops_.reserve(2*nedges);
        eprops_.reserve(nedges);
        fprops_.reserve(nfaces);
    }

      void resize(size_type nvertices,
                 size_type nedges,
                 size_type nfaces )
    {
        vprops_.resize(nvertices);
        hprops_.resize(2*nedges);
        eprops_.resize(nedges);
        fprops_.resize(nfaces);
    }
  
  bool join(const Surface_mesh& other)
  {
    size_type nv = num_vertices(), nh = num_halfedges(), nf = num_faces();
    resize(num_vertices()+  other.num_vertices(),
            num_edges()+  other.num_edges(),
            num_faces()+  other.num_faces());

    vprops_.transfer(other.vprops_);
    hprops_.transfer(other.hprops_);
    fprops_.transfer(other.fprops_);
    eprops_.transfer(other.eprops_);

    for(size_type i = nv; i < nv+other.num_vertices(); i++){
      Vertex_index vi(i);
      if(vconn_[vi].halfedge_ != null_halfedge()){
        vconn_[vi].halfedge_ = Halfedge_index(size_type(vconn_[vi].halfedge_)+nh);
      }
    }
    for(size_type i = nf; i < nf+other.num_faces(); i++){
      Face_index fi(i);
      if(fconn_[fi].halfedge_ != null_halfedge()){
        fconn_[fi].halfedge_ = Halfedge_index(size_type(fconn_[fi].halfedge_)+nh);
      }
    }
    for(size_type i = nh; i < nh+other.num_halfedges(); i++){
      Halfedge_index hi(i);
      if(hconn_[hi].face_ != null_face()){
        hconn_[hi].face_ = Face_index(size_type(hconn_[hi].face_)+nf);
      }
      if( hconn_[hi].vertex_ != null_vertex()){
        hconn_[hi].vertex_ = Vertex_index(size_type(hconn_[hi].vertex_)+nv);
      }
      if(hconn_[hi].next_halfedge_ != null_halfedge()){
        hconn_[hi].next_halfedge_ = Halfedge_index(size_type(hconn_[hi].next_halfedge_)+nh);
      }
      if(hconn_[hi].prev_halfedge_ != null_halfedge()){
        hconn_[hi].prev_halfedge_ = Halfedge_index(size_type(hconn_[hi].prev_halfedge_)+nh);
      }
    }
    size_type nil = -1;
    if(other.vertices_freelist_ != nil){
      if(vertices_freelist_ != nil){
        Vertex_index vi(nv+other.vertices_freelist_);
        Halfedge_index inf(-1);
        while(vconn_[vi].halfedge_ != inf){
          vi = Vertex_index(size_type(vconn_[vi].halfedge_));
        }
        vconn_[vi].halfedge_ = Halfedge_index(vertices_freelist_);
      }
      vertices_freelist_ = nv + other.vertices_freelist_; 
    }
    if(other.faces_freelist_ != nil){
      if(faces_freelist_ != nil){
        Face_index fi(nf+other.faces_freelist_);
        Halfedge_index inf(-1);
        while(fconn_[fi].halfedge_ != inf){
          fi = Face_index(size_type(fconn_[fi].halfedge_));
        }
        fconn_[fi].halfedge_ = Halfedge_index(faces_freelist_);
      }
      faces_freelist_ = nf + other.faces_freelist_; 
    }
    if(other.edges_freelist_ != nil){
      if(edges_freelist_ != nil){
        Halfedge_index hi((nh>>1)+other.edges_freelist_);
        Halfedge_index inf(-1);
        while(hconn_[hi].next_halfedge_ != inf){
          hi = hconn_[hi].next_halfedge_;
        }
        hconn_[hi].next_halfedge_ = Halfedge_index(edges_freelist_);
      }
      edges_freelist_ = (nh>>1) + other.edges_freelist_; 
    }
    garbage_ = garbage_ || other.garbage_;
    removed_vertices_ += other.removed_vertices_;
    removed_edges_ += other.removed_edges_;
    removed_faces_ += other.removed_faces_;
    return true;
  }

    ///@}

    
    /// \name Garbage Collection
    ///
    /// While removing elements only marks them as removed
    /// garbage collection really removes them.
    /// The API in this section allows to check whether 
    /// an element is removed, to get the number of
    /// removed elements, and to collect garbage.
    /// The number of elements together with the number of  removed elements is
    /// an upperbound on the index, and is needed
    /// by algorithms that temporarily store a 
    /// property in a vector of the appropriate size.
    /// Note however that by garbage collecting elements get new indices.
    /// In case you store indices in an auxiliary data structure
    /// or in a property these indices are potentially no longer 
    /// refering to the right elements. 


    ///@{
#ifndef DOXYGEN_RUNNING
   /// returns the number of used and removed vertices in the mesh.
    size_type num_vertices() const { return (size_type) vprops_.size(); }

    /// returns the number of used and removed halfedges in the mesh.
    size_type num_halfedges() const { return (size_type) hprops_.size(); }

    /// returns the number of used and removed edges in the mesh.
    size_type num_edges() const { return (size_type) eprops_.size(); }

    /// returns the number of used and removed faces in the mesh.
    size_type num_faces() const { return (size_type) fprops_.size(); }

#endif

    /// returns the number of vertices in the mesh which are marked removed.
    size_type number_of_removed_vertices() const { return removed_vertices_; }

    /// returns the number of halfedges in the mesh which are marked removed.
    size_type number_of_removed_halfedges() const { return 2*removed_edges_; }

    /// returns the number of edges in the mesh which are marked removed.
    size_type number_of_removed_edges() const { return removed_edges_; }

    /// returns the number offaces in the mesh which are marked removed.
    size_type number_of_removed_faces() const { return removed_faces_; }



    /// returns whether vertex `v` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Vertex_index v) const
    {
        return vremoved_[v];
    }
    /// returns whether halfedge `h` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Halfedge_index h) const
    {
        return eremoved_[edge(h)];
    }
    /// returns whether edge `e` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Edge_index e) const
    {
        return eremoved_[e];
    }
    /// returns whether face `f` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Face_index f) const
    {
        return fremoved_[f];
    }

    /// checks if any vertices, halfedges, edges, or faces are marked as removed.
    /// \sa collect_garbage
    bool has_garbage() const { return garbage_; }

    /// really removes vertices, halfedges, edges, and faces which are marked removed.
    /// \sa `has_garbage()`
    /// \attention By garbage collecting elements get new indices.
    /// In case you store indices in an auxiliary data structure
    /// or in a property these indices are potentially no longer 
    /// refering to the right elements. 
    void collect_garbage();


    /// @cond CGAL_DOCUMENT_INTERNALS
    /// removes unused memory from vectors. This shrinks the storage
    /// of all properties to the minimal required size.
    /// \attention Invalidates all existing references to properties.

    void shrink_to_fit()
    {
        vprops_.shrink_to_fit();
        hprops_.shrink_to_fit();
        eprops_.shrink_to_fit();
        fprops_.shrink_to_fit();
    }
    /// @endcond

    ///@}

    /// @cond CGAL_DOCUMENT_INTERNALS
    ///
    /// \name Simple Validity Checks
    ///
    /// Functions in this group check if the index is valid, that is between
    /// `0` and the currently allocated maximum amount of the
    /// elements. They do not check if an element is marked as removed.
    ///@{

    /// returns whether the index of vertex `v` is valid, that is within the current array bounds.
    bool has_valid_index(Vertex_index v) const
    {
      return ((size_type)v < num_vertices());
    }
    /// returns whether the index of halfedge `h` is valid, that is within the current array bounds.
    bool has_valid_index(Halfedge_index h) const
    {
      return ((size_type)h < num_halfedges());
    }
    /// returns whether the index of edge `e` is valid, that is within the current array bounds.
    bool has_valid_index(Edge_index e) const
    {
      return ((size_type)e < num_edges());
    }
    /// returns whether the index of face `f` is valid, that is within the current array bounds.
    bool has_valid_index(Face_index f) const
    {
        return ((size_type)f < num_faces());
    }

    /// @}
    /// @endcond
    
    /// \name Validity Checks
    ///
    /// Functions in this group perform checks for structural
    /// consistency of a complete surface mesh, or an individual element. 
    /// They are expensive and should only be used in debug configurations.

    ///@{

    /// perform an expensive validity check on the data structure and
    /// print found errors to `std::cerr` when `verbose == true`.
  bool is_valid(bool verbose = true) const
    {
        bool valid = true;
        size_type vcount = 0, hcount = 0, fcount = 0;
        for(Halfedge_iterator it = halfedges_begin(); it != halfedges_end(); ++it) { 
            ++hcount;
            valid = valid && next(*it).is_valid();
            valid = valid && opposite(*it).is_valid();
            if(!valid) {
                if (verbose)
                  std::cerr << "Integrity of halfedge " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && (opposite(*it) != *it);
            valid = valid && (opposite(opposite(*it)) == *it);
            if(!valid) {
              if (verbose)
                std::cerr << "Integrity of opposite halfedge of " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && (next(prev(*it)) == *it);
            if(!valid) {
                if (verbose)
                  std::cerr << "Integrity of previous halfedge of " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && (prev(next(*it)) == *it);
            if(!valid) {
                if (verbose)
                  std::cerr << "Integrity of next halfedge of " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && target(*it).is_valid();
            if(!valid) {
                if (verbose)
                  std::cerr << "Integrity of vertex of halfedge " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && (target(*it) == target(opposite(next(*it))));
            if(!valid) {
                if (verbose)
                  std::cerr << "Halfedge vertex of next opposite is not the same for " << *it << "."  << std::endl;
                break;
            }
        }

        for(Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it) {
          ++vcount;
            if(halfedge(*it).is_valid()) {
                // not an isolated vertex
                valid = valid && (target(halfedge(*it)) == *it);
                if(!valid) {
                    if (verbose)
                      std::cerr << "Halfedge of " << *it << " is not an incoming halfedge." << std::endl;
                    break;
                }
            }
        }
        for(Face_iterator it = faces_begin(); it != faces_end(); ++it) {
          ++fcount;
        }
        
        valid = valid && (vcount == number_of_vertices());
        if(!valid && verbose){
          std::cerr << "#vertices: iterated: " << vcount << " vs number_of_vertices(): " << number_of_vertices()<< std::endl;
        }

        valid = valid && (hcount == number_of_halfedges());
        if(!valid && verbose){
          std::cerr << "#halfedges: iterated: " << hcount << " vs number_of_halfedges(): " << number_of_halfedges()<< std::endl;
        }

        valid = valid && (fcount == number_of_faces());
        if(!valid && verbose){
          std::cerr << "#faces: iterated: " << fcount << " vs number_of_faces(): " << number_of_faces()<< std::endl;
        }
        return valid;
    }

    /// performs a validity check on a single vertex.
    bool is_valid(Vertex_index v) const {
        Halfedge_index h = vconn_[v].halfedge_;
        if(h!= null_halfedge() && (!has_valid_index(h) || is_removed(h))) {
          std::cerr << "Vertex connectivity halfedge error in " << (size_type)v
                    << " with " << (size_type)h << std::endl;
            return false;
        }
        return true;
    }

    /// performs a validity check on a single halfedge.
    bool is_valid(Halfedge_index h) const {
        Face_index f = hconn_[h].face_;
        Vertex_index v = hconn_[h].vertex_;
        Halfedge_index hn = hconn_[h].next_halfedge_;
        Halfedge_index hp = hconn_[h].prev_halfedge_;

        bool valid = true;
        // don't validate the face if this is a border halfedge
        if(!is_border(h)) {
            if(!has_valid_index(f) || is_removed(f)) {
                std::cerr << "Halfedge connectivity Face "
                          << (!has_valid_index(f) ? "invalid" : "removed")
                          << " in " << (size_type)h << std::endl;
                valid = false;
            }
        }

        if(!has_valid_index(v) || is_removed(v)) {
            std::cerr << "Halfedge connectivity Vertex "
                      << (!has_valid_index(v) ? "invalid" : "removed")
                      << " in " << (size_type)h << std::endl;
            valid = false;
        }

        if(!has_valid_index(hn) || is_removed(hn)) {
            std::cerr << "Halfedge connectivity hnext "
                      << (!has_valid_index(hn) ? "invalid" : "removed")
                      << " in " << (size_type)h << std::endl;
            valid = false;
        }
        if(!has_valid_index(hp) || is_removed(hp)) {
            std::cerr << "Halfedge connectivity hprev "
                      << (!has_valid_index(hp) ? "invalid" : "removed")
                      << " in " << (size_type)h << std::endl;
            valid = false;
        }
        return valid;
    }


    /// performs a validity check on a single ede.
    bool is_valid(Edge_index e) const {
      Halfedge_index h = halfedge(e);
      return is_valid(h) && is_valid(opposite(h));
    }


    /// performs a validity check on a single face.
    bool is_valid(Face_index f) const {
        Halfedge_index h = fconn_[f].halfedge_;
        if(!has_valid_index(h) || is_removed(h)) {
          std::cerr << "Face connectivity halfedge error in " << (size_type)f
                      << " with " << (size_type)h << std::endl;
            return false;
        }
        return true;
    }

    ///@}



    /// \name Low-Level Connectivity
    ///@{

    /// returns the vertex the halfedge `h` points to.
    Vertex_index target(Halfedge_index h) const
    {
        return hconn_[h].vertex_;
    }

    /// sets the vertex the halfedge `h` points to to `v`.
    void set_target(Halfedge_index h, Vertex_index v)
    {
        hconn_[h].vertex_ = v;
    }

    /// returns the face incident to halfedge `h`.
    Face_index face(Halfedge_index h) const
    {
        return hconn_[h].face_;
    }

    /// sets the incident face to halfedge `h` to `f`.
    void set_face(Halfedge_index h, Face_index f)
    {
        hconn_[h].face_ = f;
    }

    /// returns the next halfedge within the incident face.
    Halfedge_index next(Halfedge_index h) const
    {
        return hconn_[h].next_halfedge_;
    }

    /// returns the previous halfedge within the incident face.
    Halfedge_index prev(Halfedge_index h) const
    {
        return hconn_[h].prev_halfedge_;
    }

    /// @cond CGAL_DOCUMENT_INTERNALS
    // sets the next halfedge of `h` within the face to `nh`.
    void set_next_only(Halfedge_index h, Halfedge_index nh)
    {
      hconn_[h].next_halfedge_ = nh;
    }

    // sets previous halfedge of `h` to `nh`.
    void set_prev_only(Halfedge_index h, Halfedge_index nh)
    {
      if(h != null_halfedge()){
        hconn_[h].prev_halfedge_ = nh;
      }
    }
    /// @endcond  

    /// sets the next halfedge of `h` within the face to `nh` and
    /// the previous halfedge of `nh` to `h`.
    void set_next(Halfedge_index h, Halfedge_index nh)
    {
      set_next_only(h, nh);
      set_prev_only(nh, h);
    }

    /// returns an incoming halfedge of vertex `v`.
    /// If `v` is a border vertex this will be a border halfedge.
    /// \invariant `target(halfedge(v)) == v`
    Halfedge_index halfedge(Vertex_index v) const
    {
        return vconn_[v].halfedge_;
    }

    /// sets the incoming halfedge of vertex `v` to `h`.
    void set_halfedge(Vertex_index v, Halfedge_index h)
    {
        vconn_[v].halfedge_ = h;
    }


    /// returns a halfedge of face `f`.
    Halfedge_index halfedge(Face_index f) const
    {
        return fconn_[f].halfedge_;
    }

    /// sets the halfedge of face `f` to `h`.
    void set_halfedge(Face_index f, Halfedge_index h)
    {
        fconn_[f].halfedge_ = h;
    }

    /// returns the opposite halfedge of `h`. Note that there is no function `set_opposite()`.
    Halfedge_index opposite(Halfedge_index h) const
    {
        return Halfedge_index(((size_type)h & 1) ? (size_type)h-1 : (size_type)h+1);
    }

    ///@}

    /// \name Low-Level Connectivity Convenience Functions
    ///@{

    /// returns the vertex the halfedge `h` emanates from.
    Vertex_index source(Halfedge_index h) const
    {
        return target(opposite(h));
    }

    /// returns `opposite(next(h))`, that is the next halfedge \ref SurfaceMeshOrientation 
    /// "clockwise" around the target vertex of `h`.
    Halfedge_index next_around_target(Halfedge_index h) const
    {
        return opposite(next(h));
    }

    /// returns `prev(opposite(h))`, that is the previous halfedge \ref SurfaceMeshOrientation
    /// "clockwise" around the target vertex of `h`.
    Halfedge_index prev_around_target(Halfedge_index h) const
    {
        return prev(opposite(h));
    }

    /// returns `next(opposite(h))`, that is the next halfedge \ref SurfaceMeshOrientation 
    /// "clockwise" around the source vertex of `h`.
    Halfedge_index next_around_source(Halfedge_index h) const
    {
        return next(opposite(h));
    }

    /// returns `opposite(prev(h))`, that is the previous halfedge \ref SurfaceMeshOrientation
    /// "clockwise" around the source vertex of `h`.
    Halfedge_index prev_around_source(Halfedge_index h) const
    {
        return opposite(prev(h));
    }

    /// returns the i'th vertex of edge `e`, for `i=0` or `1`.
    Vertex_index vertex(Edge_index e, unsigned int i) const
    {
        CGAL_assertion(i<=1);
        return target(halfedge(e, i));
    }

    /// finds a halfedge between two vertices. Returns a default constructed
    /// `Halfedge_index`, if  `source` and  `target` are not connected.
    Halfedge_index halfedge(Vertex_index source, Vertex_index target) const;

    ///@}


    /// \name Switching between Halfedges and Edges
    ///@{

    /// returns the edge that contains halfedge `h` as one of its two halfedges.
    Edge_index edge(Halfedge_index h) const
    {
        return Edge_index(h);
    }

    /// returns the halfedge corresponding to the edge `e`.
    Halfedge_index halfedge(Edge_index e) const
    {
        return Halfedge_index(e.halfedge());
    }

    /// returns the i'th halfedge of edge `e`, for `i=0` or `1`.
    Halfedge_index halfedge(Edge_index e, unsigned int i) const
    {
        CGAL_assertion(i<=1);
        return Halfedge_index(((size_type)e << 1) + i);
    }

    ///@}


    /// \name Degree Functions
    ///@{

    /// returns the number of incident halfedges of vertex `v`.
    size_type degree(Vertex_index v) const;

    /// returns the number of incident halfedges of face `f`.
    size_type degree(Face_index f) const;

    ///@}



    /// \name Borders
    ///
    ///  A halfedge, or edge is on the border of a surface mesh
    /// if it is incident to a `null_face()`.  A vertex is on a border
    /// if it is incident to a border halfedge. While for a halfedge and
    /// edge this is a constant time operation, for a vertex it means
    /// to look at all incident halfedges.  If algorithms operating on a 
    /// surface mesh maintain that the halfedge associated to a border vertex is
    /// a border halfedge, this is a constant time operation too.  
    /// This section provides functions to check if an element is on a 
    /// border and to change the halfedge associated to a border vertex.
    ///@{

    /// returns whether `v` is a border vertex. 
    /// \cgalAdvancedBegin
    /// With the default value for
    /// `check_all_incident_halfedges` the function iteratates over the incident halfedges.
    /// With `check_all_incident_halfedges == false` the function returns `true`, if the incident
    /// halfedge associated to vertex `v` is a border halfedge.
    /// \cgalAdvancedEnd
  bool is_border(Vertex_index v, bool check_all_incident_halfedges = true) const
    {
        Halfedge_index h(halfedge(v));
        if (h == null_halfedge()){
          return false;
        }
        if(check_all_incident_halfedges){
          Halfedge_around_target_circulator hatc(h,*this), done(hatc);
          do {
            if(is_border(*hatc)){
              return true;
            }
          }while(++hatc != done);
          return false;
        }
        return (!(is_valid(h) && is_border(h)));
    }

    /// returns whether `h` is a border halfege, that is if its incident face is `sm.null_face()`.
    bool is_border(Halfedge_index h) const
    {
        return !face(h).is_valid();
    }


    /// returns whether `e` is a border edge, i.e., if any 
    /// of its two halfedges is a border halfedge.
    bool is_border(Edge_index e) const
    {
      return is_border(e.halfedge()) || is_border(opposite(e.halfedge()));
    }

  /// iterates over the incident halfedges and sets the incident halfedge
  /// associated to vertex `v` to a border halfedge and returns `true` if it exists.
  bool set_vertex_halfedge_to_border_halfedge(Vertex_index v)
  {
    if(halfedge(v) == null_halfedge()){
      return false;
    }
    Halfedge_around_target_circulator hatc(halfedge(v),*this), done(hatc);
    do {
      if(is_border(*hatc)){
        set_halfedge(v,*hatc);
        return true;
      }
    }while(++hatc != done);
    return false;
  }

  /// applies `set_vertex_halfedge_to_border_halfedge(Vertex_index)` on all vertices 
  /// around the face associated to `h`.
  void set_vertex_halfedge_to_border_halfedge(Halfedge_index h)
  {
    if(is_border(h)){
      Halfedge_around_face_circulator hafc(h,*this),done(hafc);
      do {
        set_halfedge(target(*hafc),*hafc);
      }while(++hafc != done);
    } else {
       Vertex_around_face_circulator vafc(h,*this),done(vafc);
      do {
        set_vertex_halfedge_to_border_halfedge(*vafc);
      }while(++vafc != done);
    }
  }

  /// applies `set_vertex_halfedge_to_border_halfedge(Vertex_index)` on all vertices 
  /// of the surface mesh.
  void set_vertex_halfedge_to_border_halfedge()
  {
    BOOST_FOREACH(Halfedge_index h, halfedges()){
      if(is_border(h)){
          set_halfedge(target(h),h);
        }
    }
  }


    /// returns whether `v` is isolated, i.e., incident to `Surface_mesh::null_halfedge()`.
    bool is_isolated(Vertex_index v) const
    {
        return !halfedge(v).is_valid();
    }

    ///@}


private: //--------------------------------------------------- property handling

    /// @cond BROKEN_DOC
    typedef boost::fusion::map<
      boost::fusion::pair< typename Surface_mesh::Vertex_index, Property_container<Vertex_index> (Surface_mesh::*)>,
      boost::fusion::pair< typename Surface_mesh::Halfedge_index, Property_container<Halfedge_index> (Surface_mesh::*)>,
      boost::fusion::pair< typename Surface_mesh::Edge_index, Property_container<Edge_index> (Surface_mesh::*)>,
      boost::fusion::pair< typename Surface_mesh::Face_index, Property_container<Face_index> (Surface_mesh::*)> >
        map_type;

    map_type pmap_;
    /// @endcond

    public:
 

 /*! \name Property Handling

 A `Property_map<I,T>` allows to associate properties of type `T` to a vertex, halfdge, edge, or face index type I.
 Properties can be added, and looked up with a string, and they can be removed at runtime.
 The \em point property of type `P` is associated to the string "v:point". 

    */
    ///@{

  /// Model of `LvaluePropertyMap` with `I` as key type and `T` as value type, where `I`
  /// is either a vertex, halfedge, edge, or face index type.
#ifdef DOXYGEN_RUNNING
  template <class I, class T>
  using Property_map = unspecified_type;

#else


#endif

    /// adds a property map named `name` with value type `T` and default `t`
    /// for index type `I`. Returns the property map together with a Boolean 
    /// that is `true` if a new map was created. In case it already exists
    /// the existing map together with `false` is returned.

  
    template<class I, class T>
      std::pair<Property_map<I, T>, bool>
    add_property_map(const std::string& name, const T t=T()) {
        return (this->*boost::fusion::at_key<I>(pmap_)).template add<T>(name, t);
    }

 
    /// returns a property map named `name` with key type `I` and value type `T`, 
    /// and a Boolean that is `true` if the property exists. 
    /// In case it does not exist the Boolean is `false` and the behavior of
    /// the property map is undefined.
    template <class I, class T>
    std::pair<Property_map<I, T>,bool> property_map(const std::string& name) const
    {
      return (this->*boost::fusion::at_key<I>(pmap_)).template get<T>(name);
    }


    /// removes property map `p`. The memory allocated for that property map is
    /// freed.
    template<class I, class T>
    void remove_property_map(Property_map<I, T>& p)
    {
        (this->*boost::fusion::at_key<I>(pmap_)).remove(p);
    }

    /// @cond CGAL_DOCUMENT_INTERNALS
    /// returns the std::type_info of the value type of the
    /// property identified by `name`.  `typeid(void)` if `name`
    /// does not identify any property.
    ///
    /// @tparam I The key type of the property. 

    template<class I>
    const std::type_info& property_type(const std::string& name)
    {
        return (this->*boost::fusion::at_key<I>(pmap_)).get_type(name);
    }
  /// @endcond

    /// returns a vector with all strings that describe properties with the key type `I`.
    /// @tparam I The key type of the properties.
    template<class I>
    std::vector<std::string> properties() const
    {
        return (this->*boost::fusion::at_key<I>(pmap_)).properties();
    }

    /// returns the property for the string "v:point".
    Property_map<Vertex_index, Point>
    points() const { return vpoint_; }

    /// returns the point associated to vertex `v`.
    const Point&
    point(Vertex_index v) const { return vpoint_[v]; }

    /// returns the point associated to vertex `v`.
    Point&
    point(Vertex_index v) { return vpoint_[v]; }

    /// @cond CGAL_DOCUMENT_INTERNALS
    /// prints property statistics to the stream `out`. The output is human-readable but
    /// not machine-friendly.  
    ///
    void property_stats(std::ostream& out = std::cout) const;
    /// @endcond
    ///@}


 /// \name Null Elements
    ///@{

  /// returns `Vertex_index(-1)`.
  static Vertex_index null_vertex()
  {
    return vertex_index(-1);
  }

  /// returns `Edge_index(-1)`.
  static Edge_index null_edge()
  {
    return edge_index(-1);
  }
  /// returns `Halfedge_index(-1)`.
  static Halfedge_index null_halfedge()
  {
    return halfedge_index(-1);
  }
  /// returns `Face_index(-1)`.
  static Face_index null_face()
  {
    return face_index(-1);
  }
  /// @}

  
private: //--------------------------------------------------- helper functions


    /// make sure that the incoming halfedge of vertex v is a border halfedge
    /// if `v` is a border vertex.
    void adjust_incoming_halfedge(Vertex_index v);

private: //------------------------------------------------------- private data
    Property_container<Vertex_index> vprops_;
    Property_container<Halfedge_index> hprops_;
    Property_container<Edge_index> eprops_;
    Property_container<Face_index> fprops_;

    Property_map<Vertex_index, Vertex_connectivity>      vconn_;
    Property_map<Halfedge_index, Halfedge_connectivity>  hconn_;
    Property_map<Face_index, Face_connectivity>          fconn_;

    Property_map<Vertex_index, bool>  vremoved_;
    Property_map<Edge_index, bool>    eremoved_;
    Property_map<Face_index, bool>    fremoved_;

    Property_map<Vertex_index, Point>   vpoint_;

    size_type removed_vertices_;
    size_type removed_edges_;
    size_type removed_faces_;

    size_type vertices_freelist_;
    size_type edges_freelist_;
    size_type faces_freelist_;
    bool garbage_;
};

  /*! \addtogroup PkgSurface_mesh
   *
   * @{
   */

  /// \relates Surface_mesh
  /// Inserts `other` into `sm`. 
  /// Shifts the indices of vertices of `other` by `sm.number_of_vertices() + sm.number_of_removed_vertices()`
  /// and analoguously for halfedges, edges, and faces.
  /// Copies entries of all property maps which have the same name in `sm` and `other`. 
  /// that is, property maps which are only in `other` are ignored.
  /// Also copies elements which are marked as removed, and concatenates the freelists of `sm` and `other`. 

  template <typename P>
  Surface_mesh<P>& operator+=(Surface_mesh<P>& sm, const Surface_mesh<P>& other)
  {
    sm.join(other);
    return sm;
  }

  /// \relates Surface_mesh
  /// Inserts the surface mesh in an output stream in Ascii OFF format. 
  /// Only the \em point property is inserted in the stream.
  /// \pre `operator<<(std::ostream&,const P&)` must be defined.
  template <typename P>
  std::ostream& operator<<(std::ostream& os, const Surface_mesh<P>& sm)
  {
    typedef Surface_mesh<P> Mesh;
    typedef typename Mesh::Vertex_index Vertex_index;
    typedef typename Mesh::Face_index Face_index;

    os << "OFF\n" << sm.number_of_vertices() << " " << sm.number_of_faces() << " 0\n";
    std::vector<int> reindex;
    reindex.resize(sm.num_vertices());
    int n = 0;
    BOOST_FOREACH(Vertex_index v, sm.vertices()){
      os << sm.point(v) << std::endl;
      reindex[v]=n++;
    }

    BOOST_FOREACH(Face_index f, sm.faces()){
      os << sm.degree(f);
      BOOST_FOREACH(Vertex_index v, CGAL::vertices_around_face(sm.halfedge(f),sm)){
        os << " " << reindex[v];
      }
      os << "\n";
    }
    return os;
  }

/// @cond CGAL_DOCUMENT_INTERNALS
  inline std::istream& sm_skip_comments( std::istream& in) {
      char c;
      in >> c;
      if (c == '#')
        in.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
      else
        in.putback(c);
      return in;
  }
/// @endcond


  /// \relates Surface_mesh
  /// Extracts the surface mesh from an input stream in Ascii OFF format.
  /// The operator only reads the point property and does not read files 
  /// with vertex normals or textures.
  /// \pre `operator>>(std::istream&,const P&)` must be defined.
  template <typename P>
  std::istream& operator>>(std::istream& is, Surface_mesh<P>& sm)
  {
    typedef Surface_mesh<P> Mesh;
    typedef typename Mesh::size_type size_type;
    sm.clear();
    int n, f, e;
    std::string off;
    is >> off;
    assert(off == "OFF" || off == "COFF");
    is >> n >> f >> e;
    sm.reserve(n,2*f,e);
    P p;
    for(int i=0; i < n; i++){
      is >> sm_skip_comments;
      is >> p;
      sm.add_vertex(p);
      is.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
    }
    std::vector<size_type> vr;
    std::size_t d;
    for(int i=0; i < f; i++){
      is >> sm_skip_comments;
      is >> d;
      vr.resize(d);
      for(std::size_t j=0; j<d; j++){
        is >> vr[j];
      }
      sm.add_face(vr);
    }
    return is;
  }

 
 /*! @} */

template <typename P>
Surface_mesh<P>::
Surface_mesh()
  : pmap_(boost::fusion::make_pair< typename Surface_mesh::Vertex_index >(&Surface_mesh::vprops_)
          , boost::fusion::make_pair< typename Surface_mesh::Halfedge_index >(&Surface_mesh::hprops_)
          , boost::fusion::make_pair< typename Surface_mesh::Edge_index >(&Surface_mesh::eprops_)
          , boost::fusion::make_pair< typename Surface_mesh::Face_index >(&Surface_mesh::fprops_))
{
    // allocate standard properties
    // same list is used in operator=() and assign()
    vconn_    = add_property_map<Vertex_index, Vertex_connectivity>("v:connectivity").first;
    hconn_    = add_property_map<Halfedge_index, Halfedge_connectivity>("h:connectivity").first;
    fconn_    = add_property_map<Face_index, Face_connectivity>("f:connectivity").first;
    vpoint_   = add_property_map<Vertex_index, Point>("v:point").first;
    vremoved_ = add_property_map<Vertex_index, bool>("v:removed", false).first;
    eremoved_ = add_property_map<Edge_index, bool>("e:removed", false).first;
    fremoved_ = add_property_map<Face_index, bool>("f:removed", false).first;

    removed_vertices_ = removed_edges_ = removed_faces_ = 0;
    vertices_freelist_ = edges_freelist_ = faces_freelist_ = -1;
    garbage_ = false;
}


//-----------------------------------------------------------------------------
template <typename P>
Surface_mesh<P>&
Surface_mesh<P>::
operator=(const Surface_mesh<P>& rhs)
{
    if (this != &rhs)
    {
        pmap_ = rhs.pmap_;

        // deep copy of property containers
        vprops_ = rhs.vprops_;
        hprops_ = rhs.hprops_;
        eprops_ = rhs.eprops_;
        fprops_ = rhs.fprops_;

        // property handles contain pointers, have to be reassigned
        vconn_    = property_map<Vertex_index, Vertex_connectivity>("v:connectivity").first;
        hconn_    = property_map<Halfedge_index, Halfedge_connectivity>("h:connectivity").first;
        fconn_    = property_map<Face_index, Face_connectivity>("f:connectivity").first;
        vremoved_ = property_map<Vertex_index, bool>("v:removed").first;
        eremoved_ = property_map<Edge_index, bool>("e:removed").first;
        fremoved_ = property_map<Face_index, bool>("f:removed").first;
        vpoint_   = property_map<Vertex_index, P>("v:point").first;

        // how many elements are removed?
        removed_vertices_  = rhs.removed_vertices_;
        removed_edges_     = rhs.removed_edges_;
        removed_faces_     = rhs.removed_faces_;
        vertices_freelist_ = rhs.vertices_freelist_;
        edges_freelist_    = rhs.edges_freelist_;
        faces_freelist_    = rhs.faces_freelist_;
        garbage_           = rhs.garbage_;
    }

    return *this;
}


//-----------------------------------------------------------------------------
template <typename P>
Surface_mesh<P>&
Surface_mesh<P>::
assign(const Surface_mesh<P>& rhs)
{
    if (this != &rhs)
    {
        // clear properties
        vprops_.clear();
        hprops_.clear();
        eprops_.clear();
        fprops_.clear();

        // allocate standard properties
        vconn_    = add_property_map<Vertex_index, Vertex_connectivity>("v:connectivity").first;
        hconn_    = add_property_map<Halfedge_index, Halfedge_connectivity>("h:connectivity").first;
        fconn_    = add_property_map<Face_index, Face_connectivity>("f:connectivity").first;
        vpoint_   = add_property_map<Vertex_index, P>("v:point").first;
        vremoved_ = add_property_map<Vertex_index, bool>("v:removed", false).first;
        eremoved_ = add_property_map<Edge_index, bool>("e:removed", false).first;
        fremoved_ = add_property_map<Face_index, bool>("f:removed", false).first;

        // copy properties from other mesh
        vconn_.array()     = rhs.vconn_.array();
        hconn_.array()     = rhs.hconn_.array();
        fconn_.array()     = rhs.fconn_.array();
        vpoint_.array()    = rhs.vpoint_.array();
        vremoved_.array()  = rhs.vremoved_.array();
        eremoved_.array()  = rhs.eremoved_.array();
        fremoved_.array()  = rhs.fremoved_.array();

        // resize (needed by property containers)
        vprops_.resize(rhs.num_vertices());
        hprops_.resize(rhs.num_halfedges());
        eprops_.resize(rhs.num_edges());
        fprops_.resize(rhs.num_faces());

        // how many elements are removed?
        removed_vertices_  = rhs.removed_vertices_;
        removed_edges_     = rhs.removed_edges_;
        removed_faces_     = rhs.removed_faces_;
        vertices_freelist_ = rhs.vertices_freelist_;
        edges_freelist_    = rhs.edges_freelist_;
        faces_freelist_    = rhs.faces_freelist_;
        garbage_           = rhs.garbage_;
    }

    return *this;
}

//-----------------------------------------------------------------------------
template <typename P>
void
Surface_mesh<P>::
clear()
{
    vprops_.resize(0);
    hprops_.resize(0);
    eprops_.resize(0);
    fprops_.resize(0);

    vprops_.shrink_to_fit();
    hprops_.shrink_to_fit();
    eprops_.shrink_to_fit();
    fprops_.shrink_to_fit();

    removed_vertices_ = removed_edges_ = removed_faces_ = 0;
    vertices_freelist_ = edges_freelist_ = faces_freelist_ = -1;
    garbage_ = false;
}

//-----------------------------------------------------------------------------
/// @cond CGAL_DOCUMENT_INTERNALS
template <typename P>
void
Surface_mesh<P>::
property_stats(std::ostream& out) const
{
    std::vector<std::string> props;

    out << "vertex properties:\n";
    props = properties<Vertex_index>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;

    out << "halfedge properties:\n";
    props = properties<Halfedge_index>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;

    out << "edge properties:\n";
    props = properties<Edge_index>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;

    out << "face properties:\n";
    props = properties<Face_index>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;
}
/// @endcond

//-----------------------------------------------------------------------------
template <typename P>
typename Surface_mesh<P>::Halfedge_index
Surface_mesh<P>::
halfedge(Vertex_index source, Vertex_index target) const
{
    CGAL_assertion(has_valid_index(source) && has_valid_index(target));

    Halfedge_index h  = halfedge(target);
    const Halfedge_index hh = h;

    if (h.is_valid())
    {
        do
        {
            if (this->source(h) == source)
              return h;
            h = next_around_target(h);
        }
        while (h != hh);
    }

    return Halfedge_index();
}


//-----------------------------------------------------------------------------
template <typename P>
void
Surface_mesh<P>::
adjust_incoming_halfedge(Vertex_index v)
{
    Halfedge_index h  = halfedge(v);
    Halfedge_index hh = h;

    if (h.is_valid())
    {
        if (target(h) != v)
        {
            // wrong target, flip
            h = opposite(h);
            hh = h;
            set_halfedge(v, h);
        }

        do
        {
            if (is_border(h))
            {
                set_halfedge(v, h);
                return;
            }
            h = next_around_target(h);
        }
        while (h != hh);
    }
}

//-----------------------------------------------------------------------------

 /// @cond CGAL_DOCUMENT_INTERNALS

template <typename P>
template <typename Range>
typename Surface_mesh<P>::Face_index
Surface_mesh<P>::add_face(const Range& r)
{
  std::vector<Vertex_index> vertices(r.begin(), r.end()); // quick and dirty copy
    Vertex_index                   v;
    unsigned int             i, ii, n((int)vertices.size()), id;
    std::vector<Halfedge_index>    halfedges(n);
    std::vector<bool>        is_new(n), needs_adjust(n, false);
    Halfedge_index                 inner_next, inner_prev,
    outer_next, outer_prev,
    border_next, border_prev,
    patch_start, patch_end;

    // cache for set_next and vertex' set_halfedge
    typedef std::pair<Halfedge_index, Halfedge_index>  NextCacheEntry;
    typedef std::vector<NextCacheEntry>    NextCache;

    NextCache    next_cache;
    next_cache.reserve(3*n);


    // don't allow degenerated faces
    CGAL_assertion (n > 2);

    // test for topological errors
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
      if ( !(is_isolated(vertices[i]) || is_border(vertices[i]) ) )
        {
            std::cerr << "Surface_meshT::add_face: complex vertex " << vertices[i] << std::endl;
            return Face_index();
        }
        halfedges[i] = halfedge(vertices[i], vertices[ii]);
        is_new[i]    = !halfedges[i].is_valid();
        if (!is_new[i] && !is_border(halfedges[i]))
        {
            std::cerr << "Surface_meshT::add_face: complex edge " << halfedges[i] << std::endl;
            std::cerr << std::boolalpha << is_border(halfedges[i]) << std::endl;
            std::cerr << target(halfedges[i]) << std::endl;

            return Face_index();
        }
    }

    // re-link patches if necessary
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        if (!is_new[i] && !is_new[ii])
        {
            inner_prev = halfedges[i];
           inner_next = halfedges[ii];

            if (next(inner_prev) != inner_next)
            {
                // here comes the ugly part... we have to relink a whole patch

                // search a free gap
                // free gap will be between border_prev and border_next
                outer_prev = opposite(inner_next);
                outer_next = opposite(inner_prev);
                border_prev = outer_prev;
                do
                    border_prev = opposite(next(border_prev));
                while (!is_border(border_prev) || border_prev==inner_prev);
                border_next = next(border_prev);
                CGAL_assertion(is_border(border_prev));
                CGAL_assertion(is_border(border_next));


                // ok ?
                if (border_next == inner_next)
                {
                    std::cerr << "Surface_meshT::add_face: patch re-linking failed\n";
                    return Face_index();
                }

                // other halfedges' indices
                patch_start = next(inner_prev);
                patch_end   = prev(inner_next);

                // relink
                next_cache.push_back(NextCacheEntry(border_prev, patch_start));
                next_cache.push_back(NextCacheEntry(patch_end, border_next));
                next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
            }
        }
    }
    // create missing edges
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
      if (is_new[i]){
            halfedges[i] = add_edge(vertices[i], vertices[ii]);
            set_face(opposite(halfedges[i]), null_face()); // as it may be recycled we have to reset it  
            assert(source(halfedges[i]) == vertices[i]);
      }
    // create the face
    Face_index f = add_face();
    set_halfedge(f, halfedges[n-1]);

    // setup halfedges
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        v          = vertices[ii];
        inner_prev = halfedges[i];
        inner_next = halfedges[ii];

        id = 0;
        if (is_new[i])  id |= 1;
        if (is_new[ii]) id |= 2;

        if (id)
        {
            outer_prev = opposite(inner_next);
            outer_next = opposite(inner_prev);

            // set outer links
            switch (id)
            {
                case 1: // prev is new, next is old

                    border_prev = prev(inner_next);
                    next_cache.push_back(NextCacheEntry(border_prev, outer_next));
                    set_halfedge(v, border_prev);
                    break;

                case 2: // next is new, prev is old

                    border_next = next(inner_prev);
                    next_cache.push_back(NextCacheEntry(outer_prev, border_next));
                    set_halfedge(v, outer_prev);
                    break;

                case 3: // both are new

                    if (!halfedge(v).is_valid())
                    {
                        set_halfedge(v, outer_prev);
                        next_cache.push_back(NextCacheEntry(outer_prev, outer_next));
                    }
                    else
                    {
                        border_prev = halfedge(v);
                        border_next = next(border_prev);
                        next_cache.push_back(NextCacheEntry(border_prev, outer_next));
                        next_cache.push_back(NextCacheEntry(outer_prev, border_next));
                    }
                    break;
            }

            // set inner link
            next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
        }
        else {
            needs_adjust[ii] = true; // (halfedge(v) == inner_next); the code is over adjusting
        }

        // set face index
        set_face(halfedges[i], f);
    }



    // process next halfedge cache
    typename NextCache::const_iterator ncIt(next_cache.begin()), ncEnd(next_cache.end());
    for (; ncIt != ncEnd; ++ncIt)
        set_next(ncIt->first, ncIt->second);



    // adjust vertices' halfedge index
    for (i=0; i<n; ++i)
      if (true) //(needs_adjust[i])
            adjust_incoming_halfedge(vertices[i]);

    return f;
}

  /// @endcond

//-----------------------------------------------------------------------------
template <typename P>
typename Surface_mesh<P>::size_type
Surface_mesh<P>::
degree(Vertex_index v) const
{
    size_type count(0);
    if(halfedge(v) == null_halfedge()){
      return 0;
    }
    Vertex_around_target_circulator vvit(halfedge(v), *this);
    Vertex_around_target_circulator vvend = vvit;
    if(vvit) do
    {
        ++count;
    } while (++vvit != vvend);

    return count;
}


//-----------------------------------------------------------------------------
template <typename P>
typename Surface_mesh<P>::size_type
Surface_mesh<P>::
degree(Face_index f) const
{
    size_type count(0);
    if(halfedge(f) == null_halfedge()){
      return 0;
    }
    Vertex_around_face_circulator fvit(halfedge(f),*this);
    Vertex_around_face_circulator fvend = fvit;
    if(fvit) do {
        ++count;
    } while (++fvit != fvend);

    return count;
}

template <typename P>
void
Surface_mesh<P>::
collect_garbage()
{
    int  i, i0, i1,
    nV(num_vertices()),
    nE(num_edges()),
    nH(num_halfedges()),
    nF(num_faces());

    Vertex_index    v;
    Halfedge_index  h;
    Face_index      f;


    // setup index mapping%
    Property_map<Vertex_index, Vertex_index>      vmap = add_property_map<Vertex_index, Vertex_index>("v:garbage-collection").first;
    Property_map<Halfedge_index, Halfedge_index>  hmap = add_property_map<Halfedge_index, Halfedge_index>("h:garbage-collection").first;
    Property_map<Face_index, Face_index>          fmap = add_property_map<Face_index, Face_index>("f:garbage-collection").first;
    for (i=0; i<nV; ++i)
        vmap[Vertex_index(i)] = Vertex_index(i);
    for (i=0; i<nH; ++i)
        hmap[Halfedge_index(i)] = Halfedge_index(i);
    for (i=0; i<nF; ++i)
        fmap[Face_index(i)] = Face_index(i);



    // really remove vertices
    if (nV > 0)
    {
        i0=0;  i1=nV-1;

        while (1)
        {
            // find first removed and last un-removed
            while (!vremoved_[Vertex_index(i0)] && i0 < i1)  ++i0;
            while ( vremoved_[Vertex_index(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            vprops_.swap(i0, i1);
        };

        // remember new size
        nV = vremoved_[Vertex_index(i0)] ? i0 : i0+1;
    }

    // really remove edges
    if (nE > 0)
    {
        i0=0;  i1=nE-1;

        while (1)
        {
            // find first removed and last un-removed
            while (!eremoved_[Edge_index(i0)] && i0 < i1) ++i0;
            while ( eremoved_[Edge_index(i1)] && i0 < i1) --i1;
            if (i0 >= i1) break;

            // swap
            eprops_.swap(i0, i1);
            hprops_.swap(2*i0,   2*i1);
            hprops_.swap(2*i0+1, 2*i1+1);
        };

        // remember new size
        nE = eremoved_[Edge_index(i0)] ? i0 : i0+1;
        nH = 2*nE;
    }


    // really remove faces
    if (nF > 0)
    {
        i0=0;  i1=nF-1;

        while (1)
        {
            // find 1st removed and last un-removed
            while (!fremoved_[Face_index(i0)] && i0 < i1)  ++i0;
            while ( fremoved_[Face_index(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            fprops_.swap(i0, i1);
        };

        // remember new size
        nF = fremoved_[Face_index(i0)] ? i0 : i0+1;
    }


    // update vertex connectivity
    for (i=0; i<nV; ++i)
    {
        v = Vertex_index(i);
        if (!is_isolated(v))
            set_halfedge(v, hmap[halfedge(v)]);
    }


    // update halfedge connectivity
    for (i=0; i<nH; ++i)
    {
        h = Halfedge_index(i);
        set_target(h, vmap[target(h)]);
        set_next(h, hmap[next(h)]);
        if (!is_border(h))
            set_face(h, fmap[face(h)]);
    }


    // update indices of faces
    for (i=0; i<nF; ++i)
    {
        f = Face_index(i);
        set_halfedge(f, hmap[halfedge(f)]);
    }

    // remove index maps
    remove_property_map<Vertex_index>(vmap);
    remove_property_map<Halfedge_index>(hmap);
    remove_property_map<Face_index>(fmap);

    // finally resize arrays
    vprops_.resize(nV); vprops_.shrink_to_fit();
    hprops_.resize(nH); hprops_.shrink_to_fit();
    eprops_.resize(nE); eprops_.shrink_to_fit();
    fprops_.resize(nF); fprops_.shrink_to_fit();

    removed_vertices_ = removed_edges_ = removed_faces_ = 0;
    vertices_freelist_ = edges_freelist_ = faces_freelist_ = -1;
    garbage_ = false;
}

} // CGAL

#endif /* CGAL_SURFACE_MESH_H */

