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
//                 Andreas von Dziegielewski (dziegiel@uni-mainz.de)
//
// ================================================================================


#ifndef HASHSET_VOXEL_H_
#define HASHSET_VOXEL_H_

#include <CGAL/tuple.h>
#include <CGAL/Real_timer.h> 
#include <CGAL/ipower.h> 
#include <CGAL/Handle_with_policy.h>


#include <boost/unordered_set.hpp>
#include <boost/foreach.hpp>

namespace SV {
namespace internal {
template <typename INT = short> 
class HashSet_Voxel_rep {
//public:
//  mutable CGAL::Real_timer insertion_timer;
//  mutable CGAL::Real_timer compression_timer;
//  mutable CGAL::Real_timer erase_timer;
//  mutable CGAL::Real_timer contains_timer;

public:
  typedef CGAL::cpp0x::tuple <INT, INT, INT> Voxel; 
  // Hashing scheme (should become a template parameter?)
  typedef CGAL::cpp0x::tuple <INT, INT, INT> Key;
  typedef unsigned int Hash_value;
  struct Hash_funtion : public std::unary_function <Key, Hash_value> {

    template <class T>
    Hash_value hash(const T& x) const {
      return boost::hash<T>()(x);
    }
      
    Hash_value
    operator()(const Key& key) const
    {
      return Hash_value(    
          CGAL::cpp0x::get<0>(key)+
          CGAL::cpp0x::get<1>(key)*1000+
          CGAL::cpp0x::get<2>(key)*100000);
        
    }
  };
  
  struct Key_equal : public std::binary_function <Key, Key, bool> {
    bool
    operator()(const Key& k1, const Key&k2) const
    {
      if(k1 == k2){
        return true; 
      }else{
        return false;
      }
    }
  };

public:
  typedef boost::unordered_set <Key, Hash_funtion, Key_equal> USet; 
  typedef typename USet::size_type size_type; 

  typedef typename USet::iterator iterator; 
  iterator begin(){return uset.begin();}
  iterator end(){return uset.end();}
  

  // Make it a container: 
  typedef typename USet::const_iterator  const_iterator; 
  typedef typename USet::const_reference const_reference; 
  const_iterator begin() const {return uset.begin();}
  const_iterator end() const {return uset.end();}

private: // member
  
  USet uset;    
  
public: // make protected once debugging is completed 
  USet& set() { return uset; }
  const USet& set() const { return uset; }

public:  
  
  size_type     size()       const { return uset.size();}
  bool          empty()      const { return uset.empty();}
  void          clear()            { uset.clear();}
 
public:
  
  bool
  contains(INT x, INT y, INT z) const
  {
    return this->contains(Key(x,y,z));
  }
  
  bool contains(Key k) const
  {
    if (uset.count(k) == 1) return true;
    return false;
  }

  bool insert(const CGAL::cpp0x::tuple<INT,INT,INT>& voxel)
  {
    return insert(
        CGAL::cpp0x::get<0>(voxel), 
        CGAL::cpp0x::get<1>(voxel), 
        CGAL::cpp0x::get<2>(voxel));
  }
  
  bool
  insert(INT x, INT y, INT z)
  {
    
    assert( x >= 0 );
    assert( y >= 0 );
    assert( z >= 0 );

    if (this->contains(x, y, z)) {
      return false;
    }
    uset.insert(Key(x, y, z) );
    return true;
  }

  bool erase(const CGAL::cpp0x::tuple<INT,INT,INT>& voxel)
  {
    return erase(
        CGAL::cpp0x::get<0>(voxel), 
        CGAL::cpp0x::get<1>(voxel), 
        CGAL::cpp0x::get<2>(voxel));
  }
  
  bool
  erase(INT x, INT y, INT z)
  {
    if (!this->contains(x, y, z)) {
      return false;
    }
    uset.erase(Key(x, y, z) );
    return true;
  }
};
} // namespace internal 

template <typename INT_ > 
class HashSet_Voxel
  : public CGAL::Handle_with_policy< internal::HashSet_Voxel_rep<INT_> >    
{
public:
  typedef INT_ INT; 
  typedef internal::HashSet_Voxel_rep<INT> Rep;
  typedef CGAL::Handle_with_policy< Rep > Base;
  typedef HashSet_Voxel<INT> Self; 
  
  typedef typename Rep::USet USet; 
  typedef typename Rep::size_type size_type; 
  typedef typename Rep::Voxel Voxel;


  HashSet_Voxel():Base(){}
  HashSet_Voxel(const Self& oct) : Base(static_cast<const Base&>(oct)){}
  template <class InputIterator> 
  HashSet_Voxel(InputIterator begin, InputIterator end): Base(){
    for(InputIterator it=begin; it != end; it++){
       insert(*it);
    }
  }

  
public:
  struct Inserter{
    HashSet_Voxel* m_hset; 
    Inserter(HashSet_Voxel* hset):m_hset(hset){};
    Inserter& operator++(){return *this;}
    Inserter& operator++(int){return *this;}
    Inserter& operator*(){return *this;} 
    Inserter& operator=( const Voxel& v) {
      this->m_hset->insert(v);
      return *this; 
    }
  };  

  struct Eraser{
    HashSet_Voxel* m_hset; 
    Eraser(HashSet_Voxel* hset):m_hset(hset){};
    Eraser& operator++(){return *this;}
    Eraser& operator++(int){return *this;}
    Eraser& operator*(){return *this;} 
    Eraser& operator=( const Voxel& v) {
      this->m_hset->erase(v);
      return *this; 
    }
  };
  
public:
  Inserter inserter(){return Inserter(this);}
  Eraser   eraser()  {return Eraser(this);}

  
  
  typedef typename USet::const_reference const_reference; 
  
  USet& set() { this->copy_on_write(); return this->ptr()->set(); }
  const USet& set() const { return this->ptr()->set(); }

  typedef typename USet::iterator iterator; 
  iterator begin(){ this->copy_on_write(); return this->ptr()->begin();}
  iterator end()  { this->copy_on_write(); return this->ptr()->end();}
  typedef typename USet::const_iterator const_iterator; 
  const_iterator begin() const { return this->ptr()->begin();}
  const_iterator end()   const { return this->ptr()->end();}

  size_type size()  const { return this->ptr()->size(); }
  bool      empty() const { return this->ptr()->empty();}
  void      clear()       { this->copy_on_write(); return this->ptr()->clear();}
  
  bool insert(const CGAL::cpp0x::tuple<INT,INT,INT>& voxel){
    this->copy_on_write();
    return this->ptr()->insert(voxel);
  }
  bool insert(INT x, INT y, INT z){
    this->copy_on_write();
    return this->ptr()->insert(x,y,z);
  }
  bool contains(const Voxel& v) const{
    return this->ptr()->contains(v);
  }
  bool contains(INT x, INT y, INT z, int on_level) const{
    return this->ptr()->contains(x,y,z,on_level);
  }
  bool contains(INT x, INT y, INT z) const{
    return this->ptr()->contains(x,y,z);
  }
  
  bool erase(const Voxel& v){
    using CGAL::cpp0x::get; 
    return erase(get<0>(v),get<1>(v),get<2>(v));
  }
  bool erase(INT x, INT y, INT z){
    this->copy_on_write();
    return this->ptr()->erase(x,y,z);
  }  
};




} // namespace SV

#endif /* HASH_SET_VOXEL_H_ */
