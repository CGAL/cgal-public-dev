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


#ifndef SV_OCTREE_3_H
#define SV_OCTREE_3_H

#include <CGAL/tuple.h>
#include <CGAL/Real_timer.h> 
#include <CGAL/ipower.h> 
#include <CGAL/Handle_with_policy.h>


#include <boost/unordered_set.hpp>
#include <boost/foreach.hpp>

namespace SV {

namespace internal {

template <typename INT> 
class Octree_3_rep {
  
public:
  mutable CGAL::Real_timer insertion_timer;
  mutable CGAL::Real_timer compression_timer;
  mutable CGAL::Real_timer erase_timer;
  mutable CGAL::Real_timer contains_timer;
  
public:
  // Hashing scheme (should become a template parameter?)
  typedef CGAL::cpp0x::tuple <INT, INT, INT, int> Key;
  typedef CGAL::cpp0x::tuple <INT, INT, INT>      Voxel;
  typedef unsigned int Hash_value;
  struct Hash_funtion : public std::unary_function <Key, Hash_value> {

    template <class T>
    Hash_value hash(const T& x) const {
      return boost::hash<T>()(x);
    }
      
    Hash_value
    operator()(const Key& key) const
    {
      using CGAL::cpp0x::get;
      // do we need to improve this function? 
#if 0
//      Hash_value result (    
//             get<3>(key)+
//             ((get<0>(key)&255)<<2)+
//             ((get<1>(key)&255)<<10)+
//             ((get<2>(key)&255)<<18));
      //std::cerr << result << " " <<  std::endl;    
      Hash_value result (    
          get<3>(key)+
          ((get<0>(key)&1023)<<2)+
          ((get<1>(key)&1023)<<12)+
          ((get<2>(key)&1023)<<22));
      //std::cerr << result << " " <<  std::endl; 
      return result; 
#else
      return Hash_value(    
          get<0>(key)+
          get<1>(key)*100+
          get<2>(key)*10000+
          get<3>(key)*1000000);
        
#endif 
    }
  };
  struct Key_equal : public std::binary_function <Key, Key, bool> {
    bool
    operator()(const Key& k1, const Key&k2) const
    {
      if(k1 == k2){
        //std::cerr <<"T"; 
        return true; 
      }else{
        // std::cerr <<"C"; 
        return false;
      }
    }
  };

public:
  typedef boost::unordered_set <Key, Hash_funtion, Key_equal> USet; 
  typedef typename USet::size_type size_type; 

private: // member
  
  USet uset;    
  int m_resolution;
  double m_scale;
  int m_first_level; 
  
public: // make protected once debugging is completed 
  USet& set() { return uset; }
  const USet& set() const { return uset; }

public:  
 
  
  size_type     size()        const { return uset.size();}
  bool          empty()       const { return uset.empty();}
  int           resolution()  const { return m_resolution;} 
  double        scale ()      const { return m_scale;}
  void          clear()             { uset.clear();}
  int           first_level() const { return m_first_level;}
 
  
  void set_resolution(int resolution){
    m_resolution = resolution;
    m_scale = (1 << resolution);
    if(this->empty()) m_first_level = resolution-1; 
  }
public:
  // constructors
  Octree_3_rep(int resolution = 1){
    set_resolution(resolution);
    m_first_level = resolution-1;
  }
//   template<class Iterator> 
//   Octree_3_rep(int resolution, Iterator begin, Iterator end){
//     set_resolution(resolution);
//     for(Iterator it = begin; it != end; it++){
//       this.insert(*it);
//     }
//   }
  
private:
  void compress(INT xx, INT yy, INT zz){
    using CGAL::cpp0x::get;

    CGAL_precondition(this->contains(xx, yy, zz)); 
    // std::cerr << x << " " << y << " " << z << std::endl;
    // start compression:
    for (int i = 1; i <= m_resolution; i++) {
      int x = xx >> i;  
      int y = yy >> i; 
      int z = zz >> i;
      Key new_key(x,y,z, m_resolution - 1 - i );
      x <<= 1;  y <<= 1; z <<= 1;
      
      //std::cerr << x << " " << y << " " << z << std::endl;
      int level = m_resolution - i; 
      int xp1   = x + 1; 
      int yp1   = y + 1;
      int zp1   = z + 1; 
      if (uset.count(Key(x    , y    , z    , level)) == 1 && 
          uset.count(Key(x    , y    , zp1  , level)) == 1 && 
          uset.count(Key(x    , yp1  , z    , level)) == 1 && 
          uset.count(Key(x    , yp1  , zp1  , level)) == 1 && 
          uset.count(Key(xp1  , y    , z    , level)) == 1 && 
          uset.count(Key(xp1  , y    , zp1  , level)) == 1 && 
          uset.count(Key(xp1  , yp1  , z    , level)) == 1 && 
          uset.count(Key(xp1  , yp1  , zp1  , level)) == 1) {
     
        uset.erase(Key(x    , y    , z    , level));
        uset.erase(Key(x    , y    , zp1  , level));
        uset.erase(Key(x    , yp1  , z    , level));
        uset.erase(Key(x    , yp1  , zp1  , level));
        uset.erase(Key(xp1  , y    , z    , level));
        uset.erase(Key(xp1  , y    , zp1  , level));
        uset.erase(Key(xp1  , yp1  , z    , level));
        uset.erase(Key(xp1  , yp1  , zp1  , level));

        uset.insert(new_key);
        m_first_level = (std::min)(first_level(),get<3>(new_key));
      }else{
        break;
      }
    }

    // std::cerr << uset.count(Key(x,y,z));
  }

public:
  bool insert(const CGAL::cpp0x::tuple<INT,INT,INT>& voxel)
  {
    using CGAL::cpp0x::get;
    return insert(get<0>(voxel), get<1>(voxel), get<2>(voxel));
  }
  
  bool
  insert(INT x, INT y, INT z)
  {

#if 0
    if (!(x >= 0)) return false;
    if (!(y >= 0)) return false;
    if (!(z >= 0)) return false;
#else 
    assert( x >= 0 );
    assert( y >= 0 );
    assert( z >= 0 );
#endif

    if (this->contains(x, y, z)) {
      // std::cerr<< 0; 
      return false;
    }
    //std::cerr<< 1; 
    // Key is not covered by the OCtree, mark new cell
    uset.insert(Key(x, y, z, m_resolution - 1) );
    this->compress(x,y,z);
    return true; 
  }

public:

  // test whether voxel of on a certain level is covered by the octree 
  bool
  contains(INT x, INT y, INT z, int on_level) const
  {
    if(on_level > first_level()) return false; 
    
    int i = m_resolution - 1 - first_level() - on_level;
    while (i >= 0 ) {
      Key key(
          x >> i, 
          y >> i, 
          z >> i, 
          m_resolution - 1 - i + on_level);
      if (uset.count(key) == 1) break;
      i--;
    }
    if (i >= 0 ) return true;
    return false;
  }


  template<class NT> 
  bool contains(const CGAL::cpp0x::tuple<NT,NT,NT>& v) const {
    using CGAL::cpp0x::get; 
    return contains(get<0>(v),get<1>(v),get<2>(v));
  }
  
  bool
  contains(INT x, INT y, INT z) const
  {
    if(empty()) return false; 

    int i = m_resolution - 1 - first_level();
    while (i >= 0 ) {
      Key key( x >> i, y >> i, z >> i, m_resolution - 1 - i );
      if (uset.count(key) == 1) break;
      i--;
    }
    if (i >= 0) return true;
    return false;
  }
  
  
  bool erase(INT xx, INT yy, INT zz){
    using CGAL::cpp0x::get; 
    
    // find the level of the to be erased voxel 
    int i = m_resolution - 1 - first_level();
    Key key( xx >> i, yy >> i, zz >> i, m_resolution - 1 - i ) ; 
    bool result; 
    while (uset.count(key)!=1 && i >= 0) {
      i--;
      key = Key( xx >> i, yy >> i, zz >> i, m_resolution - 1 - i );
    }
    if (i < 0) return false; 
   
    uset.erase(key);
    while(i>0){
      // insert all voxel on finer level 
      int x     = get<0>(key)<<1;
      int y     = get<1>(key)<<1;
      int z     = get<2>(key)<<1;
      int new_level = get<3>(key)+1; 
      int xp1   = x + 1; 
      int yp1   = y + 1;
      int zp1   = z + 1; 
      uset.insert(Key(x    , y    , z    , new_level));
      uset.insert(Key(x    , y    , zp1  , new_level));
      uset.insert(Key(x    , yp1  , z    , new_level));
      uset.insert(Key(x    , yp1  , zp1  , new_level));
      uset.insert(Key(xp1  , y    , z    , new_level));
      uset.insert(Key(xp1  , y    , zp1  , new_level));
      uset.insert(Key(xp1  , yp1  , z    , new_level));
      uset.insert(Key(xp1  , yp1  , zp1  , new_level));
      
      i--;
      key = Key( xx >> i, yy >> i, zz >> i, m_resolution - 1 - i );
      uset.erase(key);
    }
    CGAL_postcondition(!this->contains(xx,yy,zz));
    return true;
  }
};

} // namespace internal 


template <typename INT_ > 
class Octree_3
  : public CGAL::Handle_with_policy< internal::Octree_3_rep<INT_> >    
{
public:
  typedef INT_ INT; 
  typedef internal::Octree_3_rep<INT> Rep;
  typedef CGAL::Handle_with_policy< Rep > Base;
  typedef Octree_3<INT> Self; 
  
  typedef typename Rep::USet USet; 
  typedef typename Rep::size_type size_type; 
  typedef typename Rep::Voxel Voxel;
  
private:
  struct Inserter{
    Octree_3* m_octree; 
    Inserter(Octree_3* octree):m_octree(octree){};
    Inserter& operator++(){return *this;}
    Inserter& operator++(int){return *this;}
    Inserter& operator*(){return *this;} 
    Inserter& operator=( const Voxel& v) {
      this->m_octree->insert(v);
      return *this; 
    }
  };  

  struct Eraser{
    Octree_3* m_octree; 
    Eraser(Octree_3* octree):m_octree(octree){};
    Eraser& operator++(){return *this;}
    Eraser& operator++(int){return *this;}
    Eraser& operator*(){return *this;} 
    Eraser& operator=( const Voxel& v) {
      this->m_octree->erase(v);
      return *this; 
    }
  };
 
public:
  Inserter inserter(){return Inserter(this);}
  Eraser   eraser()  {return Eraser(this);}

public:
  Octree_3():Base(){}
  Octree_3(int resolution):Base(resolution){}
  Octree_3(const Self& oct) : Base(static_cast<const Base&>(oct)){}
  
  USet& set() { this->copy_on_write(); return this->ptr()->set(); }
  const USet& set() const { return this->ptr()->set(); }
  
  size_type size() const { return this->ptr()->size(); }
  void set_resolution(int res){ 
    this->copy_on_write();
    this->ptr()->set_resolution(res);
  }
  
  bool          empty()       const { return this->ptr()->empty();}
  int           resolution()  const { return this->ptr()->resolution();} 
  double        scale ()      const { return this->ptr()->scale();}
  void          clear()             { this->copy_on_write(); return this->ptr()->clear();}
  int           first_level() const { return this->ptr()->first_level();}
  
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

  Octree_3 offset() const {
    typedef typename Rep::Key Key; 
    
    Octree_3 offset = *this;
    offset.copy_on_write(); 

    int  res = resolution() - 1;

    using CGAL::cpp0x::get;
    
    BOOST_FOREACH(Key key, this->set()){
      int range = CGAL::ipower(2,(res - get<3>(key)));
      INT x= (get<0>(key)<<(res - get<3>(key)))-1;
      INT y= (get<1>(key)<<(res - get<3>(key)))-1;
      INT z= (get<2>(key)<<(res - get<3>(key)))-1;

      // inserting the smaller and larger xy wall 
      for(int i = 0; i <= range+1; i++){
        for(int j = 0; j <= range+1; j++){
          offset.insert(x+i,y+j,z); 
          offset.insert(x+i,y+j,z+range+1); 
        }
      }
        
      for(int j = 0; j <= range+1; j++){
        for(int k = 1; k <= range; k++){
          offset.insert(x        ,y+j,z+k); 
          offset.insert(x+range+1,y+j,z+k); 
        }
      }
        
      for(int i = 1; i <= range; i++){
        for(int k = 1; k <= range; k++){
          offset.insert(x+i,y        ,z+k); 
          offset.insert(x+i,y+range+1,z+k); 
        }
      }
    }
    return offset; 
  }

#if 1
  template < class OutputIterator>
  OutputIterator hull(OutputIterator oit) const {
    using CGAL::cpp0x::get;
    
    typedef typename Rep::Key Key; 
    typedef CGAL::cpp0x::tuple<INT,INT,INT> Voxel;
    
    Octree_3 offset = *this; 
    offset.copy_on_write(); 
    
    int  res = resolution() - 1;

#define IN_HULL(x,y,z)                                  \
    if(offset.insert(x,y,z)) *oit++ = Voxel(x,y,z)
    
    BOOST_FOREACH(Key key, this->set()){
      int range = CGAL::ipower(2,(res - get<3>(key)));
      INT x= (get<0>(key)<<(res - get<3>(key)))-1;
      INT y= (get<1>(key)<<(res - get<3>(key)))-1;
      INT z= (get<2>(key)<<(res - get<3>(key)))-1;

      // inserting the smaller and larger xy wall 
      for(int i = 0; i <= range+1; i++){
        for(int j = 0; j <= range+1; j++){
          IN_HULL(x+i,y+j,z); 
          IN_HULL(x+i,y+j,z+range+1); 
        }
      }
        
      for(int j = 0; j <= range+1; j++){
        for(int k = 1; k <= range; k++){
          IN_HULL(x        ,y+j,z+k); 
          IN_HULL(x+range+1,y+j,z+k); 
        }
      }
        
      for(int i = 1; i <= range; i++){
        for(int k = 1; k <= range; k++){
          IN_HULL(x+i,y        ,z+k); 
          IN_HULL(x+i,y+range+1,z+k); 
        }
      }
    }
    return oit;

#undef IN_HULL
  }  
#endif 

  Octree_3 inverse_offset() const {
    typedef typename Rep::Key Key;
    using CGAL::cpp0x::get; 
    
    Octree_3 ioffset = *this;
    ioffset.copy_on_write(); 
    
    std::vector<Voxel> hull; 
    this->hull(std::back_inserter(hull)); 
    BOOST_FOREACH(Voxel v, hull){
      for(int i = -1; i <= 1; i++){
        for(int j = -1; j <= 1; j++){
          for(int k = -1; k <= 1; k++){
            ioffset.erase(get<0>(v)+i,get<1>(v)+j,get<2>(v)+k);
          }
        }
      }
    }
    return ioffset;
  }
  
};

  

} // namespace SV
#endif // SV_OCTREE_3_H 
