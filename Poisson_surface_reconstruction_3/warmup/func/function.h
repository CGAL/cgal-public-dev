#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;

template <typename K, typename Point, typename T>
class Func{
private:
  T* m_tr;
  double m_isovalue;
public:
  Func(T* t, double isovalue): m_tr(t), m_isovalue(isovalue){}
  ~Func(){}

  FT operator()(Point query) const{
    return m_tr->compute_func_value(query) - m_isovalue;
  }
};

template <typename K, typename Point, typename T>
class FuncSmooth{
private:
  T* m_tr;
  double m_isovalue;
public:
  FuncSmooth(T* t, double isovalue): m_tr(t), m_isovalue(isovalue){}
  ~FuncSmooth(){}

  FT operator()(Point query) const{
    return m_tr->compute_func_value_BB(query) - m_isovalue;
  }
};
