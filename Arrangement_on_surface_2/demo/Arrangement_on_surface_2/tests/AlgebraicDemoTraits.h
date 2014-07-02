#ifndef ALGEBRAIC_DEMO_TRAITS_H
#define ALGEBRAIC_DEMO_TRAITS_H
#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#error "CGAL needs CORE."
#endif

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include "ArrangementGraphicsItem.h"
#include "ArrangementDemoTab.h"
#include "DcelTypes.h"

/**
Provided types
-----
#. ArrTraitsType
#. ArrangementType

Provided constants
#. Name - human-readable string naming this traits type
*/
struct AlgebraicDemoTraits
{
  typedef CGAL::CORE_algebraic_number_traits Nt_traits;
  typedef Nt_traits::Integer Coefficient;
  typedef CGAL::Arr_algebraic_segment_traits_2<Coefficient> Algebraic_traits_2;
  typedef Dcel<Algebraic_traits_2> Algebraic_dcel;
  typedef CGAL::Arrangement_with_history_2<Algebraic_traits_2, Algebraic_dcel>
    Algebraic_arrangement_2;
  typedef CGAL::Qt::ArrangementGraphicsItem< Algebraic_arrangement_2 >
    AlgebraicArrangementGraphicsItem;
  typedef ArrangementDemoTab< Algebraic_arrangement_2 > AlgebraicTab;

  // exports
  typedef Algebraic_traits_2 ArrTraitsType;
  typedef Algebraic_arrangement_2 ArrangementType;
  static const std::string Name;
}; // struct AlgebraicDemoTraits
#endif //ALGEBRAIC_DEMO_TRAITS_H
