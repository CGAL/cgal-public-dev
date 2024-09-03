#ifndef CGAL_ENVVOR_VORONOI_DIAGRAM_2_H
#define CGAL_ENVVOR_VORONOI_DIAGRAM_2_H

#include <CGAL/Envelope_3/Envelope_diagram_on_surface_2.h>
#include <CGAL/Voronoi_2_to_Envelope_3_adaptor.h>

#if defined(FILTERED_ZONE)
#include <CGAL/Envelope_surface_id_traits_2.h>
#endif

namespace CGAL {

namespace Envelope_voronoi_2 {

// The ifdef is temporary until we'll make the envelope diagram to contain
// the xy-monotone surfaces inside and then we can compare pointers to them.
// the Envelope_surface_id_traits_2 is used to set an id to each of the
// surfaces so we can compare them without the need of expensively comparing
// them.

/*! \todo The Envelope_surface_id_traits_2 should be removed and a regular
  comparison should be used */

template <typename GeomTraits_,
          class TopTraits_ = typename Default_planar_topology
          <GeomTraits_, Arr_default_dcel<GeomTraits_> >::Traits>
class Voronoi_diagram_2 :
    public Envelope_diagram_on_surface_2<
#if defined(FILTERED_ZONE)
    Envelope_surface_id_traits_2< Voronoi_2_to_Envelope_3_adaptor<GeomTraits_>>,
#else
    Voronoi_2_to_Envelope_3_adaptor< GeomTraits_ >,
#endif
    typename TopTraits_::template rebind
    <
#if defined(FILTERED_ZONE)
      Envelope_surface_id_traits_2< Voronoi_2_to_Envelope_3_adaptor<GeomTraits_> >,
#else
      Voronoi_2_to_Envelope_3_adaptor<GeomTraits_>,
#endif
      Envelope_3::Envelope_pm_dcel
      <
#if defined(FILTERED_ZONE)
        Envelope_surface_id_traits_2<Voronoi_2_to_Envelope_3_adaptor<GeomTraits_ > >,
        typename Envelope_surface_id_traits_2< Voronoi_2_to_Envelope_3_adaptor<
          GeomTraits_>>::
          Xy_monotone_surface_3
#else
        Voronoi_2_to_Envelope_3_adaptor< GeomTraits_>,
        typename Voronoi_2_to_Envelope_3_adaptor<GeomTraits_>::
          Xy_monotone_surface_3
#endif
      >
    >::other
  >
{
  using Geom_traits_2 = GeomTraits_;
  using Top_traits = TopTraits_;
  using Adapt_geom_traits_2 = Voronoi_2_to_Envelope_3_adaptor< Geom_traits_2 >;
  using Adapt_top_traits =
    typename Top_traits::template rebind
    <Adapt_geom_traits_2, Envelope_3::Envelope_pm_dcel
     <Adapt_geom_traits_2,
      typename Adapt_geom_traits_2::Xy_monotone_surface_3>>::other;

public:
  using Base =
    Envelope_diagram_on_surface_2<Adapt_geom_traits_2, Adapt_top_traits>;
};

template <class GeomTraits_>
class Spherical_voronoi_diagram_2 :
    public Voronoi_diagram_2< GeomTraits_,
                              Arr_spherical_topology_traits_2< GeomTraits_>>
{};

/*
    Voronoi_2_to_Envelope_3_adaptor< GeomTraits_ >,

    Envelope_3::Envelope_pm_dcel< Voronoi_2_to_Envelope_3_adaptor< GeomTraits_ >,
      typename Voronoi_2_to_Envelope_3_adaptor< GeomTraits_ >::
        Xy_monotone_surface_3 >,
      typename GeomTraits_::Has_boundary_category >::Traits >

 */
}

} //namespace CGAL

#endif
