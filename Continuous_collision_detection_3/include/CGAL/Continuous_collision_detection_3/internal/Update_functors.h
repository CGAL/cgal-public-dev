// Copyright (c) 2023 GeometryFactory (France).
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef UPDATE_FUNCTORS_H
#define UPDATE_FUNCTORS_H

namespace CGAL {
namespace Collisions {
namespace internal {

template <class CollisionScene>
struct Swap_current_next_functor {

    using Vertex_index  = typename CollisionScene::Scene_vertex_index;
    using Mesh          = typename CollisionScene::Mesh;


    Swap_current_next_functor() {}

    void operator() (Mesh* mesh, Vertex_index svi) {
        swap(
          mesh->point(svi.local_index()),
          mesh->next_point(svi.local_index())
        );
    }
};

/// \ingroup PkgCollisions3Ref
/// \cgalModels `CollisionSceneUpdateFunctor`
/// blablabla
template <class CollisionScene>
struct Translate_functor {

    using K             = typename CollisionScene::K;
    using Mesh_index    = typename CollisionScene::Mesh_index;
    using Vertex_index  = typename CollisionScene::Scene_vertex_index;
    using Vector        = typename CollisionScene::Vector;
    using Mesh          = typename CollisionScene::Mesh;
    using Transform     = typename K::Aff_transformation_3;

    Transform translation;
    Mesh_index mi{0};

    Translate_functor(Mesh_index mi, Vector translation_vector)
        : translation(CGAL::TRANSLATION, translation_vector)
        , mi{mi}
        {}

    void operator() (Mesh* mesh, Vertex_index svi) {
        if( mi == svi.mesh_index())
            {
            mesh->point( svi.local_index() ) = mesh->point( svi.local_index() ).transform(this->translation),
            swap(
                mesh->point( svi.local_index() ),
                mesh->next_point( svi.local_index() )
            );
        }
    }
};

template <class CollisionScene>
struct Contraction_functor {

    using K             = typename CollisionScene::K;
    using Mesh_index    = typename CollisionScene::Mesh_index;
    using Vertex_index  = typename CollisionScene::Scene_vertex_index;
    using Point         = typename CollisionScene::Point;
    using Vector        = typename CollisionScene::Vector;
    using Mesh          = typename CollisionScene::Mesh;
    using Transform     = typename K::Aff_transformation_3;

    Point contraction_point;
    Mesh_index mi{0};
    double t{1};

    Contraction_functor(Mesh_index mi, Point contraction_point, double t)
        : contraction_point{contraction_point}
        , mi{mi}
        , t{t}
        {}

    void operator() (Mesh* mesh, Vertex_index svi) {
        if( mi == svi.mesh_index())
        {
            Point  p = mesh->point( svi.local_index() );
            Vector v = t*(contraction_point - p);
            mesh->next_point( svi.local_index() ) = p+v;
        }
    }
};


}
}
}

#endif