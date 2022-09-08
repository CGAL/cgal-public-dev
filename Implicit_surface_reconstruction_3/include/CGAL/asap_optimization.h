// Copyright (c) 2007  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Pierre Alliez, Tong Zhao, Hongyi Liu

#ifndef CGAL_IMPLICIT_RECONSTRUCTION_ASAP_OPTIMIZATION_H
#define CGAL_IMPLICIT_RECONSTRUCTION_ASAP_OPTIMIZATION_H

#include <CGAL/license/Implicit_surface_reconstruction_3.h>

#include <boost/shared_ptr.hpp>

#include <unordered_map>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>



namespace CGAL{

template<class Gt, class Tr>
class AsapOptimizer{
public: // types
    typedef Gt Geom_traits;
    typedef Tr Triangulation;
    typedef typename Geom_traits::FT        FT;
    typedef typename Geom_traits::Point_3   Point;
    typedef typename Geom_traits::Vector_3  Vector;

    // Triangulation types
    typedef typename Triangulation::Vertex_handle                Vertex_handle;
    typedef typename Triangulation::Cell_handle                  Cell_handle;
    typedef typename Triangulation::Finite_cells_iterator        Finite_cells_iterator;
    typedef typename Triangulation::Finite_vertices_iterator     Finite_vertices_iterator;

    // Eigen types
    typedef typename Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>     EMatrix;
    typedef typename Eigen::JacobiSVD<EMatrix>                             JacobiSVD;
    typedef typename Eigen::VectorXd                                       EVector;

    typedef Eigen::SparseMatrix<FT>                                                 ESMatrix;
    typedef typename std::vector<Eigen::Triplet<FT> >                               ESTripleList;
    typedef typename Eigen::ConjugateGradient<ESMatrix, Eigen::Lower|Eigen::Upper>  ESSolver;

    typedef std::unordered_map<int, double>         VolumeMap;
    typedef std::unordered_map<int, std::set<int>>  LabelMap;
    typedef LabelMap::iterator                      LabelIter;
    typedef std::unordered_set<int>                 IntSet;

private: // members
    boost::shared_ptr<Triangulation> &  m_tr;
    EMatrix                             m_ref;

public: // functions
    AsapOptimizer(boost::shared_ptr<Triangulation> & tr): m_tr(tr)
    {
        // init ref tetrahedron
        m_ref.resize(3, 4);
        m_ref << -0.5     , 0.5      , 0.       , 0., 
                 -0.288675, -0.288675, 0.57735  , 0.,
                 -0.204125, -0.204125, -0.204125, 0.61237;

        // init triangulation indices
        m_tr->index_all_vertices();
        m_tr->index_all_cells();
    }

    void optimize_asap_linear_solver(IntSet& slivers, 
                                     VolumeMap& vol_map)
    {
        double radius = sqrt((m_tr->bounding_sphere()).squared_radius());
        int num_vertices = m_tr->number_of_vertices();
        int num_cells = m_tr->number_of_finite_cells();
        
        int num_asap_equation = (num_cells - slivers.size()) * 12;
        int num_loc_equation = num_vertices * 3;
        int num_equation = num_asap_equation;

        ESMatrix A(num_equation, 3 * num_vertices);
        EVector  B(num_equation), X(3 * num_vertices);
        B.setZero();
        X.setZero();
        ESTripleList A_triplets(num_equation);
        
        // ASAP equation
        int m = 0; // cell index
        for(Finite_cells_iterator cit = m_tr->finite_cells_begin(); 
            cit != m_tr->finite_cells_end(); 
            cit++)
        {
            // skip slivers
            auto got = slivers.find(cit->info());
            if(got != slivers.end())
                continue;

            // find rotation
            EMatrix pts(3, 4);
            Point center = CGAL::centroid(  cit->vertex(0)->point(), 
                                            cit->vertex(1)->point(), 
                                            cit->vertex(2)->point(),
                                            cit->vertex(3)->point());

            for(int i = 0; i < 4; i++)
            {
                pts(0, i) = cit->vertex(i)->point().x() - center.x();
                pts(1, i) = cit->vertex(i)->point().y() - center.y();
                pts(2, i) = cit->vertex(i)->point().z() - center.z();
            }

            EMatrix tensor_mat = m_ref * pts.transpose();
            JacobiSVD svdStruct = tensor_mat.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);

            auto got_vol = vol_map.find(cit->info());
            double vol = (got_vol == vol_map.end()) ? std::abs((m_tr->tetrahedron(cit)).volume()) : got_vol->second;

            double vol_ref = 0.11785;
            double scale = (vol / vol_ref);
            double axisScaling = std::pow(scale, 0.3333);
            EMatrix unit = EMatrix::Identity(3, 3) * axisScaling;
            EMatrix rot = svdStruct.matrixV() * unit * svdStruct.matrixU().transpose();
            
            EMatrix new_pos = rot * m_ref;
            double new_vol = tetrahedron_volume(new_pos);
            double coeff = new_vol * std::pow(radius, 3);

            for(int i = 0; i < 4; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    int local_ind = m * 12 + i * 3 + j;
                    A_triplets.emplace_back(local_ind, 3 * cit->vertex(i)->index() + j, coeff);
                    B[local_ind] = coeff * (new_pos(j, i) + center[j]);
                }
            }
            
            m++;
        }

        // Location function
        for(Finite_vertices_iterator vit = m_tr->finite_vertices_begin(); 
            vit != m_tr->finite_vertices_end(); 
            vit++)
        {
            for(int i = 0; i < 3; i++)
                X[vit->index() * 3 + i] = vit->point()[i];
        }

        // solve system
        A.setFromTriplets(A_triplets.begin(), A_triplets.end());

        A = A.unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; });
        B = B.unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; });
        
        ESMatrix mat_left = A.transpose() * A;
        EVector vec_right = A.transpose() * B;
        ESSolver solver;
        
        X = solver.compute(mat_left).solveWithGuess(vec_right, X);

        // update vertex position
        for(Finite_vertices_iterator vit = m_tr->finite_vertices_begin(); 
            vit != m_tr->finite_vertices_end(); 
            vit++)
        {
            Point new_pos(X[3 * vit->index()], X[3 * vit->index() + 1], X[3 * vit->index() + 2]);

            if(!CGAL::is_valid(new_pos))
            {
                CGAL_TRACE_STREAM << "Invalid new position!" << std::endl;
                continue;
            }

            m_tr->move(vit, new_pos);
        }

        return;
    }

    unsigned int find_sliver_indices(double angle_thresh, IntSet& sliver_indices)
    {
        double tan_ub = tan(angle_thresh / 180.0 * boost::math::constants::pi<double>());

        for(Finite_cells_iterator cit = m_tr->finite_cells_begin(); 
            cit != m_tr->finite_cells_end(); 
            cit++)
        {
            FT min_angle_tan = 1 / largest_cot(cit);
            
            if(min_angle_tan < tan_ub)
                sliver_indices.insert(cit->info());
        }

        return sliver_indices.size();
    }    

    void smooth_volume_field(int smooth_range, VolumeMap& vol_map)
    {
        LabelMap m_neighbors;
        int num_cells = m_tr->number_of_finite_cells(); 
        vol_map.reserve(num_cells);
        m_neighbors.reserve(num_cells);

        for(Finite_cells_iterator cit = m_tr->finite_cells_begin(); 
            cit != m_tr->finite_cells_end(); 
            cit++)
        {
            double volume = std::abs((m_tr->tetrahedron(cit)).volume());
            vol_map.insert({cit->info(), volume});
            m_neighbors.insert({cit->info(), std::set<int>()});

            std::vector<Cell_handle> nb_cells;
            for(int i = 0; i < 4; i++)
                m_tr->incident_cells(cit->vertex(i), std::back_inserter(nb_cells));
            
            for(auto& cell: nb_cells)
            {
                if(m_tr->is_infinite(cell))
                    continue;

                m_neighbors[cit->info()].insert(cell->info());
            }
        }

        for(int i = 0; i < smooth_range; i++)
        {
            VolumeMap new_vol_map;
            new_vol_map.reserve(num_cells);

            for(auto& elem: vol_map)
            {
                double volume = elem.second;

                auto neighbors = m_neighbors[elem.first];

                for(auto& neighbor: neighbors)
                    volume += vol_map[neighbor];
                
                volume = volume / (1 + neighbors.size());
                new_vol_map.insert({elem.first, volume});
            }

            vol_map.swap(new_vol_map);
            new_vol_map.clear();
        }
    }

private: // functions

    double tetrahedron_volume(EMatrix& tet)
    {
        Vector a(tet(0, 1) - tet(0, 0), tet(1, 1) - tet(1, 0), tet(2, 1) - tet(2, 0));
        Vector b(tet(0, 2) - tet(0, 0), tet(1, 2) - tet(1, 0), tet(2, 2) - tet(2, 0));
        Vector c(tet(0, 3) - tet(0, 0), tet(1, 3) - tet(1, 0), tet(2, 3) - tet(2, 0));

        Vector cross_ab = CGAL::cross_product(a, b);
        double volume = std::abs(cross_ab * c) / 6.;

        return volume;
    }

    FT largest_cot(Cell_handle cell)
    {
        FT max_cotan = -1e7;
        for(int i = 0; i < 3; i++)
            for (int j = i + 1; j < 4; j++)
            {
                double cotan = cotan_per_edge(cell, i, j);
                if(cotan > max_cotan) max_cotan = cotan;
            }
        return max_cotan;
    };

    FT cotan_per_edge(Cell_handle cell, int i, int j)
    {
        Vertex_handle vi = cell->vertex(i);
        Vertex_handle vj = cell->vertex(j);

        Point pi = vi->point();
        Point pj = vj->point();

        std::vector<Point> vpq;

        for(int i = 0; i < 4; i++)
            if(cell->vertex(i)->index() != vi->index() && cell->vertex(i)->index() != vj->index())
                vpq.push_back(cell->vertex(i)->point());

        Vector ni = CGAL::cross_product(pi - vpq[0], pi - vpq[1]);
        Vector nj = CGAL::cross_product(pj - vpq[0], pj - vpq[1]);

        ni = ni / std::sqrt(ni * ni);
        nj = nj / std::sqrt(nj * nj);

        Vector nij = CGAL::cross_product(ni, nj);
        FT cotan = (ni * nj) / std::sqrt(nij * nij);

        return cotan;
    }

}; // end of class AsapOptimizer

// wrapper
template<class Gt,
         class Tr>
void asap_optimization(boost::shared_ptr<Tr>& tr, double threshold, int range = 1)
{
    AsapOptimizer<Gt, Tr> optimizer(tr);

    std::unordered_set<int> sliver_indices;
    optimizer.find_sliver_indices(threshold, sliver_indices);

    std::unordered_map<int, double> vol_map;
    optimizer.smooth_volume_field(range, vol_map);

    optimizer.optimize_asap_linear_solver(sliver_indices, vol_map);

}

} //namespace CGAL

#endif