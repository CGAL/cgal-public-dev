#ifndef TRIANGLE_FIT_H
#define TRIANGLE_FIT_H

#include <CGAL/squared_distance_3.h> 
#include <CGAL/SCIP_mixed_integer_program_traits.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/enum.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <boost/property_map/property_map.hpp>

#include <unordered_set>
#include <map>

#include "types.h"
#include "candidate.h"
#include "qem.h"
#include "alpha_shape_mesh.h"




namespace qem {
class TriangleFit{
public: // template
    typedef std::pair<int, int>                                     IntPair;
    typedef std::vector<int>                                        IntList;
    typedef std::vector<IntPair>                                    IntPairList;
    typedef std::vector<IntList>                                    IntListList;
    typedef typename std::vector<QEM_metric>                        QemList;

    typedef std::pair<Point, std::size_t>                                               Point_with_index;
    typedef std::vector<Point_with_index>                                               PwiList;

    typedef std::set<int>                                                                   IntSet;
    typedef typename std::unordered_set<IntSet, HashIntSet>                                 IntSetSet; 
    typedef typename std::unordered_map<int, IntSet>                                        IntIntSetMap;
    typedef typename std::unordered_map<IntPair, IntList, HashPairIndex, EqualPairIndex>    IntPairIntListMap;

    typedef CGAL::SCIP_mixed_integer_program_traits<double>         MIP_Solver;
    typedef typename MIP_Solver::Variable			                Variable;
	typedef typename MIP_Solver::Linear_objective	                Linear_objective;
	typedef typename MIP_Solver::Linear_constraint	                Linear_constraint;

    typedef typename Kernel::Point_2			                    Point2;
    typedef CGAL::Polygon_2<Kernel>				                    Polygon;
    typedef CGAL::Surface_mesh<Point>			                    Polygon_mesh;
    typedef typename Polygon_mesh::Face_index		                Face_descriptor;
	typedef typename Polygon_mesh::Edge_index		                Edge_descriptor;
	typedef typename Polygon_mesh::Vertex_index		                Vertex_descriptor;
    typedef typename Polygon_mesh::Halfedge_index		            Halfedge_descriptor;

    typedef typename Polyhedron::HalfedgeDS                         HalfedgeDS;
    typedef typename Polyhedron::Facet_handle 	                    Facet_handle;
    typedef typename Polyhedron::Halfedge_around_facet_circulator   Halfedge_around_facet_circulator;

    typedef typename CGAL::AABB_face_graph_triangle_primitive<Polyhedron>    AABBPrimitive;
    typedef typename CGAL::AABB_traits<Kernel, AABBPrimitive>                AABBTraits;
    typedef typename CGAL::AABB_tree<AABBTraits>                             AABBTree;
    typedef typename AABBTree::Point_and_primitive_id                        AABBPoint_and_primitive_id;

    typedef typename boost::graph_traits<Polyhedron>::face_descriptor        BFace_descriptor;
    typedef typename std::map<BFace_descriptor, Vector>                      BFaceNormalMap;


private: // member
    PointList   m_points;
    QemList     m_qems;  // not yet used
    IntPairList m_edges;
    IntListList m_facets;
    Bbox        m_bbox;
    double      m_diag;
    Polyhedron  m_dual_mesh;
    IntSet      m_select_indices;

public: // function

    TriangleFit() {}

    ~TriangleFit() 
    {
        reset();
    } 

    void reset()
    {
        m_points.clear();
        m_edges.clear();
        m_facets.clear();
        m_qems.clear();
        m_dual_mesh.clear();
        m_select_indices.clear();

        m_bbox = Bbox();
        m_diag = 0.;
    }

    template <class IntPairSet>
    void initialize_adjacent_graph(PointList& dual_points, IntPairSet& dual_edges, Bbox& bbox, float diag)
    {
        reset();

        // copy points
        std::copy(dual_points.begin(), dual_points.end(), std::back_inserter(m_points)); 

        // copy edges
        for(auto &elem: dual_edges)
            m_edges.push_back(std::make_pair(std::min(elem.first, elem.second), std::max(elem.first, elem.second)));

        // init bbox
        m_bbox = CGAL::bbox_3(m_points.begin(), m_points.end());
        m_diag = std::sqrt(CGAL::squared_distance(Point(m_bbox.min(0), m_bbox.min(1), m_bbox.min(2)),
                                                  Point(m_bbox.max(0), m_bbox.max(1), m_bbox.max(2))));
    }
    void create_candidate_facets()
    {
        if(m_points.size() == 0 || m_edges.size() == 0)
            return;

        m_facets.clear();
        IntIntSetMap adjacent_verts;

        // init map
        for(int i = 0; i < m_points.size(); i++)
            adjacent_verts.insert({i, IntSet()});
        // fill map
        for(int i = 0; i < m_edges.size(); i++)
        {
            adjacent_verts[m_edges[i].first].insert(m_edges[i].second);
            adjacent_verts[m_edges[i].second].insert(m_edges[i].first);
        }

        // find unique triple sets and useless edges
        IntSetSet unique_facets;
        IntList   remove_indices;
        for(int i = 0; i < m_edges.size(); i++)
        {
            IntSet edge_union;
            IntSet nb_p1 = adjacent_verts[m_edges[i].first];
            IntSet nb_p2 = adjacent_verts[m_edges[i].second];
            std::set_intersection(nb_p1.begin(), nb_p1.end(), nb_p2.begin(), nb_p2.end(), std::inserter(edge_union, edge_union.begin()));

            for(auto &elem: edge_union)
            {
                IntSet face({m_edges[i].first, m_edges[i].second, elem});
                unique_facets.insert(face);
            }
        }

        // fill face map
        for(auto &facet_set: unique_facets)
        {
            IntList facet_list;
            for(auto &index: facet_set)
                facet_list.push_back(index);

            m_facets.push_back(facet_list);
        }
    }
    void update_adjacent_edges(std::vector<float>& adjacent_edges)
    {
        if(m_points.size() == 0 || m_edges.size() == 0)
            return;

        for(int i = 0; i < m_edges.size(); i++)
        {
            auto pair = m_edges[i];
            Point p1 = m_points[pair.first];
            Point p2 = m_points[pair.second];

            adjacent_edges.push_back((float) p1.x());
            adjacent_edges.push_back((float) p1.y());
            adjacent_edges.push_back((float) p1.z());

            adjacent_edges.push_back((float) p2.x());
            adjacent_edges.push_back((float) p2.y());
            adjacent_edges.push_back((float) p2.z());
        }
    }

    void update_candidate_facets(std::vector<float>& candidate_facets, std::vector<float>& candidate_normals)
    {
        std::cout << "Candidate facet size: " << m_facets.size() << std::endl;

        for(int i = 0; i < m_facets.size(); i++)
        {
            Vector normal = CGAL::normal(m_points[m_facets[i][0]],
                                         m_points[m_facets[i][1]],
                                         m_points[m_facets[i][2]]);

            for(int j = 0; j < 3; j++)
            {
                Point p = m_points[m_facets[i][j]];
                candidate_facets.push_back((float) p.x());
                candidate_facets.push_back((float) p.y());
                candidate_facets.push_back((float) p.z());
                candidate_normals.push_back((float) normal.x());
                candidate_normals.push_back((float) normal.y());
                candidate_normals.push_back((float) normal.z());
            }
        }
    }

    void assemble_edge_face_adjacency(IntListList& edges_adj_facets_list)
    {
        IntPairIntListMap edge_face_map;

        for(int i = 0; i < m_facets.size(); i++)
        {
            IntList facet = m_facets[i];

            for(int j = 0; j < 3; j++)
            {
                IntPair edge_pair = std::make_pair(facet[(j + 1) % 3], facet[(j + 2) % 3]);

                if(edge_face_map.find(edge_pair) == edge_face_map.end())
                    edge_face_map.insert({edge_pair, IntList()});

                edge_face_map[edge_pair].push_back(i);
            }
        }

        for(int i = 0; i < m_edges.size(); i++)
        {
            IntPair edge_pair = std::make_pair(m_edges[i].first, m_edges[i].second);
            
            if(edge_face_map.find(edge_pair) == edge_face_map.end())
            {
                std::cout << "Single edge!" << std::endl;
                edges_adj_facets_list.push_back(IntList());
            }
            else
                edges_adj_facets_list.push_back(edge_face_map[edge_pair]);
        }

        /* // check
        for(int i = 0; i < edges_adj_facets_list.size(); i++)
        {
            std::cout << "Edge: " << m_edges[i].first << " " << m_edges[i].second << std::endl;
            IntList facet_indices = edges_adj_facets_list[i];

            for(int j = 0; j < facet_indices.size(); j++)
            {
                std::cout << "Adjacent facet " << facet_indices[j] << ": " << m_facets[facet_indices[j]][0] << " " << m_facets[facet_indices[j]][1] << " " << m_facets[facet_indices[j]][2] << std::endl;
            }
        } */
    }

    int assemble_support_points_for_faces(PointList& data_points, IntListList& face_point_map, std::vector<double>& face_confidence_map, double thresh)
    {
        int count = 0;
        std::cout << "Distance threshold for support points: " << thresh << std::endl;

        for(int i  = 0; i < m_facets.size(); i++)
        {
            face_point_map.push_back(IntList());
            face_confidence_map.push_back(0.);
        }

        for(int i = 0; i < data_points.size(); i++)
        {
            for(int j = 0; j < m_facets.size(); j++)
            {
                Triangle triangle(m_points[m_facets[j][0]], m_points[m_facets[j][1]], m_points[m_facets[j][2]]);
                double dist = CGAL::squared_distance(triangle, data_points[i]);

                if(dist < thresh)
                {
                    // record point index
                    face_point_map[j].push_back(i);

                    // compute confidence
                    double confidence = (1. - dist / thresh) * 1.; // 1. -> a confidence between 0 and 1, 0 less confidence, 1 best confidence
                    face_confidence_map[j] += confidence;

                    count++;
                }
            }
        }

        return count;
    }
    void nonmanifold_reconstruct(PointList& data_points, double data_spacing, double dist_ratio, double fitting, double coverage, double complexity)
    {
        if(m_points.size() == 0 || m_edges.size() == 0 || m_facets.size() < 4)
        { 
            std::cout << "Candidate not available! Reconstruction stops!" << std::endl;
            return;
        }

        int num_edges = m_edges.size();
        int num_facets = m_facets.size();
        std::cout << "Number of edges: " << num_edges << ", number of facets: " << num_facets << std::endl;

        // initialize adjacent map
        IntListList edges_adj_facets_list;
        assemble_edge_face_adjacency(edges_adj_facets_list);

        // initialize support points and assemble confidence map
        IntListList face_point_map;
        DoubleList face_confidence_map;
        double dist_thresh = m_diag * dist_ratio;
        dist_thresh = dist_thresh * dist_thresh;
        int num_support_points = assemble_support_points_for_faces(data_points, face_point_map, face_confidence_map, dist_thresh);
        std::cout << "Number of supported points: " << num_support_points << std::endl;
        std::cout << "Face point map size: " << face_point_map.size() << std::endl;
        std::cout << "Face confidence map size: " << face_confidence_map.size() << std::endl;
        
        // assemble coverage map
        DoubleList face_coverage_map;
        assemble_coverage_map(data_points, face_point_map, face_coverage_map, data_spacing);
        std::cout << "Face coverage map size: " << face_coverage_map.size() << std::endl;
      
        // Binary variables:
		// x[0] ... x[num_facets - 1] : binary labels of all the input faces
		// x[num_facets] ... x[num_facets + num_edges - 1] : binary labels of all the intersecting edges (remain or not)
        MIP_Solver solver;

        int total_variables = num_facets + num_edges;
		const std::vector<Variable*>& variables = solver.create_variables(total_variables);

        // init variables and objective
		for (std::size_t i = 0; i < total_variables; ++i) {
			Variable* v = variables[i];
			v->set_variable_type(Variable::BINARY);
		}

        Linear_objective* objective = solver.create_objective(Linear_objective::MINIMIZE);
        
        // init weights
        double dx = m_bbox.xmax() - m_bbox.xmin();
		double dy = m_bbox.ymax() - m_bbox.ymin();
		double dz = m_bbox.zmax() - m_bbox.zmin();
		double box_area = 2. * (dx * dy + dy * dz + dz * dx);
        double coeff_data_fitting = fitting;
		double coeff_coverage = num_support_points * coverage / num_facets;

        // add coverage and fitting terms
        for(int i = 0; i < m_facets.size(); i++) {
			double confidence = -coeff_data_fitting * std::abs(face_confidence_map[i]);
            double area_ratio = -coeff_coverage * std::abs(face_coverage_map[i]);

            //std::cout << "Confidence: " << confidence << ", area ratio: " << area_ratio << ", number of supported points: " << face_point_map[i].size() << std::endl;

			objective->add_coefficient(variables[i], confidence);
			objective->add_coefficient(variables[i], area_ratio);
		}



        // Add soft constraints: the number of faces associated with an edge must be either 2 or 0
		for (int i = 0; i < edges_adj_facets_list.size(); i++) 
        {
            IntList neighbors = edges_adj_facets_list[i];
            //Linear_constraint* c = solver.create_constraint(-Linear_constraint::infinity(), 0);

			for (int j = 0; j < neighbors.size(); j++) 
            {
                objective->add_coefficient(variables[neighbors[j]], complexity);
                //c->add_coefficient(variables[neighbors[j]], -1.0);
            }
                
			if (neighbors.size() > 1) 
            {
                objective->add_coefficient(variables[num_facets + i], -2.0 * complexity);  
                //c->add_coefficient(variables[num_facets + i], 1.0);  
            }
		}


        if(solver.solve()) 
        {
			m_dual_mesh.clear();
            m_select_indices.clear();

			const DoubleList& X = solver.solution();
            IntListList selected_facets;
            int count = 0;

            for(int i = 0; i < m_facets.size(); i++)
            {
                if(static_cast<int>(std::round(X[i])) == 1)
                {
                    m_select_indices.insert(i);
                    selected_facets.push_back(m_facets[i]);
                    count++;
                }
            }

            std::cout << "Point size: " << m_points.size() << ", selected facet size: " << count << std::endl;

            CGAL::Polygon_mesh_processing::orient_polygon_soup(m_points, selected_facets);
            CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(m_points, selected_facets, m_dual_mesh);

            std::cout << "Generated mesh with " << m_dual_mesh.size_of_vertices() << " points and " << m_dual_mesh.size_of_facets() << " facets!" << std::endl;
		}
		else 
        {
			std::cout << "solving the binary program failed!" << std::endl;
			return;
		}
    }
    void assemble_coverage_map(PointList& data_points, IntListList& face_point_map, std::vector<double>& face_coverage_map, double spacing)
    {
        for(int i = 0; i < m_facets.size(); i++)
        {
            if(face_point_map[i].size() < 3)
            {
                face_coverage_map.push_back(0.);
                continue;
            }

            IntList point_indices = face_point_map[i];
            PointList points;

            for(int j = 0; j < point_indices.size(); j++)
                points.push_back(data_points[point_indices[j]]);

            assert(points.size() == point_indices.size());

            Plane supporting_plane(m_points[m_facets[i][0]], m_points[m_facets[i][1]], m_points[m_facets[i][2]]);
            Point center = m_points[m_facets[i][0]];
            Vector base_1 = supporting_plane.base1();
            Vector base_2 = supporting_plane.base2();
            base_1 = base_1 / std::sqrt(base_1.squared_length());
            base_2 = base_2 / std::sqrt(base_2.squared_length());
            Alpha_shape_mesh alpha_mesh(points.begin(), points.end(), center, base_1, base_2);

            double alpha_radius = spacing * 5.;
            double covered_area = 0.;
            Polygon_mesh covering_mesh;

            if(alpha_mesh.extract_mesh(alpha_radius * alpha_radius, covering_mesh)) {
                // We cannot use the area of the 3D faces, because the alpha shape mesh is not perfectly planar
                const typename Polygon_mesh::template Property_map<Vertex_descriptor, Point>& coords = covering_mesh.points();
                
                for(auto face : covering_mesh.faces()) {    // We have to use the projected version
                    Polygon plg;                            // the projection of the face onto it supporting plane
                    CGAL::Halfedge_around_face_circulator<Polygon_mesh> cir(covering_mesh.halfedge(face), covering_mesh), done(cir);
                    do {
                        Halfedge_descriptor hd = *cir;
                        Vertex_descriptor vd = covering_mesh.target(hd);
                        const Point&  p = coords[vd];
                        Point2 q = alpha_mesh.project_point_on_plane(p);
                        plg.push_back(q);
                        ++cir;
                    } while (cir != done);
                    covered_area += std::abs(plg.area());
                }
            }
            //else
            //    std::cout << "Alpha extract mesh fails" << std::endl;

            double area = std::sqrt(CGAL::squared_area(m_points[m_facets[i][0]], m_points[m_facets[i][1]], m_points[m_facets[i][2]]));
            double area_ratio = std::min(covered_area, area) / area;

            //std::cout << "Number of supported points: " << points.size() << std::endl;
            //std::cout << "Cover mesh size: " << covering_mesh.number_of_vertices() << std::endl;
            //std::cout << "Cover area: " << covered_area << std::endl;
            //std::cout << "Area: " << area << std::endl;

            face_coverage_map.push_back(area_ratio);
        }
    }

    double compute_normal_orientation(int facet_ind, int edge_ind)
    {
        int p1 = m_edges[edge_ind].first;
        int p2 = m_edges[edge_ind].second;
        IntList facet = m_facets[facet_ind];

        //std::cout << "p1 p2: " << p1 << " " << p2 << std::endl;

        int local_ind = -1;

        for(int i = 0; i < 3; i++)
        {
            //std::cout << "Facet " << i << ": " << facet[i] << std::endl;

            if(facet[i] == p1)
            {
                local_ind = i;
                break;
            }
        }

        if(local_ind == -1)
        {
            std::cout << "Wrong local index!!!!!!" << std::endl;
            return 0.;
        }

        if(facet[(local_ind + 1) % 3] == p2)
            return 1.;
        
        return -1.;
    }



    void reconstruct(PointList& data_points, double data_spacing, double dist_ratio, double fitting, double coverage, double complexity)
    {
        if(m_points.size() == 0 || m_edges.size() == 0 || m_facets.size() < 4)
        { 
            std::cout << "Candidate not available! Reconstruction stops!" << std::endl;
            return;
        }

        int num_edges = m_edges.size();
        int num_facets = m_facets.size();
        std::cout << "Number of edges: " << num_edges << ", number of facets: " << num_facets << std::endl;

        // initialize adjacent map
        IntListList edges_adj_facets_list;
        assemble_edge_face_adjacency(edges_adj_facets_list);

        // initialize support points and assemble confidence map
        IntListList face_point_map;
        std::vector<double> face_confidence_map;
        double dist_thresh = m_diag * dist_ratio;
        dist_thresh = dist_thresh * dist_thresh;
        int num_support_points = assemble_support_points_for_faces(data_points, face_point_map, face_confidence_map, dist_thresh);
        std::cout << "Number of supported points: " << num_support_points << std::endl;
        std::cout << "Face point map size: " << face_point_map.size() << std::endl;
        std::cout << "Face confidence map size: " << face_confidence_map.size() << std::endl;
        
        // assemble coverage map
        std::vector<double> face_coverage_map;
        assemble_coverage_map(data_points, face_point_map, face_coverage_map, data_spacing);
        std::cout << "Face coverage map size: " << face_coverage_map.size() << std::endl;
      
        // Binary variables:
		// x[0] ... x[num_facets - 1] : binary labels of all the input faces
		// x[num_facets] ... x[num_facets + num_edges - 1] : binary labels of all the intersecting edges (remain or not)
        MIP_Solver solver;

        int total_variables = num_facets + num_edges;
		const std::vector<Variable*>& variables = solver.create_variables(total_variables);

        // init variables and objective
		for (std::size_t i = 0; i < total_variables; ++i) {
			Variable* v = variables[i];
			v->set_variable_type(Variable::BINARY);
		}

        Linear_objective* objective = solver.create_objective(Linear_objective::MINIMIZE);
        
        // init weights
        double dx = m_bbox.xmax() - m_bbox.xmin();
		double dy = m_bbox.ymax() - m_bbox.ymin();
		double dz = m_bbox.zmax() - m_bbox.zmin();
		double box_area = 2. * (dx * dy + dy * dz + dz * dx);
        double coeff_data_fitting = fitting;
		double coeff_coverage = num_support_points * coverage / num_facets;
        
        // add coverage and fitting terms
        for(int i = 0; i < m_facets.size(); i++) {
			double confidence = -coeff_data_fitting * face_confidence_map[i];
            double area_ratio = -coeff_coverage * face_coverage_map[i];

            //std::cout << "Confidence: " << confidence << ", area ratio: " << area_ratio << ", number of supported points: " << face_point_map[i].size() << std::endl;

			objective->add_coefficient(variables[i], confidence);
			objective->add_coefficient(variables[i], area_ratio);
		}

        // Add constraints: the number of faces associated with an edge must be either 2 or 0
		for (int i = 0; i < edges_adj_facets_list.size(); i++) 
        {
			Linear_constraint* c = solver.create_constraint(0.0, 0.0);
            IntList neighbors = edges_adj_facets_list[i];
			for (int j = 0; j < neighbors.size(); j++) 
                c->add_coefficient(variables[neighbors[j]], 1.0);

			if (neighbors.size() > 1) 
                c->add_coefficient(variables[num_facets + i], -2.0);  
		}
        if(solver.solve()) 
        {
			m_dual_mesh.clear();
            m_select_indices.clear();

			const std::vector<double>& X = solver.solution();
            IntListList selected_facets;
            int count = 0;

            for(int i = 0; i < m_facets.size(); i++)
            {
                if(static_cast<int>(std::round(X[i])) == 1)
                {
                    m_select_indices.insert(i);
                    selected_facets.push_back(m_facets[i]);
                    count++;
                }
            }

            std::cout << "Point size: " << m_points.size() << ", selected facet size: " << count << std::endl;

            CGAL::Polygon_mesh_processing::orient_polygon_soup(m_points, selected_facets);
            CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(m_points, selected_facets, m_dual_mesh);

            std::cout << "Generated mesh with " << m_dual_mesh.size_of_vertices() << " points and " << m_dual_mesh.size_of_facets() << " facets!" << std::endl;
		    save_trianglefit_mesh("test.off");
            save_trianglefit_soup("test_soup.off");
            save_candidate_edge_ply("testedge.ply");
        }
		else 
        {
			std::cout << "solving the binary program failed!" << std::endl;
			return;
		}
    }
    int count_non_manifold_edges(IntListList& edges_adj_facets_list)
    {
        int total_count = 0;

        for(int i = 0; i < edges_adj_facets_list.size(); i++)
        {
            IntList nb_faces = edges_adj_facets_list[i];
            int count_face = 0;

            for(int j = 0; j < nb_faces.size(); j++)
            {
                if(m_select_indices.find(nb_faces[j]) != m_select_indices.end())
                    count_face++;
            }

            if(count_face > 0 && count_face != 2)
                total_count++;
        }

        return total_count;
    }

    void project_input_normal(PwiList& data_points, std::vector<Vector> & data_normals)
    {
        if(m_dual_mesh.empty())
            return;

        BFaceNormalMap dual_mesh_normals;
        CGAL::Polygon_mesh_processing::compute_face_normals(m_dual_mesh, boost::make_assoc_property_map(dual_mesh_normals));

        data_normals.clear();
        PwiList new_data_points;
        new_data_points.reserve(data_points.size());

        AABBTree aabb_tree(faces(m_dual_mesh).first, faces(m_dual_mesh).second, m_dual_mesh);
        aabb_tree.accelerate_distance_queries();

        for(int i = 0; i < data_points.size(); i++)
        {
            AABBPoint_and_primitive_id pp = aabb_tree.closest_point_and_primitive(data_points[i].first);
            auto f = pp.second;
            data_normals.push_back(dual_mesh_normals[f]);
            new_data_points.push_back(std::make_pair(pp.first, data_points[i].second));
        }

        data_points.clear();
        std::copy(new_data_points.begin(), new_data_points.end(), std::back_inserter(data_points));

        std::cout << "Update data points and normals!" << std::endl;
    }

    void update_fit_surface(std::vector<float>& fit_facets, std::vector<float>& fit_normals)
    {
        if(m_dual_mesh.empty())
            return;
        
        fit_facets.reserve(m_dual_mesh.size_of_facets() * 9);
        fit_normals.reserve(m_dual_mesh.size_of_facets() * 9);

        for(Facet_handle fd: faces(m_dual_mesh)) 
        {
            Vector normal = CGAL::Polygon_mesh_processing::compute_face_normal(fd, m_dual_mesh);
            Halfedge_around_facet_circulator vit = fd->facet_begin();
            do {
                Point p = vit->vertex()->point();
                fit_facets.push_back((float)p.x());
                fit_facets.push_back((float)p.y());
                fit_facets.push_back((float)p.z());
                fit_normals.push_back((float)normal.x());
                fit_normals.push_back((float)normal.y());
                fit_normals.push_back((float)normal.z());
            } while ( ++vit != fd->facet_begin());
        }
    }

    void update_fit_soup(std::vector<float>& fit_soup_facets, std::vector<float>& fit_soup_normals)
    {
        if(m_facets.size() == 0 || m_select_indices.size() == 0)
            return;
        
        fit_soup_facets.reserve(m_select_indices.size() * 9);
        fit_soup_normals.reserve(m_select_indices.size() * 9);

        for(int i = 0; i < m_facets.size(); i++) 
        {
            if(m_select_indices.find(i) == m_select_indices.end())
                continue;

            IntList face = m_facets[i];
            Vector normal = CGAL::normal(m_points[face[0]], m_points[face[1]], m_points[face[2]]);

            for(int j = 0; j < 3; j++)
            {
                Point p = m_points[face[j]];
                fit_soup_facets.push_back((float)p.x());
                fit_soup_facets.push_back((float)p.y());
                fit_soup_facets.push_back((float)p.z());
                fit_soup_normals.push_back((float)normal.x());
                fit_soup_normals.push_back((float)normal.y());
                fit_soup_normals.push_back((float)normal.z());
            }
        }
    }

    void save_trianglefit_mesh(std::string filename)
    {
        if(m_dual_mesh.empty())
            return;

        std::ofstream mesh_file;
        mesh_file.open(filename);
        CGAL::write_off(mesh_file, m_dual_mesh);
        mesh_file.close();
    }

    void save_trianglefit_soup(std::string filename)
    {
        if(m_facets.size() == 0)
            return;

        std::ofstream soup_file;
        soup_file.open(filename);

        soup_file << "OFF"  << std::endl;
        soup_file << m_points.size() << " " << m_facets.size() << " 0" << std::endl;

        for(int i = 0; i < m_points.size(); i++)
            soup_file << m_points[i].x() << " " << m_points[i].y() << " " << m_points[i].z() << std::endl;

        for(int i = 0; i < m_facets.size(); i++)
            soup_file << "3 " << m_facets[i][0] << " " << m_facets[i][1] << " " << m_facets[i][2] << std::endl;

        soup_file.close();
    }

    void save_candidate_edge_ply(std::string filename)
    {
        std::ofstream edge_file;
        edge_file.open(filename);

        edge_file << "ply\n"
                  << "format ascii 1.0\n"
                  << "element vertex " << m_points.size() << "\n"
                  << "property float x\n"
                  << "property float y\n"
                  << "property float z\n"
                  << "element edge " << m_edges.size() << "\n"
                  << "property int32 vertex1\n"
                  << "property int32 vertex2\n"
                  << "end_header\n";

        for(int i = 0; i < m_points.size(); i++)
            edge_file << m_points[i].x() << " " << m_points[i].y() << " " << m_points[i].z() << std::endl;
        
        for(int i = 0; i < m_edges.size(); i++)
            edge_file << m_edges[i].first << " " << m_edges[i].second << std::endl;

        edge_file.close();
    }
        const Polyhedron& get_mesh()
    {
        return m_dual_mesh;
    }

}; // end of class TriangleFit

} // end of namespace qem

#endif
