#include "propagation.h"
#include "geometry.h"
#include <ctime>
#include <queue>
#include "defs.h"
#include "trace.h"
#include "means.h"
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <fstream>
#include "svg.h"



void Propagation::propagate(Kinetic_Model *model)
{
    // Input data, constant
    Parameters* params = model->params;
    Matrix<uchar> & I = model->I;
    vector<Segment *> & segments = model->segments;

    // Data to be accessed and modified
    vector<Ray *> & rays = model->rays;
    IndexedEvent* & schedule = model->schedule;
    Partition* & graph = model->graph;

    // Deletes data from previous execution
    Ray::clear_rays(rays);
    IndexedEvent::clear_schedule(schedule);
    if (graph != NULL) delete graph;
    graph = new Partition(I.rows, I.cols, params);

	clock_t t_begin, t_end;
	// All segments may have been reoriented and replaced to favor geometrical relationship between them.
	// However, this might have been created groups of colinear, even overlapping segments, which may be
	// troublesome when propagating the rays stemed from them. That's why we are going to identify pairs
	// of colinear segments, and merge overlapping segments.
	t_begin = clock();

	if (!params->rega_regp_enabled) {
		Geometry::build_tree_when_regularization_disabled(segments, model->tree);
	}
	Geometry::merge_overlapped_segments(model->tree);
	Geometry::disable_segments_outside_boundaries(segments, model->I.rows, model->I.cols);

	if (params->lsd_create_additional) {
		for (uint i = 0 ; i < segments.size() ; i++) {

			if (segments[i]->is_artificial && !segments[i]->is_disabled) {
				Point2d bar_i = segments[i]->finalBarycenter;
				for (uint j = 0 ; j < segments.size() ; j++) {

					Point2d end1_j = segments[j]->finalEnd1;
					Point2d end2_j = segments[j]->finalEnd2;
					if (cv::norm(end1_j - bar_i) < 1e-9 || cv::norm(end2_j - bar_i) < 1e-9)  {
						segments[i]->disable();
						break;
					}
				}
			}
		}
	}

	// Now we build two rays per segment
    Ray::build_rays(segments, rays, params->prop_ttl);

	size_t nb_rays = size_t(rays.size());

	// Gets the image
    Size2i size = Size2i(I.cols, I.rows);

	// Defines two lists of events : the first one is a raw, complete list of possible events,
	// and by pruning it we obtain a new list that only contains the events happening for sure
	IndexedEvent* first_schedule = NULL;
	IndexedEvent* simplified_schedule = NULL;
	IndexedEvent* last_event_of_simplified_schedule = NULL;

	// Here, we define some tables that are of particular interest to handle the problem in an efficient way
	// Indeed, collisions between all the rays define a sparse matrix of events, that can be, on the whole,
	// represented by a quadruplet (i, j, t_i, t_j). It contains many more informations and in particular,
	// pointers to other events. Here, we define tables of pointers that give us a direct access to the first
	// elements of the matrices along each row and each column.
	// We also define a table that contains all the elements involving a propagating ray and an image boundary
	//t_begin = clock();
	IndexedEvent** intersectants = new IndexedEvent*[nb_rays];
	for (unsigned int i = 0; i < nb_rays; i++) {
		intersectants[i] = NULL;
	}
	IndexedEvent** intersected = new IndexedEvent*[nb_rays];
	for (unsigned int i = 0; i < nb_rays; i++) {
		intersected[i] = NULL;
	}
	IndexedEvent** events_at_boundaries = new IndexedEvent*[nb_rays];
	for (unsigned int i = 0; i < nb_rays; i++) {
		events_at_boundaries[i] = NULL;
	}

	// Over a first phase, we are going to compute all the possible intersections between the different objects
	// of the image. Later, we will prune the list of possible intersections and build a graph that actually
	// models the partitioning of the image (every vertex or so correspond to an intersection, and every edge
	// represents a position of the path actually followed by a given ray). Here, we define some structures
	// capable of handling such information.
	list<Vertex *> outer_vertices, inner_vertices;
	vector<list<Outer_Edge *> > outer_edges;
	vector<list<Inner_Edge *> > inner_edges;
	vector<Vertex *> initial_vertices_of_rays;

	// Defines the first vertices and edges of the graph :
	// Four vertices corresponding to the corners of the image, N vertices for the centers of the N segments
	// We also add edges that correspond to the boundaries of the image
    graph->init_edges(rays, outer_vertices, inner_vertices, outer_edges, inner_edges, initial_vertices_of_rays);

	// The computation of events is performed by a loop : while there remain active rays, we compute all the possible
	// intersections happening within a certain range of time [t_0, t_1] and update the graph of intersections.
	vector<double> maximal_sizes;
	vector<IndexedEvent *> events_colinear;
    schedule_events_at_boundaries(model, I.rows, I.cols, maximal_sizes, events_at_boundaries);
    schedule_events_between_colinear_segments(model, events_colinear);
	int active_rays = nb_rays;
	int iteration = 0;
    double range = params->prop_range;
	double t_0 = -FLT_MAX, t_1 = range;

	while (active_rays > 0) {
		trace(model->params->verbose_level, 10, "Active rays : " + std::to_string(active_rays));

		// Computes all events happening between t_0 and t_1
        schedule_events(model, &first_schedule, intersectants, intersected, events_at_boundaries, maximal_sizes, events_colinear, t_0, t_1);

		// Simplifies this list of events
        prune_events_and_build_edges(model, size, &first_schedule, &simplified_schedule, &last_event_of_simplified_schedule,
            intersectants, intersected, events_at_boundaries, graph->quadtree, t_0, t_1,
			outer_vertices, inner_vertices, outer_edges, inner_edges, initial_vertices_of_rays, active_rays);

		// Updates the temporal range
        ++iteration;
		t_0 = iteration * range;
		if (active_rays > 300) {
			t_1 = (iteration + 1) * range;
		} else {
			// If there are only a few rays left (less than 5%) we can afford to accept any kind of events
			t_1 = FLT_MAX;
		}
	}

	delete[] intersectants;
	delete[] intersected;
	delete[] events_at_boundaries;
	initial_vertices_of_rays.clear();

	t_end = clock();
	//trace(model->params->verbose_level, 5, "** Scheduling events and building graph : " + std::to_string(float(t_end - t_begin) / CLOCKS_PER_SEC) + " s.");

	schedule = simplified_schedule;
	//write_schedule(params, schedule);

	/*Matrix<uchar> B;
	draw_rays(B, sqrt(A.rows * A.rows + A.cols * A.cols));
	std::string B_name = params.prefix + "_" + std::to_string(params.block_index) + "_" + std::to_string(params.prt_policy) + "_intersections.tiff";
	B.write_uchar(B_name); B.release();
	
	Matrix<uchar> C;
	graph->draw_intermediate_graph(C, outer_vertices, inner_vertices, outer_edges, inner_edges);
	std::string C_name = params.prefix + "_" + std::to_string(params.block_index) + "_" + std::to_string(params.prt_policy) + "_intermediate_edges.tiff";
	C.write_uchar(C_name); C.release();*/
	
	// Last step of the algorithm : now that all the edges are defined, we loop on these objects in order
	// to get the list of pixels that compose the facets of the algorithm. But to this end we need to define
	// one last auxiliary structure

	graph->merge_containers(outer_vertices, inner_vertices, outer_edges, inner_edges);

	//t_begin = clock();
	graph->build_faces(size);
	//t_end = clock();
	//trace(model->params->verbose_level, 5, "** Built " + std::to_string(graph->faces.size()) + " facets in " + std::to_string(float(t_end - t_begin) / CLOCKS_PER_SEC) + " s.");

	if (params->merge_enabled) {
		graph->merge_thin_facets(rays, model->params->verbose_level);
	}

	graph->seal_and_remove_bivalent_vertices();
	graph->reset_indices();

	t_end = clock();
	trace(model->params->verbose_level, 5, "** Built " + std::to_string(graph->faces.size()) + " facets in " + std::to_string(float(t_end - t_begin) / CLOCKS_PER_SEC) + " s.");

#if NOT_MEASURING_PERFORMANCES
    make_layer(model);
	print_histogram(model);
#endif

    // Cleans memory
	// Ray::clear_rays(rays);

    // And enables again segments whose state has been modified while calling the merge_overlapped segments
    Segment::enable(segments);

	// Saves the graph
	graph->save_graph_definition(params->path_output_directory, model->basename);
    //graph->save_liuyuns_input(params->path_output_directory, model->basename, segments);
#if NOT_MEASURING_PERFORMANCES
	model->partition_to_svg(params->path_output_directory);
	model->harlequin_to_svg(params->path_output_directory);
#endif
}


#if NOT_MEASURING_PERFORMANCES
void Propagation::print_histogram(Kinetic_Model *model)
{
	/*std::vector<double> hist_values(360, 0);
	int clusters_with_at_least_two_colinear_segments = 0;
	for (auto it_n1 = model->tree->parallel_segments.begin() ; it_n1 != model->tree->parallel_segments.end() ; it_n1++) {
		double theta = it_n1->first;
		int theta_bin = jclamp(0, int(4 * theta), 359);
		double value = 0;
		for (auto it_n2 = it_n1->second->colinear_segments.begin() ; it_n2 != it_n1->second->colinear_segments.end() ; it_n2++) {
			for (list<Segment *>::iterator it_s = it_n2->second->colinear_segments.begin() ; it_s != it_n2->second->colinear_segments.end() ; it_s++) {
				value += (*it_s)->length;
			}
			if (it_n2->second->colinear_segments.size() > 1) {
				clusters_with_at_least_two_colinear_segments += 1;
			}
		}
		hist_values[theta_bin] += value;
	}

	std::string filename_histogram = model->params->path_output_directory + '\\' + model->basename + "_histogram.csv";
	FILE* file_histogram = fopen(filename_histogram.c_str(), "w");
	if (file_histogram != NULL) {
		for (int i = 0 ; i < 360 ; i++) {
			fprintf(file_histogram, "%f, %f\n", 0.5 * i, hist_values[i]);
		}
		fclose(file_histogram);
	}

	std::string filename_values = model->params->path_output_directory + '\\' + model->basename + "_raw_values.csv";
	FILE* file_values = fopen(filename_values.c_str(), "w");
	if (file_values != NULL) {
		for (int i = 0 ; i < model->segments.size() ; i++) {
			double theta = atan2(model->segments[i]->finalDirection[1], model->segments[i]->finalDirection[0]);
			fprintf(file_values, "%f\n", (theta < 0 ? theta + PI : theta));
		}
		fclose(file_values);
	}

	std::string filename_clusters =  model->params->path_output_directory + '\\' + model->basename + "_clusters.txt";
	FILE* file_clusters = fopen(filename_clusters.c_str(), "w");
	if (file_clusters != NULL) {
		fprintf(file_clusters, "Clusters : %i\n", clusters_with_at_least_two_colinear_segments);
		fclose(file_clusters);
	}*/

	// Computes histograms
	std::vector<int> initial_histogram(90, 0);
	std::vector<int> final_histogram(90, 0);
	for (uint i = 0 ; i < model->segments.size() ; i++) {
		Segment* s_i = model->segments[i];
		double theta_i = atan2(s_i->direction[1], s_i->direction[0]);
		double theta_f = atan2(s_i->finalDirection[1], s_i->finalDirection[0]);
		if (theta_i < 0) theta_i += PI;
		if (theta_f < 0) theta_f += PI;
		int bin_i = jclamp(0, 89 * theta_i / PI, 90);
		int bin_f = jclamp(0, 89 * theta_f / PI, 90);
		initial_histogram[bin_i] += 1;
		final_histogram[bin_f] += 1;
	}

	// Gets maximal value
	int max_f = 0;
	for (uint i = 0 ; i < final_histogram.size() ; i++) {
		if (final_histogram[i] > max_f) max_f = final_histogram[i];
	}
	double log_max_f = log10(max_f);

	int h_height = 500;
	int h_width = 2 * h_height;

	// Prints normalized histograms
	std::string filename_rose_diagram = model->params->path_output_directory + '\\' + model->basename + "_rose_diagram_initial.svg";
	std::filebuf fb;

	fb.open(filename_rose_diagram, std::ios::out);
	std::ostream os(&fb);

	// Size of the image
	SVG::markup_header(os, h_height, h_width);

	// Prints circles for each log
	for (int i = int(ceil(log_max_f)); i >= 1; i--) {
		double r = i / (ceil(log_max_f));
		SVG::markup_point(os, h_height, h_height, r * h_height, "gray", 1, "white");
	}

	// Prints triangles for each bin
	for (int i = 0; i < initial_histogram.size() ; i++) {
		int bin_i = initial_histogram[i];
		double theta_i = i * PI / initial_histogram.size();
		double theta_j = (i + 1) * PI / initial_histogram.size();
		if (bin_i > 0) {
			double M = log10(bin_i) / ceil(log_max_f);
			list<Point2d> triangle_i;
			triangle_i.push_back(Point2d(h_height, h_height));
			triangle_i.push_back(Point2d(h_height + M * h_height * cos(theta_i), h_height - M * h_height * sin(theta_i)));
			triangle_i.push_back(Point2d(h_height + M * h_height * cos(theta_j), h_height - M * h_height * sin(theta_j)));
			SVG::markup_polygon(os, triangle_i, 255, 0, 0, 1);
		}
	}

	SVG::markup_footer(os);
	fb.close();

	filename_rose_diagram = model->params->path_output_directory + '\\' + model->basename + "_rose_diagram_final.svg";
	fb.open(filename_rose_diagram, std::ios::out);
	std::ostream os_f (&fb);

	// Size of the image
	SVG::markup_header(os_f, h_height, h_width);

	// Prints circles for each log
	for (int i = int(ceil(log_max_f)); i >= 1; i--) {
		double r = i / (ceil(log_max_f));
		SVG::markup_point(os_f, h_height, h_height, r * h_height, "gray", 1, "white");
	}

	// Prints triangles for each bin
	for (int i = 0; i < final_histogram.size() ; i++) {
		int bin_i = final_histogram[i];
		double theta_i = i * PI / final_histogram.size();
		double theta_j = (i + 1) * PI / final_histogram.size();
		if (bin_i > 0) {
			double M = log10(bin_i) / ceil(log_max_f);
			list<Point2d> triangle_i;
			triangle_i.push_back(Point2d(h_height, h_height));
			triangle_i.push_back(Point2d(h_height + M * h_height * cos(theta_i), h_height - M * h_height * sin(theta_i)));
			triangle_i.push_back(Point2d(h_height + M * h_height * cos(theta_j), h_height - M * h_height * sin(theta_j)));
			SVG::markup_polygon(os_f, triangle_i, 255, 0, 0, 1);
		}
	}

	SVG::markup_footer(os_f);
	fb.close();
}
#endif

void Propagation::build_rtree(vector<Segment *> & segments, double & D, Boost_RTree & rtree_segments)
{
    double x_min, x_max, y_min, y_max;
    for (uint i = 0 ; i < segments.size() ; i++) {

        // For each segment, we create a box that represents the bounding box of this segment
        Segment* s_i = segments[i];
        if (!s_i->is_disabled) {
            if (s_i->finalEnd1.x < s_i->finalEnd2.x) {
                x_min = s_i->finalEnd1.x; x_max = s_i->finalEnd2.x;
            } else {
                x_min = s_i->finalEnd2.x; x_max = s_i->finalEnd1.x;
            }
            if (s_i->finalEnd1.y < s_i->finalEnd2.y) {
                y_min = s_i->finalEnd1.y; y_max = s_i->finalEnd2.y;
            } else {
                y_min = s_i->finalEnd2.y; y_max = s_i->finalEnd1.y;
            }

            // We insert this box in the r-tree
            Boost_Box s_box(Boost_Point(x_min - D, y_min - D), Boost_Point(x_max + D, y_max + D));
            rtree_segments.insert(std::make_pair(s_box, i));
        }
    }
}



void Propagation::search_neighborhood(Segment* s, double & D, Boost_RTree & rtree_segments, vector<Segment *> & segments, list<Segment *> & neighborhood)
{
    neighborhood.clear();

    // We search for segments whose bounding boxes' edges are located at a distance lower than D from the bounding box of s
    double x_min, x_max, y_min, y_max;
    vector<Boost_Value> possible_neighbors;

    if (s->finalEnd1.x < s->finalEnd2.x) {
        x_min = s->finalEnd1.x; x_max = s->finalEnd2.x;
    } else {
        x_min = s->finalEnd2.x; x_max = s->finalEnd1.x;
    }
    if (s->finalEnd1.y < s->finalEnd2.y) {
        y_min = s->finalEnd1.y; y_max = s->finalEnd2.y;
    } else {
        y_min = s->finalEnd2.y; y_max = s->finalEnd1.y;
    }
    Boost_Box query(Boost_Point(x_min - D, y_min - D), Boost_Point(x_max + D, y_max + D));
    rtree_segments.query(bgi::intersects(query), std::back_inserter(possible_neighbors));

    for (uint i = 0 ; i < possible_neighbors.size() ; i++) {
        Segment* s_i = segments[possible_neighbors[i].second];
        neighborhood.push_back(s_i);
    }
}




void Propagation::schedule_events_at_boundaries(Kinetic_Model* model, int rows, int cols, vector<double> & maximal_sizes, IndexedEvent** events_at_boundaries)
{
    vector<Ray *> & rays = model->rays;

    uint nb_rays = uint(rays.size());
	maximal_sizes = vector<double>(nb_rays, 0.0);
	Image_Boundary boundary;
	double t_i = FLT_MAX, t_j = FLT_MAX;

    for (uint i = 0; i < nb_rays; i++) {
		Geometry::intersection_boundary(rays[i], rows, cols, t_i, boundary, t_j);
		maximal_sizes[i] = t_i;
		events_at_boundaries[i] = new IndexedEvent(i, boundary, t_i, t_j);
	}
}




void Propagation::schedule_events_between_colinear_segments(Kinetic_Model* model, vector<IndexedEvent *> & events_colinear)
{
	Segment_Regularization_Tree* tree = model->tree;

	// Accesses to colinear segments
	for (map<double, Node_Parallel_Segments*>::iterator it_m1 = tree->parallel_segments.begin() ; it_m1 != tree->parallel_segments.end() ; it_m1++) {
		Node_Parallel_Segments* node_parallel = it_m1->second;
		for (map<double, Node_Colinear_Segments*>::iterator it_m2 = node_parallel->colinear_segments.begin() ; it_m2 != node_parallel->colinear_segments.end() ; it_m2++) {
			Node_Colinear_Segments* node_colinear = it_m2->second;

			if (node_colinear->colinear_segments.size() == 1) {
				continue;
			} else {
				list<Segment *>::iterator it_s1, it_s2;
				for (it_s1 = node_colinear->colinear_segments.begin() ; it_s1 != node_colinear->colinear_segments.end() ; it_s1++) {
					Segment* s1 = (*it_s1);
					if (s1->is_disabled) {
						continue;
					}
					it_s2 = it_s1;
					while (++it_s2 != node_colinear->colinear_segments.end()) {
						Segment* s2 = (*it_s2);
						if (s2->is_disabled) {
							continue;
						}

						bool loop = true;
						for (int i = 0; i < 2 && loop; i++) {
							Ray* r_i = (i == 0 ? s1->rays.first : s1->rays.second);
							for (int j = 0; j < 2 && loop; j++) {
								Ray* r_j = (j == 0 ? s2->rays.first : s2->rays.second);
								double t_i, t_j;
								bool exists = Geometry::intersection_colinear(r_i, r_j, t_i, t_j);
								if (exists) {
									// We choose to duplicate the event
									// Later, when events will be pruned, every time we will consider this kind of event we will
									// also has to take into account the "twin" event
									IndexedEvent* event_ij = new IndexedEvent(int(r_i->index), int(r_j->index), t_i, t_j, true);
									IndexedEvent* event_ji = new IndexedEvent(int(r_j->index), int(r_i->index), t_j, t_i, true);
									// Insert both events in the table
									events_colinear.push_back(event_ij);
									events_colinear.push_back(event_ji);
									loop = false;
								}
							}
						}
					}
				}
			}
		}
	}

	std::sort(events_colinear.begin(), events_colinear.end(), before_indexed_event);
}



void Propagation::schedule_events_between_non_colinear_segments(vector<Segment *> & segments, vector<Ray *> & rays, vector<double> & maximal_sizes, IndexedEvent** intersectants, IndexedEvent** intersected, vector<IndexedEvent *> & T, double & t_0, double & t_1)
{
    // We use a structure of sparse matrix to store indexed events : but we also know in which order the events are created
    // Indeed, each new element is added at the end of the current row, or at the end of the current column
    // That's why we use the following tables : to get a quick reference on the last element of each row and the last element of each column
    uint nb_rays = uint(rays.size());
    IndexedEvent** intersectant_tails = new IndexedEvent*[nb_rays];
    IndexedEvent** intersected_tails = new IndexedEvent*[nb_rays];
    for (unsigned int i = 0; i < nb_rays; i++) {
        intersectant_tails[i] = NULL;
        intersected_tails[i] = NULL;
    }

    // The algorithm to compute the intersection of two rays r_i and r_j takes the shape of a double loop on i and j,
    // where i and j represent the indices of the segments. We suppose that i is the intersectant ray.
    uint nb_segments = uint(segments.size());
    for (uint i = 0 ; i < nb_segments ; i++) {
        Segment* s_i = segments[i];

        // If the i-th segment has been discarded, or if both rays related to it have already been stopped, jump to the next segment
        if (s_i->is_disabled || (s_i->rays.first->has_been_stopped() && s_i->rays.second->has_been_stopped())) continue;
        Ray* r_i = s_i->rays.first;

        // For every segment likely to be intersected by the segment i-th
        for (uint j = 0 ; j < nb_segments ; j++) {
            if (i == j) continue;
			Segment* s_j = segments[j];

            // If segment j is disabled, we jump to the next segment
            if (s_j->is_disabled) continue;
            Ray* r_j = s_j->rays.first;

            uint ind_i = -1, ind_j = -1;
            double t_i = -FLT_MAX, t_j = -FLT_MAX;
            bool exists = Geometry::intersection(r_i, r_j, maximal_sizes[r_i->index], maximal_sizes[r_i->index + 1], ind_i, ind_j, t_i, t_j);
            if (exists) {
                uint index_intersectant, index_intersected;
                double t_intersectant, t_intersected;
                if (t_i >= t_j) {
                    index_intersectant = ind_i;
                    index_intersected = ind_j;
                    t_intersectant = t_i;
                    t_intersected = t_j;
                } else {
                    index_intersectant = ind_j;
                    index_intersected = ind_i;
                    t_intersectant = t_j;
                    t_intersected = t_i;
                }
                // If the intersectant ray actually belongs to the i-th segment and if it is still propagating
                if (rays[index_intersectant]->parent->index == i && !rays[index_intersectant]->has_been_stopped() && t_0 <= t_intersectant && t_intersectant < t_1) {
                    // There are two possibilites for an a priori valid intersection :
                    // - the intersected ray is still active
                    // - the intersected ray has been stopped but it is intersected in a point that actually exists (check time)
                    if (!rays[index_intersected]->has_been_stopped() || (rays[index_intersected]->has_been_stopped() && t_intersected <= rays[index_intersected]->t)) {
                        IndexedEvent* event = new IndexedEvent(index_intersectant, index_intersected, t_intersectant, t_intersected, false);
                        // Insert it in the table
                        insert_indexed_event(event, intersectants, intersected, intersectant_tails, intersected_tails);
                        // Inserts it in a list, which will be later sorted
                        T.push_back(event);
                    }
                }
            }
        }
    }
    delete[] intersectant_tails;
    delete[] intersected_tails;
}


void Propagation::schedule_events_between_non_colinear_segments_with_rtree(vector<Segment *> & segments, vector<Ray *> & rays, vector<double> & maximal_sizes, IndexedEvent** intersectants, IndexedEvent** intersected, vector<IndexedEvent *> & T, double &t_0, double & t_1)
{
    Boost_RTree rtree_segments;
    build_rtree(segments, t_1, rtree_segments);

    // Initializes tails
    uint nb_rays = uint(rays.size());
    IndexedEvent** intersectant_tails = new IndexedEvent*[nb_rays];
    IndexedEvent** intersected_tails = new IndexedEvent*[nb_rays];
    for (unsigned int i = 0; i < nb_rays; i++) {
        intersectant_tails[i] = NULL;
        intersected_tails[i] = NULL;
    }

    uint nb_segments = uint(segments.size());
    for (uint i = 0 ; i < nb_segments ; i++) {
        Segment* s_i = segments[i];

        if (s_i->is_disabled || (s_i->rays.first->has_been_stopped() && s_i->rays.second->has_been_stopped())) continue;
        Ray* r_i = s_i->rays.first;

        // Searches for all the segments which are likely to be intersected by s_i
        list<Segment *> neighbors;
        search_neighborhood(s_i, t_1, rtree_segments, segments, neighbors);
        vector<IndexedEvent*> local_events;

        // For every segment likely to be intersected by the i-th segment
        for (list<Segment *>::iterator it_s = neighbors.begin() ; it_s != neighbors.end() ; it_s++) {

            Segment* s_j = (*it_s);
            if (s_i == s_j) continue;
            Ray* r_j = s_j->rays.first;

            uint ind_i = -1, ind_j = -1;
            double t_i = -FLT_MAX, t_j = -FLT_MAX;
            bool exists = Geometry::intersection(r_i, r_j, maximal_sizes[r_i->index], maximal_sizes[r_i->index + 1], ind_i, ind_j, t_i, t_j);
            if (exists) {
                uint index_intersectant, index_intersected;
                double t_intersectant, t_intersected;
                if (t_i >= t_j) {
                    index_intersectant = ind_i;
                    index_intersected = ind_j;
                    t_intersectant = t_i;
                    t_intersected = t_j;
                } else {
                    index_intersectant = ind_j;
                    index_intersected = ind_i;
                    t_intersectant = t_j;
                    t_intersected = t_i;
                }

                if (rays[index_intersectant]->parent->index == i && !rays[index_intersectant]->has_been_stopped() && t_0 <= t_intersectant && t_intersectant < t_1) {
                    if (!rays[index_intersected]->has_been_stopped() || (rays[index_intersected]->has_been_stopped() && t_intersected <= rays[index_intersected]->t)) {
                        IndexedEvent* event = new IndexedEvent(index_intersectant, index_intersected, t_intersectant, t_intersected, false);
                        // We do not insert this event in the matrix of events for the moment
                        // Instead we put it in a vector, sorted on index_intersected
                        local_events.push_back(event);
                    }
                }
            }
        }

        // Sorts the table of events, and inserts elements in the chronological list of events, as well as in the matrix of events
        std::sort(local_events.begin(), local_events.end(), sorted_by_index_of_intersected_segment);
        for (uint l = 0 ; l < local_events.size() ; l++) {
            insert_indexed_event(local_events[l], intersectants, intersected, intersectant_tails, intersected_tails);
            T.push_back(local_events[l]);
        }
    }

    rtree_segments.clear();
    delete[] intersectant_tails;
    delete[] intersected_tails;
}



void Propagation::schedule_events(Kinetic_Model* model, IndexedEvent** schedule, IndexedEvent** intersectants, IndexedEvent** intersected, IndexedEvent** events_at_boundaries, vector<double> & maximal_sizes, vector<IndexedEvent *> & events_colinear, double t_0, double t_1)
{
    vector<Segment *> & segments = model->segments;
    vector<Ray *> & rays = model->rays;
	vector<IndexedEvent *> T = vector<IndexedEvent *>();
    uint nb_rays = uint(rays.size());

	// First step : selects the possible intersections of the rays with the boundaries
    for (uint i = 0; i < nb_rays; i++) {
		IndexedEvent* event = events_at_boundaries[i];
		if (!rays[i]->has_been_stopped() && event->t_intersectant >= t_0 && event->t_intersectant < t_1) {
			T.push_back(event);
		}
	}

	// Second, computes the intersections between all the non-colinear rays
    // Our algorithm takes the shape of a double loop on segments
    // To reduce the complexity we may use a structure of R-Tree
    if (t_1 < FLT_MAX) {
        schedule_events_between_non_colinear_segments_with_rtree(segments, rays, maximal_sizes, intersectants, intersected, T, t_0, t_1);
    } else {
        schedule_events_between_non_colinear_segments(segments, rays, maximal_sizes, intersectants, intersected, T, t_0, t_1);
    }

	// Third step : take into account all the intersections between colinear rays
	if (events_colinear.size() > 0) {
		vector<IndexedEvent *>::iterator it = events_colinear.begin();
		IndexedEvent* event = NULL;
		while (it != events_colinear.end() && (*it)->t_intersectant < t_1) {
			event = (*it);
			if (!rays[event->intersectant]->has_been_stopped() && !rays[event->intersected]->has_been_stopped()) {
				insert_indexed_event(event, intersectants, intersected);
				T.push_back(event);
			} else {
				delete event;
			}
			it++;
		}
		events_colinear.erase(events_colinear.begin(), it);
	}

	// Final step : obtains a list of chronologically sorted events and deletes the auxiliary vector
	//std::cout << T.size() << " events scheduled" << std::endl;
	std::sort(T.begin(), T.end(), before_indexed_event);

	if (T.size() > 0) {
		*schedule = T[0];
	} else {
		*schedule = NULL;
	}
	for (int i = 1; i < int(T.size()); i++) T[i]->previous = T[i - 1];
	for (int i = 0; i < int(T.size()) - 1; i++) T[i]->next = T[i + 1];
	T.clear();
}



void Propagation::insert_indexed_event(IndexedEvent* event, IndexedEvent** intersectants, IndexedEvent** intersected)
{
	int i = event->intersectant;
	int j = event->intersected;
	// Inserts this event at the appropriate location horizontally
	IndexedEvent* previous_event = NULL;
	IndexedEvent* current_event = intersectants[i];
	while (current_event != NULL && current_event->intersected < j) {
		previous_event = current_event;
		current_event = current_event->next_j;
	}
	if (current_event != NULL) {
		event->next_j = current_event;
		current_event->prev_j = event;
	}
	if (previous_event != NULL) {
		event->prev_j = previous_event;
		previous_event->next_j = event;
	} else {
		intersectants[i] = event;
	}
	// Inserts this event at the appropriate location vertically
	previous_event = NULL;
	current_event = intersected[j];
	while (current_event != NULL && current_event->intersectant < i) {
		previous_event = current_event;
		current_event = current_event->next_i;
	}
	if (current_event != NULL) {
		event->next_i = current_event;
		current_event->prev_i = event;
	}
	if (previous_event != NULL) {
		event->prev_i = previous_event;
		previous_event->next_i = event;
	} else {
		intersected[j] = event;
	}
}



void Propagation::insert_indexed_event(IndexedEvent* event, IndexedEvent** intersectants, IndexedEvent** intersected, IndexedEvent** intersectants_tails, IndexedEvent** intersected_tails) {
	int i = event->intersectant;
	int j = event->intersected;
	IndexedEvent* intersectant_tail = intersectants_tails[i];
	IndexedEvent* intersected_tail = intersected_tails[j];
	if (intersectant_tail == NULL) {
		intersectants[i] = event;
	} else {
		intersectant_tail->next_j = event;
		event->prev_j = intersectant_tail;
	}
	intersectants_tails[i] = event;
	if (intersected_tail == NULL) {
		intersected[j] = event;
	} else {
		intersected_tail->next_i = event;
		event->prev_i = intersected_tail;
	}
	intersected_tails[j] = event;
}



void Propagation::remove_indexed_event(IndexedEvent** schedule, IndexedEvent* indexed_event, IndexedEvent** intersectants, IndexedEvent** intersected, bool destroy) {
	IndexedEvent* previous = indexed_event->previous;
	IndexedEvent* next = indexed_event->next;
	int r = indexed_event->intersectant;
	int s = indexed_event->intersected;
	// Updates schedule
	if (previous != NULL) {
		previous->next = next;
	} else {
		*schedule = next;
	}
	if (next != NULL) {
		next->previous = previous;
	}
	// Updates sparse matrix
	IndexedEvent* previous_i = indexed_event->prev_i;
	IndexedEvent* previous_j = indexed_event->prev_j;
	IndexedEvent* next_i = indexed_event->next_i;
	IndexedEvent* next_j = indexed_event->next_j;
	if (previous_j == NULL) {
		intersectants[r] = next_j;
	} else {
		previous_j->next_j = next_j;
	}
	if (next_j != NULL) {
		next_j->prev_j = previous_j;
	}
	if (previous_i == NULL) {
		intersected[s] = next_i;
	} else {
		previous_i->next_i = next_i;
	}
	if (next_i != NULL) {
		next_i->prev_i = previous_i;
	}
	// Calls destructor if needed
	if (destroy) delete indexed_event;
}



void Propagation::remove_indexed_event(IndexedEvent** schedule, IndexedEvent* event_at_boundary, IndexedEvent** events_at_boundaries, double t_0, double t_1, bool destroy) {
	IndexedEvent* previous = event_at_boundary->previous;
	IndexedEvent* next = event_at_boundary->next;
	int r = event_at_boundary->intersectant;
	// Updates schedule if the event is included in it
	if (t_0 <= event_at_boundary->t_intersectant && event_at_boundary->t_intersectant < t_1) {
		if (previous != NULL) {
			previous->next = next;
		} else {
			*schedule = next;
		}
		if (next != NULL) {
			next->previous = previous;
		}
	}
	// Updates table
	events_at_boundaries[r] = NULL;
	// Calls destructor if needed
	if (destroy) delete event_at_boundary;
}



void Propagation::remove_references_to_intersectant_ray(IndexedEvent* upcoming_event, int index, IndexedEvent** schedule, IndexedEvent** intersectants, IndexedEvent** intersected, IndexedEvent** events_at_boundaries, double t_0, double t_1)
{
	// Removes all the events where ray 'ray_index' plays the role of intersectant
	IndexedEvent* current_event = intersectants[index];
	IndexedEvent* next_event = NULL;
	while (current_event != NULL) {
		next_event = current_event->next_j;
		remove_indexed_event(schedule, current_event, intersectants, intersected, current_event != upcoming_event);
		current_event = next_event;
	}

	// Removes some events where ray r plays the role of intersected
	current_event = intersected[index];
	next_event = NULL;
	while (current_event != NULL) {
		next_event = current_event->next_i;
		if (current_event->t_intersected > upcoming_event->t_intersectant) {
			remove_indexed_event(schedule, current_event, intersectants, intersected, current_event != upcoming_event);
		}
		current_event = next_event;
	}

	// Removes the event from the list of events involving boundaries
	IndexedEvent* event_at_boundary = events_at_boundaries[index];
	if (event_at_boundary != NULL) {
		remove_indexed_event(schedule, event_at_boundary, events_at_boundaries, t_0, t_1, event_at_boundary != upcoming_event);
	}
}



void Propagation::prune_events_and_build_edges(Kinetic_Model *model, Size2i & size, IndexedEvent** first_schedule, IndexedEvent** simplified_schedule, IndexedEvent** last_event, IndexedEvent** intersectants, IndexedEvent** intersected, IndexedEvent** events_at_boundaries,
    Quadtree* & quadtree, double t_0, double t_1, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
    vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges, vector<Vertex *> & initial_vertices_of_rays, int & active_rays)
{
	Vertex* intersection_point = NULL;
	bool is_corner = false;
	bool is_new_vertex = true;

    vector<Ray *> & rays = model->rays;
    Partition* & graph = model->graph;

	// Pops the first event of the raw schedule
	// If last_event == NULL, the processed event is the first to be processed
	// If last_event != NULL, some events have already be processed at a previous iteration, which means that the behaviour
	// of the algorithm depends on the value of curr_event : if it is not NULL, we append this event to the simplified schedule
	// with the help of last_event, otherwise we do nothing.

    IndexedEvent* curr_event = pop_event(model, size, first_schedule, intersectants, intersected, events_at_boundaries, quadtree, t_0, t_1, intersection_point, is_corner, is_new_vertex, active_rays);
	if (*last_event == NULL) {
		*simplified_schedule = curr_event;
		if (*simplified_schedule != NULL) {
			(*simplified_schedule)->previous = NULL;
            graph->build_edge(rays, curr_event, intersection_point, is_corner, is_new_vertex, outer_vertices, inner_vertices, outer_edges, inner_edges, initial_vertices_of_rays);
			*last_event = curr_event;
		}
	} else if (curr_event != NULL) {
		(*last_event)->next = curr_event;
		if (curr_event != NULL) {
			curr_event->previous = *last_event;
            graph->build_edge(rays, curr_event, intersection_point, is_corner, is_new_vertex, outer_vertices, inner_vertices, outer_edges, inner_edges, initial_vertices_of_rays);
			*last_event = curr_event;
		}
	}

	// While there are events to process, we append them to the simplified schedule using last_event.
	// Naturally, we update this variable at each iteration.
	while (curr_event != NULL) {
        curr_event = pop_event(model, size, first_schedule, intersectants, intersected, events_at_boundaries, quadtree, t_0, t_1, intersection_point, is_corner, is_new_vertex, active_rays);
		(*last_event)->next = curr_event;
		if (curr_event != NULL) {
			curr_event->previous = *last_event;
			graph->build_edge(rays, curr_event, intersection_point, is_corner, is_new_vertex, outer_vertices, inner_vertices, outer_edges, inner_edges, initial_vertices_of_rays);
			*last_event = curr_event;
		}
	}
}



IndexedEvent* Propagation::pop_event(Kinetic_Model* model, Size2i & size, IndexedEvent** schedule, IndexedEvent** intersectants, IndexedEvent** intersected, IndexedEvent** events_at_boundaries,
    Quadtree* & quadtree, double t_0, double t_1, Vertex* & intersection_point, bool & is_corner, bool & is_new_vertex, int & active_rays)
{
	if (*schedule == NULL) {
		return NULL;
	} else {
        vector<Ray *> & rays = model->rays;
        Parameters* params = model->params;
        Partition* & graph = model->graph;

		// Gets the first element of the schedule
		IndexedEvent* upcoming_event = *schedule;
		int i = upcoming_event->intersectant;
		int j = upcoming_event->intersected;
		double t_i = upcoming_event->t_intersectant, t_j = upcoming_event->t_intersected;
		
		// This upcoming event defines a new vertex in the graph, unless multiple rays intersect in the same point
		// So, let's find or create such a vertex
		Ray *r_i = rays[i], *r_j = (j >= 0 ? rays[j] : NULL);
#if NOT_MEASURING_PERFORMANCES
		bool intersects_brother_segment = (params->lsd_create_additional && j >= 0 ? r_i->parent->is_brother(r_j->parent) : false);
#else 
		bool intersects_brother_segment = false;
#endif
		Point2d pt;
		is_corner = Vertex::approximate_coordinates(r_i, r_j, t_i, model->I.rows, model->I.cols, pt);

		bool meets_colinear_ray = false;
		if (is_corner) {
			list<Vertex *> seeked_corner;
            quadtree->search(pt, 0.8 * params->prop_min_edge, seeked_corner);
			intersection_point = seeked_corner.front();
			is_new_vertex = false;
		} else {
			bool should_create_vertex = true;

			// If the intersectant ray is colinear to another ray,
			// we determine is the point where is meets the intersectant ray is issued from the intersection of a colinear ray
			list<Vertex *> neighbors;
			quadtree->search(pt, 1, neighbors);
			if (!neighbors.empty()) {
				if (r_j == NULL) {
					// If the intersection point falls on a boundary
					for (list<Vertex *>::iterator it_v = neighbors.begin() ; it_v != neighbors.end() ; ++it_v) {
						Vertex* v = (*it_v);
						if (v->outer_vertex_is_very_close(pt)) {
							intersection_point = v;
							intersection_point->events.push_back(upcoming_event);
							assert(intersection_point->has_no_duplicated_event());
							is_new_vertex = false;
							should_create_vertex = false;
							break;
						}
					}
				} else {
					for (list<Vertex *>::iterator it_v = neighbors.begin() ; it_v != neighbors.end() ; ++it_v) {
						Vertex* v = (*it_v);
						if (v->inner_vertex_created_by_colinear_ray(r_i, r_j, rays, false)) {
							intersection_point = v;
							intersection_point->events.push_back(upcoming_event);
							is_new_vertex = false;
							should_create_vertex = false;
							meets_colinear_ray = true;
							break;
						} else if (intersects_brother_segment && cv::norm(v->pt - pt) < 1e-9 && v->inner_vertex_barycenter_of_brother_segment(r_i, r_j, rays)) {
							intersection_point = v;
							intersection_point->events.push_back(upcoming_event);
							is_new_vertex = false;
							should_create_vertex = false;
							meets_colinear_ray = false;
							break;
						}
					}
				}
			}

			/*if (r_i->parent->node_colinear != NULL) {
				list<Vertex *> neighbor_vertices;
				quadtree->search(pt, 1, neighbor_vertices);
				if (!neighbor_vertices.empty()) {
					intersection_point = NULL;
					for (list<Vertex *>::iterator it_v = neighbor_vertices.begin(); it_v != neighbor_vertices.end(); it_v++) {
						if ((*it_v)->created_by_colinear_ray(r_i, r_i->parent->node_colinear, rays)) {
							intersection_point = (*it_v);
							is_new_vertex = false;
							intersection_point->events.push_back(upcoming_event);
							meets_colinear_ray = true;
							break;
						}
					}
				}
			}*/

			// If none of the conditions above could be fulfilled we create a new vertex
			if (should_create_vertex) {
				intersection_point = new Vertex(graph->get_id_vertices(), upcoming_event, pt);
				quadtree->add(intersection_point);
				is_new_vertex = true;
			}
		}
        // The ray r stops propagating if one of the four conditions is filled :
		// - it has interescted an image boundary
		// - it has intersected a colinear ray
		// - it has intersected a non-colinear ray in a point where it meets a colinear ray
		// We finally define the stop condition
		bool force_stop = !upcoming_event->is_colliding_ray || upcoming_event->is_colliding_colinear || meets_colinear_ray;
		if (force_stop) {
			r_i->stop();
		} else {
			if (r_i->primary_condition) {
				evaluate_primary_condition(r_i, r_j, t_i, params);
			}
			if (params->prop_extra_enabled) {
				if (!r_i->primary_condition || r_i->secondary_condition) {
#if NOT_MEASURING_PERFORMANCES
					evaluate_secondary_condition(r_i, model->I_grad_m, model->I_grad_t, pt, params);
#else
					r_i->secondary_condition = false;
#endif
					if (r_i->secondary_condition && r_i->t_swap > FLT_MAX / 2) r_i->t_swap = t_i;
				}
			}
			r_i->should_ray_be_stopped();
		}
		if (r_i->has_been_stopped()) {
			// If the ray stops its propagation, we erase :
			// - all future events actively involving rays[i]
			// - all events passively involving rays[i] which are posterior to the current time
			// - the event on the boundary, if it still exists
			remove_references_to_intersectant_ray(upcoming_event, i, schedule, intersectants, intersected, events_at_boundaries, t_0, t_1);

			// If, in addition to this, if the upcoming event is a collision between two colinear rays,
			// we set the time to leave of the other ray to 0 and remove all references to that ray
			if (upcoming_event->is_colliding_colinear) {
				r_j->stop();
				remove_references_to_intersectant_ray(upcoming_event, j, schedule, intersectants, intersected, events_at_boundaries, t_0, t_1);
				r_j->set_time(t_j);
				active_rays--;
			}
			r_i->set_time(t_i);
			active_rays--;
		} else {
			// Otherwise, if the ray doesn't stop its propagation, we simply erase the reference to this event in the matrix of events
			remove_indexed_event(schedule, upcoming_event, intersectants, intersected, false);
		}

		// Updates the pointers of the structure, and finally returns the event
		if (*schedule != NULL) (*schedule)->previous = NULL;
		upcoming_event->prev_i = NULL;
		upcoming_event->prev_j = NULL;
		upcoming_event->next_i = NULL;
		upcoming_event->next_j = NULL;
		return upcoming_event;
	}
}



void Propagation::evaluate_primary_condition(Ray* r_i, Ray* r_j, double t_i, Parameters* params)
{
#if NOT_MEASURING_PERFORMANCES
	if (t_i < 0 || r_i->parent->is_brother(r_j->parent)) {
#else
	if (t_i < 0) {
#endif
		r_i->primary_condition = true;
	} else {
		switch (params->prop_policy) {
		case 1:
			// The ray should be stopped if it intersects another ray for the k-th time
			r_i->ttl--;
			r_i->primary_condition = (r_i->ttl != 0);
			break;
		case 2:
			// The ray should be stopped if it intersects another ray for the k-th time
			// or if it intersects another non orthogonal ray
			r_i->ttl--;
			r_i->primary_condition = (r_i->ttl != 0 && r_i->intersects_orthogonal_ray(r_j));
			break;
		case 3:
			// The ray should be stopped if it has lived more than t units of time
			r_i->primary_condition = (t_i <= params->prop_distance);
			break;
		case 4:
			// The ray should be stopped if it has lived more than t units of time
			// or if it intersects another non orthogonal ray
			r_i->primary_condition = (t_i <= params->prop_distance && r_i->intersects_orthogonal_ray(r_j));
			break;
		}
	}
}


#if NOT_MEASURING_PERFORMANCES
void Propagation::evaluate_secondary_condition(Ray* r_i, Matrix<double> & I_m, Matrix<double> & I_t, Point2d & P, Parameters* params)
{
	if (!params->prop_check_m_enabled && !params->prop_check_t_enabled) return;

	int step = 1;
	uint rows = I_m.rows;
	uint cols = I_m.cols;
	typedef pair<int, int> Pixel;

	// Over a first phase we define a region of interest
	// It is defined by lines of pixels going from P to Q, Q being located at a certain distance from P.
	// This distance is defined by the user

	int r_width = int(params->prop_region_width);
	double r_length = params->prop_region_length;
	vector<vector<Pixel> > regions;
	define_region(r_i, P, r_width, r_length, I_m, regions);

	// Over a second phase we define our subregions, lines of pixes
	// And we measure quantities on them
	int sub_region_length = int(params->prop_sub_region_length);
	for (uint i = 0 ; i < regions.size() ; i++) {

		if (regions[i].size() < sub_region_length) continue;
		for (int j = 0; j <= regions[i].size() - sub_region_length ; j += step) {
			
			// Gets a subregion, and evaluates the condition
			vector<Pixel> sub_region = vector<Pixel>(regions[i].begin() + j, regions[i].begin() + j + sub_region_length);
			bool condition_satisfied = evaluate_secondary_condition_at_region_scale(r_i, sub_region, I_m, I_t, params);

			// Returns true, or iterates
			if (condition_satisfied) {

				r_i->secondary_condition = true;
				return;
			}
		}
	}

	// Reaching this point means that all the previous tests failed
	r_i->secondary_condition = false;
}
#endif


void Propagation::define_region(Ray* r_i, Point2d & P, int width, double length, Matrix<double> & I, vector<vector<pair<int, int> > > & region)
{
	typedef pair<int, int> Pixel;
	uint rows = I.rows, cols = I.cols;

	// A region of interest is defined by lines of pixels going from P to Q, Q being located at a distance r_length from P.
	// However if r_width > 1, P and Q move along the normal axis to the ray r_i.

	Vec2d r_dir = r_i->OA;
	Vec2d r_nor = Vec2d(-r_dir[1], r_dir[0]);

	Point2d Q = Point2d(P.x + length * r_dir[0], P.y + length * r_dir[1]);

	uint p_x = uint(jclamp(0, P.x, cols - 1)), p_y = uint(jclamp(0, rows - P.y, rows - 1));
	uint q_x = uint(jclamp(0, Q.x, cols - 1)), q_y = uint(jclamp(0, rows - Q.y, rows - 1));

	uint pu_x = uint(jclamp(0, P.x + (width / 2) * r_nor[0], cols - 1));
	uint pu_y = uint(jclamp(0, rows - P.y - (width / 2) * r_nor[1], rows - 1));
	uint pl_x = uint(jclamp(0, P.x - (width / 2) * r_nor[0], cols - 1));
	uint pl_y = uint(jclamp(0, rows - P.y + (width / 2) * r_nor[1], rows - 1));

	list<Pixel> normal_axis;
	I.line(pu_y, pu_x, pl_y, pl_x, normal_axis);

	region.clear();
	region.reserve(normal_axis.size());
	for (list<Pixel>::iterator it_p = normal_axis.begin() ; it_p != normal_axis.end() ; it_p++) 
	{
		// From the translated P, finds the translated Q
		uint pc_y = it_p->first, pc_x = it_p->second;
		int dx = int(pc_x) - int(p_x), dy = int(pc_y) - int(p_y);
		uint qc_x = uint(int(q_x) + dx), qc_y = uint(int(q_y) + dy);

		// Finds the line between the two points, and eliminates pixels that do not fall in the image
		list<Pixel> line;
		I.line(pc_y, pc_x, qc_y, qc_x, line);

		list<Pixel>::iterator it_p1 = line.begin();
		Pixel front = line.front();
		if (front.first < 0 || front.first >= int(rows) || front.second < 0 || front.second >= int(cols)) {

			// If the first pixel of the line is not in the image it is needless to include this line in the region
			// since it means that, obviously, we are too far away from the initial pixel, or there is not enough pixels
			// that we can measure
			continue;

		} else {
			
			// Removes the last pixels that don't fall within the image
			Pixel last = line.back();
			if (last.first < 0 || last.first >= int(rows) || last.second < 0 || last.second >= int(cols)) {
				line.reverse();
				list<Pixel>::iterator it_p2 = line.begin();
				while (it_p2->first < 0 || it_p2->first >= int(rows) || it_p2->second < 0 || it_p2->second >= int(cols)) ++it_p2;
				line.erase(line.begin(), it_p2);
				line.reverse();
			}

			// Adds this line
			region.push_back(vector<Pixel>(line.begin(), line.end()));
		}
	}
}


#if NOT_MEASURING_PERFORMANCES
bool Propagation::evaluate_secondary_condition_dot_product(Ray* r_i, vector<pair<int, int> > & region, Matrix<double> & M, Matrix<double> & T, Parameters* params)
{
	double dp = 0;
	Vec2d & r_dir = r_i->OA;
	for (uint r = 0 ; r < region.size() ; r++) {
		int i = region[r].first;
		int j = region[r].second;
		double m_ij = M(i, j);
		double t_ij = T(i, j);
		Vec2d g_dir = Vec2d(m_ij * cos(t_ij), m_ij * sin(t_ij));
		//dp += r_dir[0] * g_dir[0] + r_dir[1] * g_dir[1];
		dp += r_dir[0] * (-g_dir[1]) + r_dir[1] * g_dir[0];
	}
	dp /= region.size();
	return (fabs(dp) >= params->prop_dot_th);
}


bool Propagation::evaluate_secondary_condition_double_detector(Ray* r_i, vector<pair<int, int> > & region, Matrix<double> & M, Matrix<double> & T, Parameters* params)
{
	// Gets a reference value for the gradient's magnitude...
	double magnitude_ref;
	if (params->prop_m_compare_to == 0) {
		magnitude_ref = params->prop_m_factor * r_i->parent->grad_magnitude;
	} else if (params->prop_m_compare_to == 1) {
		magnitude_ref = params->prop_m_fixed_magnitude;
	}

	// ... and the gradient's orientation
	double theta_ref;
	if (params->prop_t_compare_to == 0) {
		theta_ref = r_i->parent->grad_theta;
	}

	// Vector revealing values taken by the gradient's magnitude and orientation in all points included in the region
	vector<double> v_magnitude(region.size(), 0);
	vector<double> v_theta(region.size(), 0);
	for (uint r = 0; r < region.size(); r++) {
		int i = region[r].first, j = region[r].second;
		v_magnitude[r] = M(i, j);
		v_theta[r] = T(i, j);
	}


	// Evaluates the condition for not disabling the ray

	if (params->prop_compared == 0) {
		// Case when we compare the mean gradient of the region
		bool check_m = true, check_t = true;

		// Tests if the mean magnitude is high enough
		if (params->prop_check_m_enabled) {
			double magnitude_mean = Means::mean(v_magnitude);
			check_m = (magnitude_mean >= magnitude_ref);
		}

		// Tests if the average angle is close enough to the reference angle
		if (params->prop_check_t_enabled) {
			check_t = false;
			double theta_mean;
			Means::circular_mean_180(v_theta, theta_mean);
			for (int k = -1; k <= 1; k++) {
				if (fabs(theta_mean - theta_ref + 180 * k) < params->prop_t_tolerance) {
					check_t = true;
					break;
				}
			}
		}

		// Returns the result
		return (check_m && check_t);

	} else if (params->prop_compared == 1) {
		// Case when we compare a ratio of pixels
		uint n = 0;

		for (uint r = 0; r < region.size(); r++) {
			// For each pixel of the region
			// We test if the local magnitude is high enough
			bool check_m_r = (params->prop_check_m_enabled ? (v_magnitude[r] >= magnitude_ref) : true);

			// We test if the local angle is close enough to the reference angle
			bool check_t_r = true;
			if (params->prop_check_t_enabled) {
				check_t_r = false;
				for (int k = -1; k <= 1; k++) {
					double theta_r = 180 * v_theta[r] / PI;
					if (theta_r < 0) theta_r += 180;
					if (fabs(theta_r - theta_ref + 180 * k) < params->prop_t_tolerance) {
						check_t_r = true; break;
					}
				}
			}

			// If the test (or the tests) are successful, the increment the number of successes
			if (check_m_r && check_t_r) ++n;
		}

		// The quality of the region is expressed as a ratio
		double ratio = double(n) / v_magnitude.size();
		if (ratio >= params->prop_ratio) {
			return true;
		} else {
			return false;
		}

	} else {
		return false;
	}
}


bool Propagation::evaluate_secondary_condition_at_region_scale(Ray* r_i, vector<pair<int, int> > & region, Matrix<double> & M, Matrix<double> & T, Parameters* params)
{
	if (params->prop_compared == 0 || params->prop_compared == 1) {
		return evaluate_secondary_condition_double_detector(r_i, region, M, T, params);
	} else if (params->prop_compared == 2) {
		return evaluate_secondary_condition_dot_product(r_i, region, M, T, params);
	}
	return false;
}
#endif


void Propagation::draw_rays(Kinetic_Model *model, Matrix<uchar> & J, double t_lim)
{
    vector<Ray *> & rays = model->rays;
    IndexedEvent* & schedule = model->schedule;
    Matrix<uchar> & I = model->I;

	// Counts the number of unfinished rays
	int n = 0;
	for (unsigned int r = 0; r < rays.size(); r++) {
		if (!rays[r]->has_been_stopped()) {
			rays[r]->set_time(t_lim);
			n++;
		}
	}
	if (n > 0) std::cout << "Warning : there are " << n << " unfinished rays" << std::endl;

	// We are going to loop on the list of events, starting by the last one
	IndexedEvent* current_event = schedule;
	IndexedEvent* last_event = NULL;
	while (current_event != NULL) {
		// Iterates, and leaves a pointer on the last event
		if (current_event->next == NULL) {
			last_event = current_event;
		}
		current_event = current_event->next;
	}

	// Reverses the schedule and finds the dates when other rays terminate
	current_event = last_event;
	while (current_event != NULL) {
		int r = current_event->intersectant;
		if (!rays[r]->is_time_set()) {
			rays[r]->set_time(current_event->t_intersectant);
		}
		if (current_event->is_colliding_colinear) {
			int s = current_event->intersected;
			rays[s]->set_time(current_event->t_intersected);
		}
		current_event = current_event->previous;
	}
	if (last_event != NULL) std::cout << "** Final t : " << last_event->t_intersectant << std::endl;

	uchar blue[3] = {0, 0, 255};
    J = Matrix<uchar>(I.rows, I.cols, 3);
	for (unsigned int i = 0 ; i < J.rows ; i++) {
		for (unsigned int j = 0 ; j < J.cols ; j++) {
			for (unsigned int c = 0 ; c < J.channels ; c++) {
                J(i, j, c) = I(i, j, c);
			}
		}
	}
	for (unsigned int i = 0 ; i < rays.size() ; i++) {
		Ray* r = rays[i];
		Point2d p1 = r->O;
		Point2d p2 = r->A + Point2d(r->t * r->OA);
        Point2i p1i = Point2i(int(round(jclamp(0, p1.x, J.cols - 1))), int(round(jclamp(0, J.rows - p1.y, J.rows - 1))));
        Point2i p2i = Point2i(int(round(jclamp(0, p2.x, J.cols - 1))), int(round(jclamp(0, J.rows - p2.y, J.rows - 1))));
		J.line(p1i.y, p1i.x, p2i.y, p2i.x, blue);
	}
}



void Propagation::write_schedule(Parameters* params, IndexedEvent* schedule)
{
	std::string filename = "schedule.txt";
	FILE* file = fopen(filename.c_str(), "w");
	if (file != NULL) {
		IndexedEvent* curr_event = schedule;
		while (curr_event != NULL) {
			fprintf(file, "%i %i %i %i %lf %lf\n", int(curr_event->is_colliding_ray), int(curr_event->is_colliding_colinear), curr_event->intersectant, curr_event->intersected, 
				curr_event->t_intersectant, curr_event->t_intersected);
			curr_event = curr_event->next;
		}
		fclose(file);
	}
}


#if NOT_MEASURING_PERFORMANCES
void Propagation::make_layer(Kinetic_Model *model)
{
	Matrix<uchar> & I = model->I;
	model->clear_line_items(model->L_prop);
	Edge* e = model->graph->edges_head;

	while (e != NULL) {
		Point2d pt_1 = e->v1->pt;
		Point2d pt_2 = e->v2->pt;
		Point2d pc_1 = Point2d(jclamp(0, pt_1.x, I.cols - 1), jclamp(0, I.rows - pt_1.y, I.rows - 1));
		Point2d pc_2 = Point2d(jclamp(0, pt_2.x, I.cols - 1), jclamp(0, I.rows - pt_2.y, I.rows - 1));
		model->add_line(model->L_prop, pc_1.x, pc_1.y, pc_2.x, pc_2.y, 255, 0, 0);
		e = e->e_next;
	}
}
#endif