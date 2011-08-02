// the functions from this file will be copied to the TDS and Tri_3

#include<cmath>

CGAL_BEGIN_NAMESPACE

#define TDS_FUNC(type) \
	template < class Vb, class Cb >\
	type\
	Triangulation_data_structure_3<Vb,Cb>::

#define TRI_FUNC(type) \
	template < class Vb, class Cb >\
	type\
	Triangulation_3<Vb,Cb>::

#define DT_FUNC(type) \
	template < class Vb, class Cb >\
	type\
	Delaunay_triangulation_3<Vb,Cb>::

/////////////////
// TDS structs //
/////////////////

template < class Vb, class Cb >
template < class Cell_handle > 
struct Triangulation_data_structure_3<Vb,Cb>::
less_Cell_handle
{    
    void get_id(const Cell_handle& cell, Vertex_handle& a0, Vertex_handle& a1, Vertex_handle& a2, Vertex_handle& a3) const
    {
        a0 = cell->vertex(0);
        a1 = cell->vertex(1);
        a2 = cell->vertex(2);
        a3 = cell->vertex(3);
    }
    
    bool operator()(const Cell_handle& a, const Cell_handle& b) const
    {
        Vertex_handle a0, a1, a2, a3;
        get_id(a, a0, a1, a2, a3);
        
        Vertex_handle b0, b1, b2, b3;
        get_id(b, b0, b1, b2, b3);
        
        if (a0 < b0) return true;
        if (a0 > b0) return false;
        
        if (a1 < b1) return true;
        if (a1 > b1) return false;
        
        if (a2 < b2) return true;
        if (a2 > b2) return false;
        
        if (a3 < b3) return true;
        return false;
    }
};

template < class Vb, class Cb >
template < class Vertex_handle > 
struct Triangulation_data_structure_3<Vb,Cb>::
less_Vertex_handle
{
	bool operator()(const Vertex_handle &s1, const Vertex_handle &s2) const
	{
		return s1 < s2;
	}
};

template < class Vb, class Cb >
template < class Edge >
struct Triangulation_data_structure_3<Vb,Cb>::
less_Edge
{
    void get_id(const Edge& edge, Vertex_handle& a, Vertex_handle& b) const
    {
        a = edge.first->vertex(edge.second);
        b = edge.first->vertex(edge.third );
        if ( !(a < b) ) std::swap(a, b);
    }
    
    bool operator()(const Edge& a, const Edge& b) const
    {
        std::pair<Vertex_handle, Vertex_handle> a_handle, b_handle;

        get_id(a, a_handle.first, a_handle.second);
        get_id(b, b_handle.first, b_handle.second);
        
        return a_handle < b_handle;
    }
};

/////////////////
// TDS members //
/////////////////
/*
template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
facet_vertex(Facet facet, int index)
{
	Cell_handle cell = facet.first;
	return cell->vertex(index);
}
*/
template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
get_source_vertex(const Edge& edge) const
{
	return edge.first->vertex(edge.second);
}

template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
get_target_vertex(const Edge& edge) const
{
	return edge.first->vertex(edge.third);
}

template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
get_any_other_vertex(
		Cell_handle cell,
		Vertex_handle va,
		Vertex_handle vb) const
{
	int i;
	for(i=0;i<4;i++)
		if(cell->vertex(i) != va
		&& cell->vertex(i) != vb)
			return cell->vertex(i);

	assert(false); // never come here
	return Vertex_handle();
}

// return the vertex of cell different from va, vb and vc
template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
get_remaining_vertex(Cell_handle cell,
	Vertex_handle va,
	Vertex_handle vb,
	Vertex_handle vc) const
{
	int i;
	for(i=0;i<4;i++)
		if(cell->vertex(i) != va
		&& cell->vertex(i) != vb
		&& cell->vertex(i) != vc)
			return cell->vertex(i);

	assert(false); // never come here
	return Vertex_handle();
}

template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
get_vertex_from_facet(const Facet& facet, int index) const
{
        return facet.first->vertex( (facet.second + (index%3) + 1) % 4 );
}

TDS_FUNC(void)
get_vertices_from_facet(const Facet& facet, 
                                Vertex_handle& va,
                                Vertex_handle& vb,
                                Vertex_handle& vc) const
{
	va = get_vertex_from_facet(facet, 0);
	vb = get_vertex_from_facet(facet, 1);
	vc = get_vertex_from_facet(facet, 2);
	
	if (facet.second % 2 == 0)
		std::swap(vb, vc);
}

TDS_FUNC(void)
get_vertices_from_edge(const Edge& edge, Vertex_handle& va, Vertex_handle& vb) const
{
	va = edge.first->vertex( edge.second );
	vb = edge.first->vertex( edge.third );
}

// GET LINK //
    
TDS_FUNC(void)
get_vertices_from_edge_link(const Edge& edge, Vertex_handle_set& vertices) const
{
	Vertex_handle s = get_source_vertex(edge);
	Vertex_handle t = get_target_vertex(edge);
	Cell_circulator ccirc = incident_cells(edge); 
	Cell_circulator cend  = ccirc;

	CGAL_For_all(ccirc, cend)
	{
		Cell_handle cell = ccirc;
		Vertex_handle u = get_any_other_vertex(cell, s, t);
		Vertex_handle v = get_remaining_vertex(cell, s, t, u);
		vertices.insert(u);
		vertices.insert(v);
	}
}
    
TDS_FUNC(void)
get_edges_from_edge_link(const Edge& edge, Edge_set& edges) const
{
	Vertex_handle s = get_source_vertex(edge);
	Vertex_handle t = get_target_vertex(edge); 
	Cell_circulator ccirc = incident_cells(edge); 
	Cell_circulator cend  = ccirc;
	CGAL_For_all(ccirc, cend)
	{
		Cell_handle cell = ccirc;
		Vertex_handle u = get_any_other_vertex(cell, s, t);
		Vertex_handle v = get_remaining_vertex(cell, s, t, u);
		Edge e(cell, cell->index(u), cell->index(v));
		edges.insert(e);
	}
}
    
TDS_FUNC(void)
get_vertices_from_vertex_link(Vertex_handle vertex, Vertex_handle_set& vertices) const
{
	std::list<Cell_handle> cells;		
	incident_cells(vertex, std::back_inserter(cells));
	
	typename std::list<Cell_handle>::const_iterator it;
	for (it = cells.begin(); it != cells.end(); it++) {
		Cell_handle cell = *it;
		Vertex_handle v0, v1, v2;
		Facet facet(cell, cell->index(vertex));
		get_vertices_from_facet(facet, v0, v1, v2);
		vertices.insert(v0);
		vertices.insert(v1);
		vertices.insert(v2);
	}
}
        
TDS_FUNC(void)
get_edges_from_vertex_link(Vertex_handle vertex, Edge_set& edges) const
{
	std::list<Cell_handle> cells;		
	incident_cells(vertex, std::back_inserter(cells));
	typename std::list<Cell_handle>::const_iterator it;
	for (it = cells.begin(); it != cells.end(); it++) {
		Cell_handle cell = *it;
		Vertex_handle v0, v1, v2;
		Facet facet(cell, cell->index(vertex));
		get_vertices_from_facet(facet, v0, v1, v2);
		int i0 = cell->index(v0);
		int i1 = cell->index(v1);
		int i2 = cell->index(v2);
		edges.insert( Edge(cell, i0, i1) );
		edges.insert( Edge(cell, i1, i2) );
		edges.insert( Edge(cell, i2, i0) );
	}
}

TDS_FUNC(void)
get_facets_from_link(Vertex_handle vertex, Vertex_handle dont_include, Facet_list& hull) const
{
	CGAL_triangulation_precondition( dimension() == 3 );

	std::list<Cell_handle> cells;
	incident_cells(vertex, std::back_inserter(cells));
	
	typename std::list<Cell_handle>::const_iterator it;
	for (it = cells.begin(); it != cells.end(); ++it) {
		Cell_handle cell = *it;
		if (!cell->has_vertex(dont_include)) {
			Facet facet(cell, cell->index(vertex));
			hull.push_back(facet);
		}
	}
}

TDS_FUNC(void)
get_edges_from_link(Vertex_handle vertex, Vertex_handle dont_include, Edge_list& hull) const
{
	CGAL_triangulation_precondition( dimension() == 2 );

	std::list<Cell_handle> cells;
	//incident_faces
	incident_cells(vertex, std::back_inserter(cells));
	
	typename std::list<Cell_handle>::const_iterator it;
	for (it = cells.begin(); it != cells.end(); ++it) {
		Cell_handle cell = *it;
		if (!cell->has_vertex(dont_include)) {
			int id = cell->index(vertex);
			Edge edge(cell, cw(id), ccw(id));
			hull.push_back(edge);
		}
	}
}

// O( |star(edge)|*log|star(edge)| ) – TODO: use hash?
TDS_FUNC(bool)
is_collapsible(const Edge& edge) const
{
	CGAL_triangulation_precondition( dimension() >= 1 );
	CGAL_triangulation_precondition( number_of_vertices() >= 2 );
	
	switch( dimension() ) {
		case 3: {
			bool for_vertices = is_collapsible_for_vertices(edge);
			bool for_edges = is_collapsible_for_edges(edge);

			std::cout << "top_for_vertices: " << for_vertices << std::endl;
			std::cout << "top_for_edges" << for_edges << std::endl;

			return for_vertices && for_edges;		
	
			return 	is_collapsible_for_vertices(edge)	
				&& is_collapsible_for_edges(edge);
		}
		case 2: {
			return is_collapsible_for_vertices(edge);
		}
		case 1: {
			CGAL_triangulation_expensive_precondition( is_edge(edge.first, edge.second, edge.third) );
			return true;
		}
	}

	CGAL_triangulation_assertion( false );
}
    
TDS_FUNC(bool)
is_collapsible_for_vertices(const Edge& edge) const
{
	CGAL_triangulation_expensive_precondition( is_cell(edge.first) );
	CGAL_triangulation_expensive_precondition( is_edge(edge.first, edge.second, edge.third) );

	Cell_handle cell = edge.first;
	Vertex_handle source = get_source_vertex(edge);
	Vertex_handle target = get_target_vertex(edge);

	Vertex_handle_set svertices;
	get_vertices_from_vertex_link(source, svertices);

	Vertex_handle_set tvertices;
	get_vertices_from_vertex_link(target, tvertices);

        Vertex_handle_set ivertices;
	typename Vertex_handle_set::const_iterator vit;
	for (vit = svertices.begin(); vit != svertices.end(); ++vit)
		if (tvertices.find(*vit) != tvertices.end())
			ivertices.insert(*vit);
                
	Vertex_handle_set evertices;
	get_vertices_from_edge_link(edge, evertices);

	if (ivertices.size() != evertices.size())
		return false;
 
	for (vit = ivertices.begin(); vit != ivertices.end(); ++vit)
		if (evertices.find(*vit) == evertices.end())
			return false;
        
	return true;
}
    
TDS_FUNC(bool)
is_collapsible_for_edges(const Edge& edge) const
{
	CGAL_triangulation_expensive_precondition( is_cell(edge.first) );
	CGAL_triangulation_expensive_precondition( is_edge(edge.first, edge.second, edge.third) );

	Cell_handle cell = edge.first;

	Vertex_handle source, target;
	source = get_source_vertex(edge);
	target = get_target_vertex(edge);

	Edge_set sedges, tedges;
	get_edges_from_vertex_link(source, sedges);
	get_edges_from_vertex_link(target, tedges);

	Edge_set iedges;
	typename Edge_set::const_iterator eit;
	for (eit = sedges.begin(); eit != sedges.end(); ++eit)
		if (tedges.find(*eit) != tedges.end()) {
			iedges.insert(*eit);

			CGAL_triangulation_assertion (eit->first->vertex( eit->second ) != edge.first->vertex( edge.second ));	
			CGAL_triangulation_assertion (eit->first->vertex( eit->second ) != edge.first->vertex( edge.third ));	
			CGAL_triangulation_assertion (eit->first->vertex( eit->third ) != edge.first->vertex( edge.second ));	
			CGAL_triangulation_assertion (eit->first->vertex( eit->third ) != edge.first->vertex( edge.third ));	
		}
        
	Edge_set eedges;
	get_edges_from_edge_link(edge, eedges);
 
	std::cout << iedges.size() << " <-> " << eedges.size() << std::endl; 
       
	if (iedges.size() != eedges.size())
		return false;

	for (eit = iedges.begin(); eit != iedges.end(); ++eit)
		if (eedges.find(*eit) == eedges.end())
			return false;

        return true;        
}

///////////////////
// COLLAPSE EDGE //
///////////////////

// O( |star(edge)| )
TDS_FUNC(bool)
collapse_edge(Edge& edge)
{
	CGAL_triangulation_precondition( is_collapsible(edge) );
	CGAL_triangulation_precondition( dimension() >= 1 );
	
	// two infinite and one finite
	CGAL_triangulation_precondition( number_of_vertices() >= 3 );

	// edge source-target: source is removed
	Cell_handle seed_cell = edge.first;
	Vertex_handle source = get_source_vertex(edge);
				//seed_cell->vertex(edge.second);
	Vertex_handle target = get_target_vertex(edge);
				//seed_cell->vertex(edge.third);
	
	CGAL_triangulation_precondition( source != target );
	CGAL_triangulation_precondition( is_vertex(source) );
	CGAL_triangulation_precondition( is_vertex(target) );
	CGAL_triangulation_precondition( is_simplex(seed_cell) );
	CGAL_triangulation_precondition( is_edge(edge.first, edge.second, edge.third) );

	switch ( dimension() ) {
		case 3: {
			// collect all cells incident to source - as this vertex will disappear
			// we will have to update source to target later in these cells.
			std::list<Cell_handle> scells;
			incident_cells(source, std::back_inserter(scells));

			// stores all cells incident to edge pq (they will be deleted)
			std::list<Cell_handle> rcells; 

			// revolve on cells around edge pq and update
			// incidence between mate cells.
			Cell_handle cell = seed_cell;
			Vertex_handle u = get_any_other_vertex(cell,target,source);

			// merge the two circle sets of cells
			do {
				// add cell to revolving cells
				rcells.push_back(cell);

				// get mate cells
				Vertex_handle v = get_remaining_vertex(cell,target,source,u);
				Cell_handle cell_suv = cell->neighbor(cell->index(target));
				Cell_handle cell_tuv = cell->neighbor(cell->index(source));
				Vertex_handle mt = mirror_vertex(cell,cell->index(target));
				Vertex_handle ms = mirror_vertex(cell,cell->index(source));

				// update incidence relationships
				int index_mt = cell_suv->index(mt);
				int index_ms = cell_tuv->index(ms);

				CGAL_triangulation_assertion(!cell_suv->has_neighbor(cell_tuv));
				CGAL_triangulation_assertion(!cell_tuv->has_neighbor(cell_suv));

				cell_suv->set_neighbor(index_mt, cell_tuv);
				cell_tuv->set_neighbor(index_ms, cell_suv);

				// set incident cell to ensure that vertices do not point
				// to a deleted cell (among scells).
				u->set_cell(cell_tuv);
				target->set_cell(cell_tuv);

				// revolves to neighboring cell
				cell = cell->neighbor(cell->index(u));
				u = v;
			} while(cell != seed_cell);

			// update cells incident to source
			typename std::list<Cell_handle>::iterator it;
			for(it = scells.begin(); it != scells.end(); it++) {
				Cell_handle cell = *it;
				cell->set_vertex(cell->index(source),target);
			}

			// delete cells and vertex source
			delete_cells(rcells.begin(),rcells.end());
			delete_vertex(source);

			// less than a tetrahedral
			if ( number_of_vertices() <= 4 ) 
				set_dimension(2);
			
			CGAL_triangulation_expensive_postcondition( is_valid() );
			return true;
		}
		case 2: {
			// The following drawing corrsponds to the variables
			// used in this part...
			// The vertex target is returned...
			//
			//        itl       i=v3      itr
			//         *---------*---------*
			//          \       / \       /
			//           \  tl /   \  tr /
			//            \   /  f  \   /
			//             \ /       \ /
			//target=ccw(i) *---------* cw(i)=source
			//             / \       / \
			//            /   \  g  /   \
			//           /  bl \   /  br \
			//          /       \ /       \
			//         *---------*---------*
			//        ibl       j=v4      ibr
			//                                                           
			// The situation after the "join"-operation is as follows:
			//
			//                 i
			//           *-----*-----*
			//            \    |    /
			//             \ tl|tr /
			//              \  |  /
			//               \ | /
			//                \|/
			//                 *  target
			//                /|\
			//               / | \
			//              /  |  \
			//             / bl|br \
			//            /    |    \
			//           *-----*-----*
			//

			int deg2 = degree(source);

			//CGAL_triangulation_precondition( deg2 >= 3 );

			//if ( deg2 == 3 ) {
			//	remove_degree_3(source, f->neighbor(ccw(i)));
			//	return target;
			//}
			  
			// first we register all the needed info
			Cell_handle f = seed_cell;
			int i = 3 - (edge.second + edge.third);
			Cell_handle g = f->neighbor(i);
			int j = mirror_index(f, i);

			Cell_handle tl = f->neighbor( cw(i)  );
			Cell_handle tr = f->neighbor( ccw(i) );

			int itl = mirror_index(f, cw(i)  );
			int itr = mirror_index(f, ccw(i) );

			Cell_handle bl = g->neighbor( ccw(j) );
			Cell_handle br = g->neighbor( cw(j)  );

			int ibl = mirror_index(g, ccw(j) );
			int ibr = mirror_index(g, cw(j)  );

			// we need to store the faces adjacent to source as well as the
			// indices of source w.r.t. these faces, so that afterwards we can set 
			// target to be the vertex for these faces
			std::list<Cell_handle> scells;
			//std::list<int> star_indices_of_source;

			incident_cells(source, std::back_inserter(scells));
			CGAL_triangulation_assertion( int(scells.size()) == deg2 );

			// from this point and on we modify the values
			
			// first set the neighbors
			tl->set_neighbor(itl, tr);
			tr->set_neighbor(itr, tl);
			bl->set_neighbor(ibl, br);
			br->set_neighbor(ibr, bl);
			
			// make sure that all the faces containing source
			// now contain target
			typename std::list<Cell_handle>::iterator cit;
			for (cit = scells.begin(); cit != scells.end(); cit++) {
				Cell_handle cell = *cit;
				int id = cell->index(source);
				CGAL_triangulation_assertion( cell->vertex(id) == source );
				cell->set_vertex( id, target );
			}

			// then make sure that all the vertices have correct pointers to 
			// faces
			Vertex_handle v3 = f->vertex(i);
			Vertex_handle v4 = g->vertex(j);
			//if ( v3->cell() == f )
				v3->set_cell(tr);
			//if ( v4->cell() == g )
				v4->set_cell(br);
			//if ( target->cell() == f || target->cell() == g )
				target->set_cell(tl);

			#ifndef CGAL_NO_ASSERTIONS
			  for (Cell_iterator cit = cells_begin(); cit != cells_end(); ++cit)
			    CGAL_triangulation_assertion( !cit->has_vertex(source) );
			#endif

			// memory management
			//star_faces_of_source.clear();
			//star_indices_of_source.clear();

			delete_cell(f);
			delete_cell(g);

			delete_vertex(source);

			// one finite plus the infinite one
			if (number_of_vertices() == 3)
				set_dimension(1);

			CGAL_triangulation_expensive_postcondition( is_valid() );
			return true;
		}
		case 1: {
			// before the collapse:
			//       prev_cell      seed_cell      next_cell
			// ---x--------------x------------>x----------------
			//                 source        target

			// after the collapse:
			//              prev_cell              next_cell
			// ---x----------------------------x----------------
			//                               target

			CGAL_triangulation_assertion( seed_cell->index(target) <= 1 && seed_cell->index(source) <= 1 );
			Cell_handle prev_cell = seed_cell->neighbor( seed_cell->index(target) );
			Cell_handle next_cell = seed_cell->neighbor( seed_cell->index(source) );
		
			//CGAL_triangulation_assertion( prev_cell->has_vertex(source) );
			//CGAL_triangulation_assertion( !prev_cell->has_vertex(target) );
	
			int source_index = prev_cell->index(source);
			prev_cell->set_vertex( source_index, target );
			prev_cell->set_neighbor( 1-source_index, next_cell );
			next_cell->set_neighbor( 1-next_cell->index(target), prev_cell );
			target->set_cell(next_cell);

			delete_cell(seed_cell);
			delete_vertex(source);

			// one finite plus the infinite one
			if (number_of_vertices() == 2)
				set_dimension(0);
			
			CGAL_triangulation_expensive_postcondition( is_valid() );
			return true;
		}
	}

	CGAL_triangulation_assertion( false );
}

////////////////////////////////////
// DELAUNAY_TRIANGULATION members //
////////////////////////////////////

DT_FUNC(bool)
is_collapsible(const Edge& edge, const Point& point) const
{
	return 	Tr_Base::is_collapsible(edge, point)
		&& check_delaunay_property(edge, point);
}

DT_FUNC(bool)
is_collapsible(const Edge& edge) const
{
	return 	is_collapsible(edge, Tr_Base::_tds.get_target_vertex(edge)->point());
}

// protected:
DT_FUNC(bool)
check_delaunay_property(const Edge& edge, const Point& point) const
{
	CGAL_triangulation_precondition( dimension() >= 1 );
	CGAL_triangulation_precondition( number_of_vertices() >= 2 );

	switch ( dimension() ) {
		case 3: {	
			Cell_handle seed_cell = edge.first;
			Vertex_handle source = seed_cell->vertex(edge.second);
			Vertex_handle target = seed_cell->vertex(edge.third);

			CGAL_triangulation_precondition( dimension() == 3 );
			CGAL_triangulation_precondition( !is_infinite(edge) );

			std::list<Cell_handle> scells;
			incident_cells(source, std::back_inserter(scells)); // precondition: dimension()==3

			typename std::list<Cell_handle>::iterator it;
			for(it = scells.begin(); it != scells.end(); it++) {
				Cell_handle c = *it;
				CGAL_assertion( c->has_vertex(source) );

				// shperes of infinite cells are not defined
				if (is_infinite(c)) continue;
				
				// the cells incident to both source and target are deleted due to the collapse
				if (c->has_vertex(target)) continue;

				// TODO: optimize without the array
				Point p[4] = {
					c->vertex(0)->point(),
					c->vertex(1)->point(),
					c->vertex(2)->point(),
					c->vertex(3)->point()
				};		

				// replace the source with the target vertex
				p[ c->index(source) ] = point;

				for (int i=0; i<4; i++ ) {
					Vertex_handle query = Tr_Base::mirror_vertex(c,i);
					//Vertex_handle query = c->neighbor(i)->vertex(c->neighbor(i)->index(c));
					if (is_infinite(query)) continue;

					// 'perturb'=false
					if ( (Bounded_side) side_of_oriented_sphere(p[0], p[1], p[2], p[3], query->point(), false) == ON_BOUNDED_SIDE )
						return false;
				}
			}
		}
		case 2: {
			//TODO:
			return true;
		}
		case 1: {
			//TODO;
			return true;	
		}
	}

	return true;
}
  
/////////////////////////////
// TRIANGULATION_3 members //
/////////////////////////////

TRI_FUNC(bool)
collapse_edge(Edge& edge, const Point& point)
{
	//TODO: check for moving target to an existing point
	CGAL_triangulation_precondition( is_geom_collapsible(edge, point) );

	Cell_handle cell = edge.first;
	Vertex_handle target = cell->vertex(edge.third);
	target->set_point(point);

	return collapse_edge(edge);
}

TRI_FUNC(bool)
collapse_edge(Edge& edge)
{
	CGAL_triangulation_precondition( is_geom_collapsible(edge) );

	bool res = _tds.collapse_edge(edge);

	switch( dimension() ) {
		case 3: {
			// TODO: remove flat cells
		}
		case 2: {
			// TODO: remove flat facets
		}
		case 1: {
			// TODO: nothing?
		}
	}
	

	return res;
}

TRI_FUNC(bool)
is_collapsible(const Edge& edge) const
{
	return	is_geom_collapsible(edge)
		&& _tds.is_collapsible(edge);
}

TRI_FUNC(bool)
is_collapsible(const Edge& edge, const Point& point) const
{
	return	is_geom_collapsible(edge, point)
		&& _tds.is_collapsible(edge);
}

TRI_FUNC(bool)
is_geom_collapsible(const Edge& edge) const
{
	return is_geom_collapsible(edge, _tds.get_target_vertex(edge)->point());
}

// O( |star(edge)| )
TRI_FUNC(bool)
is_geom_collapsible(const Edge& edge, const Point& point) const
{
	CGAL_triangulation_precondition( dimension() >= 1 );
	CGAL_triangulation_expensive_precondition( is_simplex(edge.first) );
	CGAL_triangulation_expensive_precondition( is_edge(edge.first, edge.second, edge.third) );

	Vertex_handle source = _tds.get_source_vertex(edge);
	Vertex_handle target = _tds.get_target_vertex(edge);

	if( is_infinite(source) )
		return false;
	
	switch( dimension() ) {
		case 3: {
			if (is_infinite(target)) return true;

			Facet_list hull;

			// then 'point' may not be visible from the link of 'source'
			if (point != source->point())
				_tds.get_facets_from_link(source, target, hull);
	
			// then 'point' may not be visible from the link of 'target'
			if (point != target->point()) 
				_tds.get_facets_from_link(target, source, hull);

			return is_in_kernel(point, hull.begin(), hull.end());
		}
		case 2: {
			if (is_infinite(target)) return true;

			Edge_list hull;

			// then 'point' may not be visible from the link of 'source'
			if (point != source->point())
				_tds.get_edges_from_link(source, target, hull);
	
			// then 'point' may not be visible from the link of 'target'
			if (point != target->point()) 
				_tds.get_edges_from_link(target, source, hull);

			return is_in_kernel_2d(point, hull.begin(), hull.end());
		}
		case 1: {
			return true;
		}
	}

	CGAL_triangulation_assertion( false );
}  

template < class Vb, class Cb >
template < class Iterator > // value_type = Face
bool
Triangulation_3<Vb,Cb>::
is_in_kernel(Point query, Iterator begin, Iterator end) const
{
	for (Iterator it = begin; it != end; ++it) {
		Facet facet = *it;
		if (is_infinite(facet))
			continue;
		
		Vertex_handle va, vb, vc;
		_tds.get_vertices_from_facet(facet, va, vb, vc);
	
		if (CGAL::orientation(query, va->point(), vb->point(), vc->point()) != CGAL::NEGATIVE)
			return false; 
	}

	return true;
}

template < class Vb, class Cb >
template < class Iterator > // value_type = Face
bool
Triangulation_3<Vb,Cb>::
is_in_kernel_2d(Point query, Iterator begin, Iterator end) const
{
	for (Iterator it = begin; it != end; ++it) {
		Edge edge = *it;
		if (is_infinite(edge))
			continue;

		Vertex_handle va, vb;
		_tds.get_vertices_from_edge(edge, va, vb);

		if (CGAL::coplanar_orientation(query, va->point(), vb->point()) != CGAL::NEGATIVE)
			return false;
	}

	return true;
}
// TODO: add into Tri_3?
TRI_FUNC(bool)
is_simplex( Cell_handle c ) const
{
  switch(dimension()) {
    case 3 : return _tds.is_cell(c);
    case 2 : return _tds.is_facet(c, 3);
    case 1 : return _tds.is_edge(c, 0, 1);
    case 0 : return _tds.is_vertex(c->vertex(0));
    case -1 : return c == _tds.cells().begin();
  }
  return false;
}

///////////
// DEBUG //
///////////

template < class Vb, class Cb >
template < class Cont >
void
Triangulation_data_structure_3<Vb,Cb>::
print_cells(const Cont S, std::string note) const
{
	std::cout << "Cells, " << note << std::endl;

	typename Cont::iterator cit;
	for (cit = S.begin(); cit != S.end(); cit++) {
		Cell_handle c = *cit;

		Vertex_handle v0 = c->vertex(0);
		Vertex_handle v1 = c->vertex(1);
		Vertex_handle v2 = c->vertex(2);
		Vertex_handle v3 = c->vertex(3);

		std::cout << "[  (" << v0->point() << ")   " << std::endl;
		std::cout << "   (" << v1->point() << ")   " << std::endl;
		std::cout << "   (" << v2->point() << ")   " << std::endl;
		std::cout << "   (" << v3->point() << ")  ]" << std::endl;
	}

	std::cout << std::endl;
}	

template < class Vb, class Cb >
template < class Cont >
void
Triangulation_data_structure_3<Vb,Cb>::
print_vertices(const Cont S, std::string note) const
{
	std::cout << "Vertices, " << note << std::endl;

	typename Cont::iterator vit;
	for (vit = S.begin(); vit != S.end(); vit++) {
		Vertex_handle v = *vit;
		std::cout << "(" << v->point() << ") ";
	}

	std::cout << std::endl;
}	

template < class Vb, class Cb >
template < class Cont >
void
Triangulation_data_structure_3<Vb,Cb>::
print_edges(const Cont S, std::string note) const
{
	std::cout << "Edges, " << note << std::endl;

	typename Cont::iterator eit;
	for (eit = S.begin(); eit != S.end(); eit++) {
		Edge e = *eit;
		std::cout << "(" << get_source_vertex(e)->point() << ") -> (" << get_target_vertex(e)->point() << ") " << std::endl;
	}

	std::cout << std::endl;
}	

CGAL_END_NAMESPACE

