// Copyright (c) 2014  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
// All rights reserved.
//

// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Fernando de Goes, Pierre Alliez

#ifndef RECONSTRUCTION_EDGE_2_H_
#define RECONSTRUCTION_EDGE_2_H_

template <class FT, class Edge, class Vertex_handle, class Face_handle>
class Reconstruction_edge_2
{
protected:
    Edge m_edge;
    Vertex_handle m_source;
    Vertex_handle m_target;

    FT  m_before_cost;
    FT  m_after_cost;

public:
    Reconstruction_edge_2()
    {
        m_edge = Edge(Face_handle(), 0);
        m_source = Vertex_handle();
        m_target = Vertex_handle();

        m_before_cost = 0.0;
        m_after_cost = 0.0;
    }

    Reconstruction_edge_2(const Reconstruction_edge_2& pedge)
    {
        m_edge = pedge.edge();
        m_source = pedge.source();
        m_target = pedge.target();

        m_before_cost = pedge.before();
        m_after_cost = pedge.after();
    }

    Reconstruction_edge_2(const Edge& edge,
           const FT before,
           const FT after)
    {
        m_edge = edge;
        get_vertices();

        m_before_cost = before;
        m_after_cost = after;
    }

    Reconstruction_edge_2(const Edge& edge,
           const FT priority = 0.0)
    {
        m_edge = edge;
        get_vertices();

        m_before_cost = 0.0;
        m_after_cost = priority;
    }

    Reconstruction_edge_2(Vertex_handle source, Vertex_handle target)
    {
        m_edge = Edge(Face_handle(), 0);
        m_source = source;
        m_target = target;

        m_before_cost = 0.0;
        m_after_cost = 0.0;
    }

    virtual ~Reconstruction_edge_2() { }

    Reconstruction_edge_2& operator = (const Reconstruction_edge_2& pedge)
    {
        m_edge = pedge.edge();
        m_source = pedge.source();
        m_target = pedge.target();

        m_before_cost = pedge.before();
        m_after_cost = pedge.after();

        return *this;
    }

    bool operator == (const Reconstruction_edge_2& pedge) const
    {
        return (m_source->id() == pedge.source()->id() &&
                m_target->id() == pedge.target()->id());
    }

    bool operator < (const Reconstruction_edge_2& pedge) const
    {
        if (m_source->id() < pedge.source()->id()) return true;
        if (m_source->id() > pedge.source()->id()) return false;

        if (m_target->id() < pedge.target()->id()) return true;
        return false;
    }

    const Edge& edge() const { return m_edge; }

    const Vertex_handle& source() const { return m_source; }

    const Vertex_handle& target() const { return m_target; }

    const FT before() const { return m_before_cost; }

    const FT after() const { return m_after_cost; }

    const FT priority() const { return after() - before(); }

protected:
    void get_vertices()
    {
        int index = m_edge.second;
        m_source = m_edge.first->vertex( (index+1)%3 );
        m_target = m_edge.first->vertex( (index+2)%3 );
    }
};

#endif


/* RECONSTRUCTION_EDGE_2_H_ */
