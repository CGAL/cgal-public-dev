// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : SoX
// File          : include/SoX/Reorder_tree.h
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Lutz Kettner  <kettner@mpi-inf.mpg.de>
// ============================================================================


/*!  \file SoX/Reorder_tree.h  \brief definition of Reorder_tree class template
*/

#ifndef SoX_REORDER_TREE_H
#define SoX_REORDER_TREE_H 1

#include <SoX/basic.h>

#include <vector>
#include <boost/optional.hpp>
#include <boost/none.hpp>

namespace SoX {

template < class CurveId >
class Reorder_tree {
    
private:
    // forward declaration
    struct Link;
    
    // forward declaration
    struct Node;

   // private members
    
    typedef std::vector< Link > Link_container;
    
    typedef typename Link_container::iterator Link_iterator;

    // n-1 nodes + one sentinel, points to this->_m_links[i],
    // this->_m_tree[i] is on the level this->_m_mults[i] 
    // (so we might have unused spots here)
    std::vector< Node > _m_tree; 

    // linked lists of children
    // For each node and each _m_mults[i] i <= n, we have one parent -> 2*n
    std::vector< Link > _m_links;
    
    // first freenode entry in _m_links[]
    int _m_freenode;

    // Stack of nodes
    // Just nodes on the stack, points to _m_tree[n]
    std::vector< int > _m_stack;  

    // topmost element on the stack
    int _m_top;

    // avoid n^2 space on calling stack
    std::vector < int > _m_reverse_leafs; 
    
    int _m_reverse_n;

    int _m_n;
    
    
    // Tree representation
    struct Link { // linear linking of the list of children of a node
        boost::optional< Link_iterator > next;  // 0 for last link
        int   index; // child index, points to a tree node if < n, 
                     // otherwise _m_curves[index - n] is a leaf
        Link() {}
        Link(int i) : index(i) {}
        Link(Link_iterator n, int i) : next(n), index(i) {}
    };
    
    struct Node {
        // start of children, if existing
        boost::optional< Link_iterator > begin; 
        // end of children list (constant time append)
        boost::optional< Link_iterator > end; 
        Node() : begin(boost::none), end(boost::none) {}
        Node(Link_iterator b) : begin(b), end(b) {}
        Node(Link_iterator b, Link_iterator e) : begin(b), end(e) {}

        Node& operator=(const Node &n) {
            if (this != &n) {
                if (n.begin) {
                    SoX_assert(n.end);
                    this->begin = *(n.begin);
                    this->end = *(n.end);
                }
            }
            return *this;
        }
    };
    
    void _add_child(int node, int i) { 
        //  i < n, child is node, otherwise Y[i-n]
        SoX_assert(this->_m_freenode < 2*this->_m_n);
        this->_m_links[this->_m_freenode] = Link(i);
        Link_iterator lit = this->_m_links.begin();
        std::advance(lit, this->_m_freenode);
        if (!this->_m_tree[node].begin) {
            this->_m_tree[node] = Node(lit);
        } else {
            (*(this->_m_tree[node].end))->next = lit;
            this->_m_tree[node] = 
                Node(*this->_m_tree[node].begin, lit);
        }
        ++this->_m_freenode;
    }
    
    
    void _build_tree() {
        const int n = this->_m_n;
        
        // init sentinel and stack
        // build sentinel node this->_m_tree[0]. 
        // Its level this->_m_mults[0] is 0.
        this->_m_tree.clear();
        this->_m_tree.resize(n);
        this->_m_tree[0] = Node();  
        
        // this->_m_freenode list starts at first element
        this->_m_freenode = 0;
        
        // sentinel on stack 
        this->_m_stack.clear();
        this->_m_stack.resize(n);
        this->_m_stack[0] = 0;
        
        this->_m_links.clear();
        this->_m_links.resize(2*n);

        this->_m_reverse_leafs.clear();
        this->_m_reverse_leafs.resize(2*n);

        // one element on stack
        this->_m_top = 0;
        
     
        // loop scanning this->_m_mults[i];
        for (int i = 0; i < n; ++i) {
            std::cout << "Adding leaf " << i 
                      << " and node " << i+1 << std::endl;
            bool y_delayed = true; // does Y[i] belong to current or next node?
            if (this->_m_mults[i] > this->_m_mults[i+1] ) { // current node
                y_delayed = false;
                this->_add_child(this->_m_stack[this->_m_top], i+n);
            } // else wait until we pushed the new node on the stack
            
            // Find first node on stack whose 
            // this->_m_mults[this->_m_stack[k]] 
            // <= this->_m_mults[i+1]. (Note the sentinel)

            // On the way, add each node as child to the next node on 
            // the way to this->_m_stack[k] (so with index one less). 
            // The sequence of nodes on the stack has strictly 
            // monotone decreasing F values. Remove the nodes from 
            // the stack on the way towards this->_m_stack[k].
            // The last node this->_m_stack[k+1] gets linked to node 
            // i+1 if this->_m_mults[this->_m_stack[k]] < this->_m_mults[i+1],
            
            // otherwise it get linked to this->_m_stack[k].
            int link_k = -1;
            
            std::cout << "multsize: " << _m_mults.size() << std::endl;
            std::cout << "stacksize: " << _m_stack.size() << std::endl;
            std::cout << "top1: " << this->_m_top << std::endl;
            std::cout << "stop: " << this->_m_stack[this->_m_top] << std::endl;
            std::cout << "mstop: " << this->_m_mults[this->_m_stack[this->_m_top]] << std::endl;
            std::cout << "i: " << i << std::endl;
            std::cout << "mi1: " << this->_m_mults[i+1] << std::endl;
                                                     
            while (this->_m_mults[this->_m_stack[this->_m_top]] > 
                   this->_m_mults[i+1]) {
                SoX_assert(this->_m_top > 0);
                --this->_m_top;
                
                SoX_assert(this->_m_mults[this->_m_stack[this->_m_top]] < 
                           this->_m_mults[this->_m_stack[this->_m_top+1]]);
                std::cout << "multsize: " << _m_mults.size() << std::endl;
                std::cout << "stacksize: " << _m_stack.size() << std::endl;
                std::cout << "top2: " << this->_m_top << std::endl;
                std::cout << "mtop: " << this->_m_stack[this->_m_top] << std::endl;
                if (this->_m_mults[this->_m_stack[this->_m_top]] >= 
                    this->_m_mults[i+1]) {
                    this->_add_child(this->_m_stack[this->_m_top], 
                                     this->_m_stack[this->_m_top+1]);
                } else {
                    link_k = this->_m_stack[this->_m_top+1];
                }
            }
            // Create new node on stack 
            // if this->_m_mults[this->_m_stack[this->_m_top]] > 
            // this->_m_mults[i+1],
            // otherwise ignore this duplicate.
            if (this->_m_mults[this->_m_stack[this->_m_top]] < 
                this->_m_mults[i+1]) {
                this->_m_tree[i+1] = Node(); // no children yet
                ++this->_m_top;
                this->_m_stack[this->_m_top] = i+1;
                // Check if there is a node waiting to get 
                // added from above loop
                if (link_k != -1) {
                    this->_add_child(i+1, link_k);
                }
            }
            if (y_delayed) {
                // now add Y[i] to the (new or previously existing) 
                // node on stack
                this->_add_child(this->_m_stack[this->_m_top], i+n);
            }
        }
        // this->_m_tree[0] should be only node on stack
        SoX_assert(this->_m_top == 0 && this->_m_stack[0] == 0); 
    }
    
#if 0    
    void print_tree( int node = 0, int level = 0) {
        const char spaces[81] = "|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   ";
        const int n = this->_m_n;
        if ( node < n) {
            std::cout << (spaces + 80 - 4 * level) 
                      << "+-- Level " << _m_curves[node] 
                      << " (node " << node << ")" << std::endl;
            while (this->_m_tree[node].begin) {
                Link_iterator l = this->_m_tree[node].begin
                print_tree( l->index, level+1);
                l = l->next;
            }
        } else {
            std::cout << (spaces + 80 - 4 * level) << "+-- Curve " 
                //<< _m_curves[node-n] 
                      << std::endl;
        }
    }

    void print_reversed_tree( int node = 0, int level = 0) {
    if ( node == 0)
        reverse_n = 0;
    if ( node < n) {
        cout << (spaces + 80 - 4 * level) << "+-- Level " << F[node] 
             << " (node " << node << ")" << endl;
        if (( F[node] & 1) == 0) { // even case as before
            Link* l = T[node].begin;
            while ( l != 0) {
                print_reversed_tree( l->index, level+1);
                l = l->next;
            }
        } else { // odd case needs reversal
            // copy leafs in reverser array
            int i = reverse_n;
            int j = i;
            Link* l = T[node].begin;
            while ( l != 0) {
                reverse_leafs[j++] = l->index;
                l = l->next;
            }
            reverse_n = j;
            for ( int k = j-1; k >= i; --k) {
                print_reversed_tree( reverse_leafs[k], level+1);
            }
        }
    } else {
        cout << (spaces + 80 - 4 * level) << "+-- Curve " << Y[node-n] <<endl;
    }
    }
#endif



    void _extract_reversed_curves(int node = 0) {
        const int n = this->_m_n;
        // could be done non-recursive later if of interest
        if (node == 0) {
            this->_m_reverse_n = 0;
            this->_m_curves_out.clear();
            this->_m_curves_out.reserve(n);
        }
        std::cout << "erc1" << std::endl;
        if (node < n) {
            std::cout << "erc2" << std::endl;
            if ((this->_m_mults[node] & 1) == 0) { // even case
                boost::optional < Link_iterator > l =
                    this->_m_tree[node].begin;
                while (l) {
                    std::cout << "erc3" << std::endl;
                    _extract_reversed_curves((*l)->index);
                    l = (*l)->next;
                }
            } else { // odd case needs reversal
                // copy leafs in reverser array
                std::cout << "erc4" << std::endl;
                int i = this->_m_reverse_n;
                int j = i;
                boost::optional < Link_iterator > l =
                    this->_m_tree[node].begin;
                while (l) {
                    std::cout << "erc5" << std::endl;
                    this->_m_reverse_leafs[j++] = (*l)->index;
                    l = (*l)->next;
                }
                this->_m_reverse_n = j;
                for (int k = j-1; k >= i; --k) {
                    std::cout << "erc6" << std::endl;
                    _extract_reversed_curves(this->_m_reverse_leafs[k]);
                }
            }
        } else {
            std::cout << "erc7" << std::endl;
            this->_m_curves_out.push_back(this->_m_curves[node-n]);
        }
    }
    
public:
    //! this instance template parameter
    typedef CurveId Curve_id;

    //! constructs the tree
    template < class CurveIdIterator, class MultIterator > 
    Reorder_tree(CurveIdIterator cbegin, CurveIdIterator cend,
                 MultIterator mbegin, MultIterator mend) {
        this->_m_n = static_cast< int >(std::distance(cbegin,cend));
        int n2 = static_cast< int >(std::distance(mbegin,mend));
        SoX_assert(this->_m_n == n2 + 1);
        
        this->_m_curves.clear();
        this->_m_curves.reserve(this->_m_n);
        this->_m_curves_out.reserve(this->_m_n);
        std::copy(cbegin, cend, std::back_inserter(this->_m_curves));
        
        std::cout << "n1: " << _m_n << std::endl;
        std::cout << "n2: " << n2 << std::endl;
        this->_m_mults.clear();
        this->_m_mults.reserve(n2 + 2);
        this->_m_mults.push_back(0);
        std::copy(mbegin, mend, std::back_inserter(this->_m_mults));
        this->_m_mults.push_back(0);
        SoX_assert(static_cast< int >(this->_m_mults.size()) == n2 + 2);

        for (int i = 0; i < n2+2; i++) {
            std::cout << "mults" << i << "=" << _m_mults[i] << std::endl;
        }

        // build tree
        this->_build_tree();
        
        //this->print_tree();


        // traverse tree and store output to _m_curves
        this->_extract_reversed_curves();
    }
    
    
    // container types
    typedef typename std::vector< Curve_id >::value_type value_type;
    typedef typename std::vector< Curve_id >::iterator iterator;
    typedef typename std::vector< Curve_id >::const_iterator const_iterator;
    typedef typename std::vector< Curve_id >::reverse_iterator 
    reverse_iterator;
    typedef typename std::vector< Curve_id >::const_reverse_iterator 
    const_reverse_iterator;
    typedef typename std::vector< Curve_id >::reference reference;
    typedef typename std::vector< Curve_id >::const_reference const_reference;
    typedef typename std::vector< Curve_id >::pointer pointer;
    typedef typename std::vector< Curve_id >::difference_type difference_type;
    typedef typename std::vector< Curve_id >::size_type size_type;
    
    iterator begin() {
        return this->_m_curves_out.begin();
    }

    iterator end() {
        return this->_m_curves_out.end();
        
    }

    const_iterator begin() const {
        return this->_m_curves_out.begin();
        
    }

    const_iterator end() const {
        return this->_m_curves_out.end();
        
    }

    reverse_iterator rbegin() {
        return this->_m_curves_out.rbegin();
    }

    reverse_iterator rend() {
        return this->_m_curves_out.rend();
        
    }

    const_reverse_iterator rbegin() const {
        return this->_m_curves_out.rbegin();
        
    }

    const_reverse_iterator rend() const {
        return this->_m_curves_out.rend();
        
    }
    
    size_type size() {
        return this->_m_curves_out.size();
    }

    size_type max_size() {
        return this->_m_curves_out.max_size();
    }
    
    bool empty() {
        return this->_m_curves_out.empty();
    }
    
private:
    //! stored curves
    std::vector < Curve_id > _m_curves;

    //! stored mults
    std::vector< int > _m_mults;
    
    //! output curves
    std::vector < Curve_id > _m_curves_out;
};




} // namespace SoX

#endif // SoX_REORDER_TREE_H
// EOF
