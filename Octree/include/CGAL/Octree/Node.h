#ifndef OCTREE_NODE_H
#define OCTREE_NODE_H

#include <array>
#include <memory>
#include <bitset>

namespace CGAL {
  namespace Octree {
    namespace Node {

      /*!
       * \ingroup PkgOctreeClasses
       *
       * \brief A class representing a single node of the tree. Alternatively referred to as a cell, octant, or subtree
       *
       * \details The role of the node isn't fully stable yet
       *
       * \tparam Value is the datatype the node will contain
       */
      template<typename Value>
      class Node {

      public:

        /// \name Types
        /// @{

        /*!
         * \brief Array for containing the children of this node
         */
        typedef std::array<Node<Value>, 8> Child_list;

        /*!
         * \brief Set of bits representing this node's relationship to its parent
         *
         * \todo This deserves a more in-depth description
         */
        typedef std::bitset<3> Index;

        /*!
         * \brief Coordinate location representing this node's relationship with its parent
         *
         * \todo This deserves a more in-depth description
         */
        typedef std::array<uint32_t, 3> Int_location;

        /// @}

      private:

        Value m_value;

        const Node<Value> *m_parent;

        uint8_t m_depth;

        Int_location m_location;

        std::unique_ptr <Child_list> m_children;

      public:

        /// \name Construction
        /// @{

        /*!
         * \brief Creates a new node, optionally as the child of a parent
         *
         * \todo This warrant further explanation
         *
         * \param parent
         * \param index
         */
        Node(Node<Value> *parent = nullptr, Index index = 0) : m_parent(parent), m_depth(0), m_location({0, 0, 0}) {

          if (parent) {

            m_depth = parent->m_depth + 1;

            for (int i = 0; i < 3; i++)
              m_location[i] = (2 * m_location[i]) + index[i];

          }
        }

        /// @}

        /// \name Mutators
        /// @{

        /*!
         * \brief Split a node into subnodes
         *
         * \todo
         */
        void split() {

          m_children = std::make_unique<Child_list>();
          for (int index = 0; index < 8; index++) {

            (*m_children)[index] = std::move(Node<Value>(this, {index}));
          }
        }

        /*!
         * \brief Eliminate this node's children, making it a leaf node
         *
         * \todo
         */
        void unsplit() {

          m_children.reset();
        }

        /// @}

        /// \name Child Accessors
        /// @{

        /*!
         * \brief Access the child nodes of this node by their indices
         *
         * \todo
         *
         * \param index
         * \return
         */
        Node<Value> &operator[](int index) {
          return (*m_children)[index];
        }

        /*!
         * \brief Read-only access the child nodes of this node by their indices
         *
         * \todo
         *
         * \param index
         * \return
         */
        const Node<Value> &operator[](int index) const {
          return (*m_children)[index];
        }

        /// @}

        /// \name Property Accessors
        /// @{

        /*!
         * \brief Read-only access to this node's parent
         *
         * \return
         */
        const Node<Value> *parent() const { return m_parent; }

        /*!
         * \brief Retrieve this node's depth in the tree
         * \return
         */
        const uint8_t &depth() const { return m_depth; }

        /*!
         * \brief Retrieve this node's location in the tree
         * \return
         */
        Int_location location() const {

          return m_location;
        }

        /*!
         * \brief Retrieve this node's index in relation to its parent
         * \return
         */
        Index index() const {

          // TODO: There's a better way of doing this

          Index result;

          result[0] = location()[0] & 1;
          result[1] = location()[1] & 1;
          result[2] = location()[2] & 1;

          return result;
        }

        /*!
         * \brief Determine whether this node is a leaf node
         * \return
         */
        bool is_leaf() const { return (!m_children); }

        /*!
         * \brief Determine whether this node is the root node
         * \return
         */
        bool is_root() const { return (!m_parent); }

        /// @}

        /// \name Value Accessors
        /// @{

        /*!
         * \brief Access to the content held by this node
         * \return
         */
        Value &value() { return m_value; }

        /*!
         * \brief Read-only access to the content held by this node
         * \return
         */
        const Value &value() const { return m_value; }

        /// @}

        /// \name Operators
        /// @{

        /*!
         * \brief Compare this node to another node
         *
         * \param rhs
         * \return
         */
        bool operator==(Node &rhs) const {

          // Compare the values they contain
          if (m_value != rhs.m_value)
            return false;

          // If one node is a leaf, and the other isn't, they're not the same
          if (is_leaf() != rhs.is_leaf())
            return false;

          // If both nodes are non-leaf nodes
          if (!is_leaf()) {

            // Check all the children
            for (int i = 0; i < 8; ++i) {

              // If any child cell is different, they're not the same
              if (!((*m_children)[i] == rhs[i]))
                return false;
            }
          }

          // If both nodes are leaf nodes they must be identical (no check necessary)

          // If all other checks pass, the two trees are identical
          return true;
        }

        /// @}
      };
    }

  }
}

#endif //OCTREE_NODE_H
