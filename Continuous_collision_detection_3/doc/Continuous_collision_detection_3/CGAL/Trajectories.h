
namespace CGAL {

/*!
\ingroup PkgCollisions3Classes

The class `Point_3_trajectory` serves as a container for the current and next positions of a `Point_3` object
during the time of interest for a collision query.
*/
template <class K>
class Point_3_trajectory : public std::tuple<typename K::Point_3, typename K::Point_3> {
    public:
        /// \name Types
        /// @{

        /*!
        the type of `Point_3` being tracked
        */
        typedef typename K::Point_3               Point;

        /*!
        the type from which `Point_3_trajectory` derives
        */
        typedef          std::tuple<Point, Point> Base;
        
        /// @}

        /// \name Creation
        /// @{

        /*!
        creates an empty trajectory
        */
        Point_3_trajectory() {}

        /*!
        creates a trajectory provided a two `Point_3` objects, corresponding to the current and next position of a point
        */
        Point_3_trajectory(const Point& current, const Point& next)
        : Base(current, next)
        {}


        /*!
        creates a trajectory provided a pair of `Point_3` objects, corresponding to the current and next position of a point
        */
        Point_3_trajectory(const Base& point_pair)
        : Base(point_pair)
        {}

        /// @}


        /// \name Methods
        /// @{

        /*!
        returns a `Point_3` object corresponding to the position of the point at the begining of the time period under consideration
        */
        Point current() const;

        /*!
        returns a `Point_3` object corresponding to the position of the point at the end of the time period under consideration
        */
        Point next() const;
        /// @}
};


/*!
\ingroup PkgCollisions3Classes

The class `Segment_3_trajectory` serves as a container for the current and next positions of a `Segment_3` object
during the time of interest for a collision query.
*/
template <class K>
class Segment_3_trajectory : public std::tuple<Point_3_trajectory<K>, Point_3_trajectory<K>> {

    public:
        /// \name Types
        /// @{

        /*!
        the type of `Segment_3` being tracked
        */
        using Segment            = typename K::Segment_3;

        /*!
        the type of `Point_3_trajectory` used to characterize the segment's trajectory
        */
        using Point_trajectory   = Point_3_trajectory<K>;
        
        /*!
        the type from which `Segment_3_trajectory` derives
        */
        using Base               = std::tuple<Point_trajectory, Point_trajectory>;

        /// @}


        /// \name Creation
        /// @{

        /*!
        creates an empty trajectory
        */
        Segment_3_trajectory();

        /*!
        creates a segment trajectory provided two `Point_3_trajectory` objects, corresponding to the trajectories of the source and target of the `Segment_3` object under consideration
        */
        Segment_3_trajectory(const Point_trajectory& source, const Point_trajectory& target);

        /*!
        creates a segment trajectory provided two `Point_3_trajectory` objects, corresponding to the trajectories of the source and target of the `Segment_3` object under consideration
        */
        Segment_3_trajectory(Point_trajectory&& source, Point_trajectory&& target);

        /*!
        creates a segment trajectory provided a tuple of two `Point_3_trajectory` objects, corresponding to the trajectories of the source and target of the `Segment_3` object under consideration
        */
        Segment_3_trajectory(const Base& point_trajectory_pair);

        /// @}


        /// \name Methods
        /// @{

        /*!
        returns a `Segment_3` object corresponding to the position of the segment at the begining of the time period under consideration
        */
        Segment current() const;

        /*!
        returns a `Segment_3` object corresponding to the position of the segment at the end of the time period under consideration
        */
        Segment next() const;

        /// @}
};

/*!
\ingroup PkgCollisions3Classes

The class `Triangle_3_trajectory` serves as a container for the current and next positions of a `Triangle_3` object
during the time of interest for a collision query.
*/
template <class K>
class Triangle_3_trajectory : public std::tuple<Point_3_trajectory<K>, Point_3_trajectory<K>, Point_3_trajectory<K>> {

    public:

        /// \name Types
        /// @{

        /*!
        the type of `Triangle_3` being tracked
        */
        using Triangle           = typename K::Triangle_3;

        /*!
        the type of `Point_3_trajectory` used to characterize the triangle's trajectory
        */
        using Point_trajectory   = Point_3_trajectory<K>;
        
        /*!
        the type from which `Segment_3_trajectory` derives
        */
        using Base               = std::tuple<Point_trajectory, Point_trajectory, Point_trajectory>;

        /// @}


        /// \name Creation
        /// @{

        /*!
        creates an empty trajectory
        */
        Triangle_3_trajectory();

        /*!
        creates a segment trajectory provided three `Point_3_trajectory` objects, corresponding to the trajectories of the three vertices of the `Triangle_3` object under consideration
        */
        Triangle_3_trajectory(const Point_trajectory& v0, const Point_trajectory& v1, const Point_trajectory& v2);

        /*!
        creates a segment trajectory provided three `Point_3_trajectory` objects, corresponding to the trajectories of the three vertices of the `Triangle_3` object under consideration
        */
        Triangle_3_trajectory(Point_trajectory&& v0, Point_trajectory&& v1, Point_trajectory&& v2);

        /*!
        creates a segment trajectory provided a tuple of three `Point_3_trajectory` objects, corresponding to the trajectories of the three vertices of the `Triangle_3` object under consideration
        */
        Triangle_3_trajectory(const Base& point_trajectory_triplet);

        /// @}


        /// \name Methods
        /// @{

        /*!
        returns a `Triangle_3` object corresponding to the position of the triangle at the begining of the time period under consideration
        */
        Triangle current() const;

        /*!
        returns a `Triangle_3` object corresponding to the position of the triangle at the end of the time period under consideration
        */
        Triangle next() const;

        /// @}
};



}