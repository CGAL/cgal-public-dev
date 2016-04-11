#ifndef ARC_OCCURENCE_H
#define ARC_OCCURENCE_H

#include <CGAL/Compact_container.h>

namespace CGAL
{
    template<class Surface_>
    class Path; 
    
    template<class Surface_>
    class Arc_occurence  : public Compact_container_base
    {
        template<class>
        friend class Path;
        
    public:
        typedef Surface_ Surface;
        typedef typename Surface::Path_handle Path_handle;
        typedef typename Surface::Halfedge_handle Halfedge_handle;
        typedef typename Surface::Arc_occurence_handle Arc_occurence_handle;
        
        //constructor
        Arc_occurence(Path_handle path, Halfedge_handle halfedge):
                mPath(path),
                mHalfedge(halfedge)
        {
        }
        
        /// return the path.
        Path_handle path(){
            return mPath;
        }
        
        /// return the associated halfedge.
        Halfedge_handle halfedge(){
            return mHalfedge;
        }
        
    private:
        Path_handle mPath;
        Halfedge_handle mHalfedge;
        Arc_occurence_handle left;
        Arc_occurence_handle right;
    };
}

#endif