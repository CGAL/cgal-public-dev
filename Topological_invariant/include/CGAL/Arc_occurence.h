#ifndef ARC_OCCURENCE_H
#define ARC_OCCURENCE_H

#include <CGAL/Compact_container.h>

namespace CGAL
{
    template<class Surface_>
    class Path; 
    
    template<class Surface_, class Path_handle_, class halfedge_handle_>
    class Arc_occurence  : public Compact_container_base
    {
        template<class>
        friend class Path;
        
    public:
        typedef Surface_ Surface;
        typedef Path_handle_ Path_handle;
        typedef halfedge_handle_ halfedge_handle;
        
        typedef typename Surface::Arc_occurence_handle Arc_occurence_handle;       
        //constructor
        Arc_occurence(Path_handle path, halfedge_handle halfedge):
                mPath(path),
                mHalfedge(halfedge)
        {
        }
        
        /// return the path.
        Path_handle path(){
            return mPath;
        }
        
        /// return the associated halfedge.
        halfedge_handle halfedge(){
            return mHalfedge;
        }
        
    private:
        Path_handle mPath;
        halfedge_handle mHalfedge;
    };
}

#endif