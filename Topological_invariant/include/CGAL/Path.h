#ifndef PATH_H
#define PATH_H

#include <CGAL/Arc_occurence.h>
#include <list>

namespace CGAL
{
    template<typename Items_, typename Alloc_>
    class Topological_surface;
    
    template<class Surface_>
    class Path : public Compact_container_base
    {
        template<class, class>
        friend class Topological_surface;
    public:
        // types
        typedef Surface_ Surface;
        typedef Path<Surface> Self;
        
        typedef typename Surface::Halfedge_handle Halfedge_handle;
        typedef typename Surface::Path_handle Path_handle;
        
        typedef Arc_occurence<Surface, Path_handle, Halfedge_handle> Arc_occurence_;
        typedef CGAL_ALLOCATOR(Arc_occurence_)    Arc_occurence_allocator;
        typedef std::list<Arc_occurence_,Arc_occurence_allocator>      Arc_occurence_container;
        typedef typename Arc_occurence_container::iterator               Arc_occurence_handle;
        typedef typename Arc_occurence_container::const_iterator         Arc_occurence_const_handle;
        
        typedef typename Arc_occurence_container::size_type size_type;
        
        //constructor
        Path()
        {
        }
        
        Path(const Path& copy):
                mHandle(copy.mHandle),
                mArcs(copy.mArcs)
        {
            
        }
        
        //access
        size_type number_of_arc_occurence() const{
            return mArcs.size();
        }
        bool is_loop() const;
        
        Halfedge_handle halfedge(Arc_occurence_handle ao)const{
            return ao->halfedge();
        }
        
        // Range 
        //Arc_occurence_range arc_occurences();
        //Arc_occurence_revert_range arc_occurences_revert();
        
        // Modifier
        Arc_occurence_handle push_front(const Path& p);
        
        Arc_occurence_handle push_front(Halfedge_handle he){
            mArcs.push_front(Arc_occurence_(mHandle, he));
            Arc_occurence_handle arc = mArcs.begin();
            return arc;
        }
        
        Arc_occurence_handle push_back(const Path& p);
        
        Arc_occurence_handle push_back(Halfedge_handle d){
            mArcs.push_back(Arc_occurence_(mHandle, d));
            Arc_occurence_handle arc = mArcs.end();
            --arc;
            return arc;
        }
        
        void pop_front(){
            mArcs.pop_front();
        }
        
        void pop_back(){
            mArcs.pop_back();
        }
        
        // opperation
        Arc_occurence_handle detour_face(Arc_occurence_handle ao); 
        Arc_occurence_handle detour_face_opposite(Arc_occurence_handle ao);    
        
      /*  // Order around an edge 
        Arc_occurence_handle left(Arc_occurence_handle a){
            Order_iterator o = a->mOrderIterator;
            --o;
            return *o;
        }
        
        Arc_occurence_handle right(Arc_occurence_handle a){
            Order_iterator o = a->mOrderIterator;
            ++o;
            return *o;
        }
        
        Arc_occurence_handle move_left(Arc_occurence_handle a);
        Arc_occurence_handle move_right(Arc_occurence_handle);*/
        
    protected:
        Path_handle mHandle;
        Arc_occurence_container mArcs;
    };
}

#endif