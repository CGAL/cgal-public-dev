#ifndef TOPOLOGICAL_SURFACE_H
#define TOPOLOGICAL_SURFACE_H

#include <CGAL/Generalized_map.h>
#include <CGAL/Generalized_map_constructors.h>
#include <CGAL/Path.h>
#include <queue>

namespace CGAL
{
    template<
            typename Items_ = CGAL::Generalized_map_min_items<2> , 
            typename Alloc_ = CGAL_ALLOCATOR(int)
            >
    class Topological_surface
    {
    public:
         //Constructor
        Topological_surface(){
            mSignature = get_new_mark();
        }

        template<class GMap>
        Topological_surface(const GMap& gmap): mGMap(gmap) {
            mSignature = get_new_mark();
        }
        
        //Basic types
        typedef Alloc_ Alloc;
        typedef Items_ Items;
        
        typedef Generalized_map<2, Items, Alloc> GMap;
        
        typedef Topological_surface<Items, Alloc> Self;
        
        //GMap types
        typedef typename GMap::Dart Dart;
        typedef typename GMap::Dart_handle Dart_handle;
        typedef typename GMap::Dart_const_handle Dart_const_handle;
        typedef typename GMap::size_type size_type;
        
        //halfedge
        struct Halfedge_handle
        {
            Halfedge_handle(Self& s, Dart_handle d):
                    dart(d),
                    surface(&s)
            {}
            
            Halfedge_handle()
            {}
            
            Halfedge_handle(const Halfedge_handle& h):
                dart(h.dart),
                surface(NULL)
            {
            }
            
            inline bool operator==(const Halfedge_handle& h){ 
                if(surface==NULL)return false;
                if(surface->is_free<2>(dart)){
                    return dart==h.dart;
                }else{
                    return dart==h.dart || h.dart==surface->alpha<2>(dart);
                }
            }
            inline bool operator!=(const Halfedge_handle& h){
                return !(dart==h.dart);
            }
            
            Dart_handle dart;
            Self* surface;
        };
        typedef const Halfedge_handle Halfedge_const_handle;
        
        //Path
        typedef Path<Self>                                          Path_;
        typedef typename Alloc::template rebind<Path_>::other       Path_allocator;
        typedef Compact_container<Path_,Path_allocator>             Path_container;
        typedef typename Path_container::iterator                   Path_handle;
        typedef typename Path_container::const_iterator             Path_const_handle;
        
        typedef typename Path_::Arc_occurence_handle Arc_occurence_handle;
        
        //Constants
        static const unsigned int dimension = GMap::dimension;
        static const size_type NB_MARKS = GMap::NB_MARKS;

        //Types for Attributes
        typedef typename GMap::Attributes Attributes;

        template<int i>
        struct Attribute_type: public GMap::template Attribute_type<i>
        {};

        template<int i>
        struct Attribute_handle: public GMap::template Attribute_handle<i>
        {};

        template<int i>
        struct Attribute_const_handle: public GMap::template Attribute_const_handle<i>
        {};

        //Range Types
        typedef typename GMap::Dart_range Dart_range;
        typedef typename GMap::Dart_const_range Dart_const_range;

        template <unsigned int i>
        struct Attribute_range: public GMap::template Attribute_range<i>
        {};

        template <unsigned int i>
        struct Attribute_const_range: public GMap::template Attribute_const_range<i>
        {};
        
        #ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
            template<unsigned int... Alpha>
            struct Dart_of_orbit_range: public GMap::template Dart_of_orbit_range<Alpha...>
            {};

            template<unsigned int ... Alpha>
            struct Dart_of_orbit_const_range: public GMap::template Dart_of_orbit_const_range<Alpha...>
            {};
        #else
            template<int B1=-1,int B2=-1,int B3=-1,int B4=-1,int B5=-1,
                int B6=-1,int B7=-1,int B8=-1,int B9=-1>
            struct Dart_of_orbit_range: public GMap::template Dart_of_orbit_range<B1, B2, B3, B4, B5, B6, B7, B8, B9>
            {
                typedef typename GMap::template Dart_of_orbit_range<B1, B2, B3, B4, B5, B6, B7, B8, B9> Base;

                Dart_of_orbit_range(GMap &amap, Dart_handle adart) : Base(amap,adart)
                {}
            };
            
            template<int B1=-1,int B2=-1,int B3=-1,int B4=-1,int B5=-1,
                int B6=-1,int B7=-1,int B8=-1,int B9=-1>
            struct Dart_of_orbit_const_range: public GMap::template Dart_of_orbit_const_range<B1, B2, B3, B4, B5, B6, B7, B8, B9>
            {
                typedef typename GMap::template Dart_of_orbit_const_range<B1, B2, B3, B4, B5, B6, B7, B8, B9> Base;

                Dart_of_orbit_const_range(GMap &amap, Dart_handle adart) : Base(amap,adart)
                {}
            };
        #endif

        template<unsigned int i,unsigned int dim=dimension>
        struct Dart_of_cell_range: public GMap::template Dart_of_cell_range<i, dim>
        {};

        template<unsigned int i,unsigned int dim=dimension>
        struct Dart_of_cell_const_range: public GMap::template Dart_of_cell_const_range<i, dim>
        {};

        template<unsigned int i,unsigned int j,unsigned int dim=dimension>
        struct One_dart_per_incident_cell_range: public GMap::template One_dart_per_incident_cell_range<i,j,dim>
        {};

        template<unsigned int i,unsigned int j,unsigned int dim=dimension>
        struct One_dart_per_incident_cell_const_range: public GMap::template One_dart_per_incident_cell_const_range<i, j, dim>
        {};

        template<unsigned int i,unsigned int dim=dimension>
        struct One_dart_per_cell_range: public GMap::template One_dart_per_cell_range<i, dim>
        {};

        template<unsigned int i,unsigned int dim=dimension>
        struct One_dart_per_cell_const_range: public GMap::template One_dart_per_cell_const_range<i, dim>{
            typedef typename GMap::template One_dart_per_cell_const_range<i, dim> Base;
            
            One_dart_per_cell_const_range(const Self &amap):
                    Base(amap.mGMap)
            {
                
            }
        };
        
        //Access Member Functions
        bool is_empty() const{
            return mGMap.is_empty();
        }

        bool is_valid() const{
            typedef Dart_const_range Range;
            typedef typename Range::const_iterator const_iterator;
            const_iterator it = darts().begin();
            for(; it!=darts().end(); ++it){
                if(!is_free<2>(it)){
                    Dart_const_handle d1 = it;
                    Dart_const_handle d2 = alpha<2>(d1);
                    
                    if(dart_signature(d1)==dart_signature(d2)){
                        std::cout<<"a"<<std::endl;
                        return false;
                    }
                }
                if(!is_free<1>(it)){
                    Dart_const_handle d1 = it;
                    Dart_const_handle d2 = alpha<1>(d1);
                    
                    if(dart_signature(d1)==dart_signature(d2)){
                        std::cout<<"a"<<std::endl;
                        return false;
                    }
                }
            }
            return mGMap.is_valid();
        }

        bool is_without_boundary(unsigned int i) const{
            return mGMap.is_without_boundary(i);
        }

        bool is_without_boundary() const{
            return mGMap.is_without_boundary();
        }

        size_type number_of_darts() const{
            return mGMap.number_of_darts();
        }

        template <unsigned int i>
        size_type number_of_attributes() const{
            return mGMap.number_of_attributes<i>();
        }

        bool is_dart_used(Dart_const_handle dh) const{
            return mGMap.is_dart_used(dh);
        }
        
        #ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
            template<typename ...Alphas>
            Dart_handle alpha(Dart_handle ADart, Alphas... alphas){
                return mGMap.alpha<Alphas>(ADart, alphas);
            }
                
            template<typename ...Alphas>
            Dart_const_handle alpha(Dart_const_handle ADart, Alphas... alphas) const{
                return mGMap.alpha<Alphas>(ADart, alphas);
            }
                
            template<int... Alphas>
            Dart_handle alpha(Dart_handle ADart){
                return mGMap.alpha<Alphas>(ADart);
            }
                
            template<int... Alphas>
            Dart_const_handle alpha(Dart_const_handle ADart) const{ 
                return mGMap.alpha<Alphas>(ADart);
            }
        #else
            Dart_handle alpha(Dart_handle ADart, int B1)
            { return mGMap.alpha(ADart, B1); }
            Dart_handle alpha(Dart_handle ADart, int B1, int B2)
            { return mGMap.alpha(ADart, B1, B2); }
            Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3)
            { return mGMap.alpha(ADart, B1, B2, B3); }
            Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                                int B4)
            { return mGMap.alpha(ADart, B1, B2, B3, B4); }
            Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                                int B4, int B5)
            { return mGMap.alpha(ADart, B1, B2, B3, B4, B5); }
            Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                                int B4, int B5, int B6)
            { return mGMap.alpha(ADart, B1, B2, B3, B4, B5, B6); }
            Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                                int B4, int B5, int B6, int B7)
            { return mGMap.alpha(ADart, B1, B2, B3, B4, B5, B6, B7); }
            Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                                int B4, int B5, int B6, int B7, int B8)
            { return mGMap.alpha(ADart, B1, B2, B3, B4, B5, B6, B7, B8); }
            Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                                int B4, int B5, int B6, int B7, int B8, int B9)
            { return mGMap.alpha(ADart, B1, B2, B3, B4, B5, B6, B7, B8, B9); }

            template<int B1>
            Dart_handle alpha(Dart_handle ADart)
            { return mGMap.alpha<B1>(ADart); }
            template<int B1, int B2>
            Dart_handle alpha(Dart_handle ADart)
            { return mGMap.alpha<B1, B2>(ADart); }
            template<int B1, int B2, int B3>
            Dart_handle alpha(Dart_handle ADart)
            { return mGMap.alpha<B1, B2, B3>(ADart); }
            template<int B1, int B2, int B3, int B4>
            Dart_handle alpha(Dart_handle ADart)
            { return mGMap.alpha<B1, B2, B3, B4>(ADart); }
            template<int B1, int B2, int B3, int B4, int B5>
            Dart_handle alpha(Dart_handle ADart)
            { return mGMap.alpha<B1, B2, B3, B4, B5>(ADart); }
            template<int B1, int B2, int B3, int B4, int B5, int B6>
            Dart_handle alpha(Dart_handle ADart)
            { return mGMap.alpha<B1, B2, B3, B4, B5, B6>(ADart); }
            template<int B1, int B2, int B3, int B4, int B5, int B6,
                    int B7>
            Dart_handle alpha(Dart_handle ADart)
            { return mGMap.alpha<B1, B2, B3, B4, B5, B6, B7>(ADart); }
            template<int B1, int B2, int B3, int B4, int B5, int B6,
                    int B7, int B8>
            Dart_handle alpha(Dart_handle ADart)
            { return mGMap.alpha<B1, B2, B3, B4, B5, B6, B7, B8>(ADart); }
            template<int B1, int B2, int B3, int B4, int B5, int B6,
                    int B7, int B8, int B9>
            Dart_handle alpha(Dart_handle ADart)
            { return mGMap.alpha<B1, B2, B3, B4, B5, B6, B7, B8, B9>(ADart); }

            Dart_const_handle alpha(Dart_const_handle ADart, int B1) const
            { return mGMap.alpha(ADart, B1); }
            Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2) const
            { return mGMap.alpha(ADart, B1, B2); }
            Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3) const
            { return mGMap.alpha(ADart, B1, B2, B3); }
            Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                                int B4) const
            { return mGMap.alpha(ADart, B1, B2, B3, B4); }
            Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                                int B4, int B5) const
            { return mGMap.alpha(ADart, B1, B2, B3, B4, B5); }
            Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                                int B4, int B5, int B6) const
            { return mGMap.alpha(ADart, B1, B2, B3, B4, B5, B6); }
            Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                                int B4, int B5, int B6, int B7) const
            { return mGMap.alpha(ADart, B1, B2, B3, B4, B5, B6, B7); }
            Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                                int B4, int B5, int B6, int B7, int B8) const
            { return mGMap.alpha(ADart, B1, B2, B3, B4, B5, B6, B7, B8); }
            Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                                int B4, int B5, int B6, int B7, int B8, int B9) const
            { return mGMap.alpha(ADart, B1, B2, B3, B4, B5, B6, B7, B8, B9); }
            
            template<int B1>
            Dart_const_handle alpha(Dart_const_handle ADart) const
            { return mGMap.alpha<B1>(ADart); }
            template<int B1, int B2>
            Dart_const_handle alpha(Dart_const_handle ADart) const
            { return mGMap.alpha<B1, B2>(ADart); }
            template<int B1, int B2, int B3>
            Dart_const_handle alpha(Dart_const_handle ADart) const
            { return mGMap.alpha<B1, B2, B3>(ADart); }
            template<int B1, int B2, int B3, int B4>
            Dart_const_handle alpha(Dart_const_handle ADart) const
            { return mGMap.alpha<B1, B2, B3, B4>(ADart); }
            template<int B1, int B2, int B3, int B4, int B5>
            Dart_const_handle alpha(Dart_const_handle ADart) const
            { return mGMap.alpha<B1, B2, B3, B4, B5>(ADart); }
            template<int B1, int B2, int B3, int B4, int B5, int B6>
            Dart_const_handle alpha(Dart_const_handle ADart) const
            { return mGMap.alpha<B1, B2, B3, B4, B5, B6>(ADart); }
            template<int B1, int B2, int B3, int B4, int B5, int B6,
                    int B7>
            Dart_const_handle alpha(Dart_const_handle ADart) const
            { return mGMap.alpha<B1, B2, B3, B4, B5, B6, B7>(ADart); }
            template<int B1, int B2, int B3, int B4, int B5, int B6,
                    int B7, int B8>
            Dart_const_handle alpha(Dart_const_handle ADart) const
            { return mGMap.alpha<B1, B2, B3, B4, B5, B6, B7, B8>(ADart); }
            template<int B1, int B2, int B3, int B4, int B5, int B6,
                    int B7, int B8, int B9>
            Dart_const_handle alpha(Dart_const_handle ADart) const
            { return mGMap.alpha<B1, B2, B3, B4, B5, B6, B7, B8, B9>(ADart); }
        #endif

        bool is_free(Dart_const_handle dh, unsigned int i) const{
            return mGMap.is_free(dh, i);
        }

        template<unsigned int i>
        bool is_free(Dart_const_handle dh) const{
            return mGMap.is_free<i>(dh);
        }

        int highest_nonfree_dimension(Dart_const_handle dh) const{
            return mGMap.highest_nonfree_dimension(dh);
        }

        Dart_handle opposite(Dart_handle dh){
            return mGMap.opposite(dh);
        }

        Dart_const_handle opposite(Dart_const_handle dh) const{
            return mGMap.opposite(dh);
        }

        Dart_handle other_extremity(Dart_handle dh){
            return mGMap.other_extremity(dh);
        }

        Dart_const_handle other_extremity(Dart_const_handle dh) const{
            return mGMap.other_extremity(dh);
        }

        std::ostream& display_characteristics(std::ostream & os) const{
            return mGMap.display_characteristics(os);
        }

        //Attributes Access Member Functions
        template <unsigned int i>
        typename Attribute_handle<i>::type attribute(Dart_handle dh){
            return mGMap.attribute<i>(dh);
        }

        template <unsigned int i>
        typename Attribute_const_handle<i>::type attribute(Dart_const_handle dh) const{
            return mGMap.attribute<i>(dh);
        }

        template<unsigned int i>
        Dart_handle dart_of_attribute(typename Attribute_handle<i>::type ah){
            return mGMap.dart_of_attribute<i>(ah);
        }

        template<unsigned int i>
        Dart_const_handle dart_of_attribute(typename Attribute_const_handle<i>::type ah) const{
            return mGMap.dart_of_attribute<i>(ah);
        }

        template <unsigned int i>
        typename Attribute_type<i>::type::Info& info_of_attribute(typename Attribute_handle<i>::type ah){
            return mGMap.info_of_attribute<i>(ah);
        }

        template <unsigned int i>
        const typename Attribute_type<i>::type::Info& info_of_attribute(typename Attribute_const_handle<i>::type ah) const{
            return mGMap.info_of_attribute<i>(ah);
        }

        template<unsigned int i>
        typename Attribute_type<i>::type::Info & info(Dart_handle adart){
            return mGMap.info<i>(adart);
        }

        template<unsigned int i>
        const typename Attribute_type<i>::type::Info & info(Dart_const_handle adart) const{
            return mGMap.info<i>(adart);
        }

        template<unsigned int i>
        Dart_handle & dart(Dart_handle adart){
            return mGMap.dart<i>(adart);
        }

        template<unsigned int i>
        Dart_const_handle dart(Dart_const_handle adart) const{
            return mGMap.dart<i>(adart);
        }

        template<unsigned int i>
        bool is_attribute_used(typename Attribute_const_handle<i>::type ah) const{
            return mGMap.is_attribute_used<i>(ah);
        }

        // Transformations Between Handles and Instances
        Dart_handle dart_handle(Dart& d){
            return mGMap.dart_handle(d);
        }

        Dart_const_handle dart_handle(const Dart& d) const{
            return mGMap.dart_handle(d);
        }

        // Range Access Member Functions
        Dart_range& darts(){
            return mGMap.darts();
        }

        Dart_const_range& darts() const{
            return mGMap.darts();
        }

        template<unsigned int i>
        typename Attribute_range<i>::type & attributes(){
            return mGMap.attributes<i>();
        }

        template<unsigned int i>
        typename Attribute_const_range<i>::type & attributes() const{
            return mGMap.attributes<i>();
        }
        
        template <unsigned int B1>
        Dart_of_orbit_range<B1> darts_of_orbit(Dart_handle adart){
            return Dart_of_orbit_range<B1>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2>
        Dart_of_orbit_range<B1,B2> darts_of_orbit(Dart_handle adart){
            return Dart_of_orbit_range<B1,B2>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3>
        Dart_of_orbit_range<B1,B2,B3> darts_of_orbit(Dart_handle adart){
            return Dart_of_orbit_range<B1,B2,B3>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4>
        Dart_of_orbit_range<B1,B2,B3,B4> darts_of_orbit(Dart_handle adart){
            return Dart_of_orbit_range<B1,B2,B3,B4>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
                unsigned int B5>
        Dart_of_orbit_range<B1,B2,B3,B4,B5> darts_of_orbit(Dart_handle adart){
            return Dart_of_orbit_range<B1,B2,B3,B4,B5>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
                unsigned int B5,unsigned int B6>
        Dart_of_orbit_range<B1,B2,B3,B4,B5,B6> darts_of_orbit(Dart_handle adart){
            return Dart_of_orbit_range<B1,B2,B3,B4,B5,B6>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
                unsigned int B5,unsigned int B6,unsigned int B7>
        Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7> darts_of_orbit(Dart_handle adart){
            return Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
                unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8>
        Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7,B8> darts_of_orbit(Dart_handle adart){
            return Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7,B8>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
                unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8,
                unsigned int B9>
        Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7,B8,B9> darts_of_orbit(Dart_handle adart){
            return Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>(mGMap, adart);
        }
        
        template <unsigned int B1>
        Dart_of_orbit_const_range<B1> darts_of_orbit(Dart_const_handle adart) const {
            return Dart_of_orbit_const_range<B1>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2>
        Dart_of_orbit_const_range<B1,B2> darts_of_orbit(Dart_const_handle adart) const {
            return Dart_of_orbit_const_range<B1,B2>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3>
        Dart_of_orbit_const_range<B1,B2,B3> darts_of_orbit(Dart_const_handle adart) const {
            return Dart_of_orbit_const_range<B1,B2,B3>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4>
        Dart_of_orbit_const_range<B1,B2,B3,B4> darts_of_orbit(Dart_const_handle adart) const {
            return Dart_of_orbit_const_range<B1,B2,B3,B4>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
                unsigned int B5>
        Dart_of_orbit_const_range<B1,B2,B3,B4,B5> darts_of_orbit(Dart_const_handle adart) const {
            return Dart_of_orbit_const_range<B1,B2,B3,B4,B5>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
                unsigned int B5,unsigned int B6>
        Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6> darts_of_orbit(Dart_const_handle adart) const {
            return Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
                unsigned int B5,unsigned int B6,unsigned int B7>
        Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7> darts_of_orbit(Dart_const_handle adart) const {
            return Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
                unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8>
        Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7,B8> darts_of_orbit(Dart_const_handle adart) const {
            return Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7,B8>(mGMap, adart);
        }
        
        template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
                unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8,
                unsigned int B9>
        Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7,B8,B9> darts_of_orbit(Dart_const_handle adart) const {
            return Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>(mGMap, adart);
        }

        template<unsigned int i> 
        Dart_of_cell_range<i> darts_of_cell(Dart_handle dh){
            return mGMap.darts_of_cell<i>(dh);
        }
        
        template<unsigned int i,unsigned int dim> 
        Dart_of_cell_range<i> darts_of_cell(Dart_handle dh){
            return mGMap.darts_of_cell<i, dim>(dh);
        }

        template<unsigned int i> 
        Dart_of_cell_const_range<i> darts_of_cell(Dart_const_handle dh) const{
            return mGMap.darts_of_cell<i>(dh);
        }
        
        template<unsigned int i,unsigned int dim> 
        Dart_of_cell_const_range<i> darts_of_cell(Dart_const_handle dh) const{
            return mGMap.darts_of_cell<i, dim>(dh);
        }

        template<unsigned int i,unsigned int j> 
        One_dart_per_incident_cell_range<i, j> one_dart_per_incident_cell(Dart_handle dh){
            return mGMap.one_dart_per_incident_cell<i, j>(dh);
        }
        
        template<unsigned int i,unsigned int j,unsigned int dim> 
        One_dart_per_incident_cell_range<i, j> one_dart_per_incident_cell(Dart_handle dh){
            return mGMap.one_dart_per_incident_cell<i, j, dim>(dh);
        }

        template<unsigned int i,unsigned int j>
        One_dart_per_incident_cell_const_range<i, j> one_dart_per_incident_cell(Dart_const_handle dh) const{
            return mGMap.one_dart_per_incident_cell<i, j>(dh);
        }
        
        template<unsigned int i,unsigned int j,unsigned int dim>
        One_dart_per_incident_cell_const_range<i, j> one_dart_per_incident_cell(Dart_const_handle dh) const{
            return mGMap.one_dart_per_incident_cell<i, j, dim>(dh);
        }

        template<unsigned int i>
        One_dart_per_cell_range<i> one_dart_per_cell(){
            return one_dart_per_cell<i>();
        }
        
        template<unsigned int i,unsigned int dim>
        One_dart_per_cell_range<i> one_dart_per_cell(){
            return one_dart_per_cell<i, dim>();
        }

        template<unsigned int i>
        One_dart_per_cell_const_range<i> one_dart_per_cell() const{
            return one_dart_per_cell<i>();
        }
        
        template<unsigned int i,unsigned int dim>
        One_dart_per_cell_const_range<i> one_dart_per_cell() const{
            return one_dart_per_cell<i, dim>();
        }

        // Modifiers
        Dart_handle create_dart(){
            Dart_handle result =  mGMap.create_dart();
            setSignature(result, true);
            return result;
        }

        void erase_dart(Dart_handle dh){
            mGMap.erase_dart(dh);
        }

        template<unsigned int i,typename T1>
        typename Attribute_handle<i>::type create_attribute(T1 t1){
            return mGMap.create_attribute<i, T1>(t1);
        }

        template <unsigned int i>
        void erase_attribute(typename Attribute_handle<i>::type ah){
            mGMap.erase_attribute<i>(ah);
        }

        template <unsigned int i>
        void set_attribute(Dart_handle dh, typename Attribute_handle<i>::type ah){
            mGMap.set_attribute<i>(dh, ah);
        }

        void clear(){
            mGMap.clear();
        }

        Self& operator= (const Self& cmap){
            return mGMap.operator=(cmap.mGMap);
        }

        void swap(Self& cmap){
            mGMap.swap(cmap.mGMap);
        }

        //Attributes management
        bool are_attributes_automatically_managed() const{
            return mGMap.are_attributes_automatically_managed();
        }

        void set_automatic_attributes_management(bool update_attributes){
            mGMap.set_automatic_attributes_management(update_attributes);
        }

        void correct_invalid_attributes(){
            mGMap.correct_invalid_attributes();
        }

        //Operations
        template <unsigned int i>
        bool is_sewable(Dart_const_handle dh1, Dart_const_handle dh2) const{
            return mGMap.is_sewable<i>(dh1, dh2);
        }

        template <unsigned int i> 
        void sew(Dart_handle dh1,Dart_handle dh2){
            mGMap.sew<i>(dh1, dh2);
            update_signature();
        }

        template <unsigned int i>
        void unsew(Dart_handle dh){
            mGMap.unsew<i>(dh);
        }
        

        template <unsigned int i>
        void link_alpha(Dart_handle dh1, Dart_handle dh2){
            mGMap.link_alpha<i>(dh1, dh2);
            if(i!=0){
                update_signature_vertex(dh1);
            }
        }

        template <unsigned int i>
        void unlink_alpha(Dart_handle dh){
            mGMap.unlink_alpha<i>(dh);
        }

        void reverse_orientation(){
            mGMap.reverse_orientation();
        }

        void reverse_orientation_connected_component(Dart_handle adart){
            mGMap.reverse_orientation_connected_component(adart);
        }

        // Dynamic Onmerge/Onsplit functors
        template<int i>
        boost::function<void(typename Attribute_type< i >::type&,
                            typename Attribute_type< i >::type&)>&
        onsplit_function(){
            return mGMap.onsplit_function<i>();
        }

        template<int i>
        const boost::function<void(typename Attribute_type< i >::type&,
                                    typename Attribute_type< i >::type&)>&
        onsplit_function() const{
            return mGMap.onsplit_function<i>();
        }

        template<int i>
        boost::function<void(typename Attribute_type< i >::type&,
                            typename Attribute_type< i >::type&)>&
        onmerge_function(){
            return mGMap.onmerge_function<i>();
        }

        template<int i>
        const boost::function<void(typename Attribute_type< i >::type&,
                                    typename Attribute_type< i >::type&)>&
        onmerge_function() const{
            return mGMap.onmerge_function<i>();
        }

        //Boolean Marks
        size_type get_new_mark() const{
            return mGMap.get_new_mark();
        }

        bool is_reserved(size_type m) const{
            return mGMap.is_reserved(m);
        }

        bool is_marked(Dart_const_handle dh, size_type m) const{
            return mGMap.is_marked(dh, m);
        }

        void mark(Dart_const_handle dh, size_type m) const{
            mGMap.mark(dh, m);
        }

        void unmark(Dart_const_handle dh, size_type m) const{
            mGMap.unmark(dh, m);
        }

        void negate_mark(size_type m) const{
            mGMap.negate_mark(m);
        }

        void unmark_all(size_type m) const{
            mGMap.unmark_all(m);
        }

        size_type number_of_marked_darts(size_type m) const{
            return mGMap.number_of_marked_darts(m);
        }

        size_type number_of_unmarked_darts(size_type m) const{
            return mGMap.number_of_unmarked_darts(m);
        }

        void free_mark(size_type m) const{
            mGMap.free_mark(m);
        }
        
        //Construction Operations
        Dart_handle make_combinatorial_hexahedron(){
            Dart_handle result = CGAL::make_combinatorial_hexahedron<GMap>(mGMap);
            update_signature();
            return result;
        }
        
        Dart_handle make_combinatorial_polygon(unsigned int lg){
            Dart_handle result = CGAL::make_combinatorial_polygon(mGMap, lg);
            update_signature();
            return result;
        }
        
        Dart_handle make_combinatorial_tetrahedron(){
            Dart_handle result = CGAL::make_combinatorial_tetrahedron<GMap>(mGMap);
            update_signature();
            return result;
        }
        
        Dart_handle make_edge(){
            return mGMap.make_edge();
        }
        
        //Modification Operations
        Dart_handle insert_cell_0_in_cell_1 (Dart_handle dh){
            return mGMap.insert_cell_0_in_cell_1(dh);
        }
 
        Dart_handle insert_cell_0_in_cell_2 (Dart_handle dh){
            return mGMap.insert_cell_0_in_cell_2(dh);
        }
        
        Dart_handle insert_cell_1_in_cell_2 (Dart_handle dh1, Dart_handle dh2){
            return mGMap.insert_cell_1_in_cell_2(dh1, dh2);
        }
        
        template<class InputIterator >
        Dart_handle insert_cell_2_in_cell_3 (InputIterator afirst, InputIterator alast){
            return mGMap.insert_cell_2_in_cell_3<InputIterator>(afirst, alast);
        }
        
        Dart_handle insert_dangling_cell_1_in_cell_2 (Dart_handle dh){
            return mGMap.insert_dangling_cell_1_in_cell_2(dh);
        }
        
        bool is_insertable_cell_1_in_cell_2 (Dart_const_handle dh1, Dart_const_handle dh2){
            return mGMap.is_insertable_cell_1_in_cell_2(dh1, dh2);
        }
        
        template<class InputIterator >
        bool is_insertable_cell_2_in_cell_3 (InputIterator afirst, InputIterator alast){
            return mGMap.is_insertable_cell_2_in_cell_3<InputIterator>(afirst, alast);
        }
        
        template<unsigned int i>
        bool is_removable (Dart_const_handle dh){
            return mGMap.is_removable<i>(dh);
        }
        
        template<unsigned int i>
        size_type remove_cell (Dart_handle dh){
            return mGMap.remove_cell<i>(dh);
        }
        
        //additional access    
        Halfedge_handle vertex_next (Halfedge_handle he){
            if(dart_signature(he.dart)){
                CGAL_precondition(!is_free<2>(he.dart) && !is_free<1>(alpha<2>(he.dart)));
                return halfedge(alpha<2, 1>(he.dart));
            }else{ 
                CGAL_precondition(!is_free<1>(he.dart));
                return halfedge(alpha<1>(he.dart));
            }
        }
        
        Halfedge_handle halfedge_opposite (Halfedge_handle he){
            CGAL_precondition(!is_free<0>(he.dart));
            return halfedge(alpha<0>(he.dart));
        }
        
        Dart_handle dart(Halfedge_handle dh){
            return dh.dart;
        }
        
        Dart_handle dart(Halfedge_handle dh, bool s){
            Dart_handle d(dh->dart);
            if(dart_signature(d)!=s){
                d = alpha<2>(d);
            }
            return d;
        }
        
        Halfedge_handle halfedge (Dart_handle d){
            return Halfedge_handle(*this, d);
        }
        
        bool halfedge_signature (Halfedge_handle he){
            CGAL_precondition(!is_free<0>(he.dart));
            return dart_signature(he.dart) != dart_signature(alpha<0>(he.dart));
        }
        
        bool dart_signature (Dart_const_handle d) const{
            return is_marked(d, mSignature);
        }
        
        Halfedge_handle halfedge_next(Halfedge_handle he){
            return halfedge_opposite(vertex_next(he));
        }
        
        //addtition modifier
        void update_signature(){
            size_type mark = get_new_mark();
            
            std::queue<Dart_handle> queue;
            
            Dart_handle start_vertex = darts().begin();
            
            queue.push(start_vertex);
            
            update_signature_vertex(start_vertex, true, mark);
            
            while(!queue.empty()){
                Dart_handle current_vertex = queue.front();
                queue.pop();
                
                for(typename Dart_of_orbit_range<1,2>::iterator it = darts_of_orbit<1,2>(current_vertex).begin();
                    it!=darts_of_orbit<1,2>(current_vertex).end();
                    ++it)
                {
                    if(!is_marked(alpha<0>(it), mark)){
                        Dart_handle neighbour = alpha<0>(it);
                        update_signature_vertex(neighbour, !dart_signature(it), mark);
                        
                        queue.push(neighbour);
                    }               
                }
            }
            
            for(typename Dart_range::iterator it = darts().begin(); it!=darts().end(); ++it){
                assert(is_marked(it, mark));
            }
        }
        
        //Path
        Path_handle create_path(){
            Path_handle p = mPaths.insert(Path_());
            p->mHandle = p;
            return p;
        }
        
        void erase_path (Path_handle p){
            mPaths.erase(p);
        }

    private:
        GMap mGMap;
        size_type mSignature;
        Path_container mPaths;
        
        void update_signature_vertex(Dart_handle d, bool sign, size_type mark_index){
            CGAL_precondition(d!=Dart_handle());
            
            Dart_handle current = d;
            bool loop2 = false;
            
            do{
                mark(current, mark_index);
                setSignature(current, sign);
                
                if(is_free<1>(current)){
                    loop2 = true;
                    break;
                }
                
                current = alpha<1>(current);
                
                mark(current, mark_index);
                setSignature(current, !sign);
                
                if(is_free<2>(current)){
                    loop2 = true;
                    break;
                }
                
                current = alpha<2>(current);
                
            }while(current!=d);
            
            if(loop2) {  
                current = d;
                do{
                    mark(current, mark_index);
                    setSignature(current, sign);
                    
                    if(is_free<2>(current)){
                        return;
                    }
                    
                    current = alpha<2>(current);
                    
                    mark(current, mark_index);
                    setSignature(current, !sign);
                    
                    if(is_free<1>(current)){
                        return;
                    }
                    
                    current = alpha<1>(current);
                    
                    
                }while(loop2);
            }
        }
        
        void update_signature_vertex(Dart_handle d){
            CGAL_precondition(d!=Dart_handle());
            
            Dart_handle current = d;
            bool loop2 = false;
            
            do{
                setSignature(current, true);
                
                if(is_free<1>(current)){
                    loop2 = true;
                    break;
                }
                
                current = alpha<1>(current);
                
                setSignature(current, false);
                
                if(is_free<2>(current)){
                    loop2 = true;
                    break;
                }
                
                current = alpha<2>(current);
                
            }while(current!=d);
            
            if(loop2) {  
                current = d;
                do{
                    setSignature(current, true);
                    
                    if(is_free<2>(current)){
                        return;
                    }
                    
                    current = alpha<2>(current);
                    
                    setSignature(current, false);
                    
                    if(is_free<1>(current)){
                        return;
                    }
                    
                    current = alpha<1>(current);
                    
                    
                }while(loop2);
            }
        }
        
        void setSignature(Dart_handle d, bool s){
            if(s){
                mark(d, mSignature);
            }else{
                unmark(d, mSignature);
            }
        }
    };
}

#endif //TOPOLOGICAL_SURFACE_H