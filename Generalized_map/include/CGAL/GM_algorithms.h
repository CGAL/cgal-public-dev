/* SRC_GM_algorithms.h
 *
 *  Created on: 10 avr. 2017
 *
 */

#ifndef SRC_GM_algorithms_H_
#define SRC_GM_algorithms_H_

#include <math.h>
#include<CGAL/Generalized_map.h>
#include <CGAL/Cell_attribute.h>
#include <boost/unordered_map.hpp>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include "My_cell_attribute_with_point.h"
#include <CGAL/Cartesian.h>
#include <CGAL/tags.h>
#include <algorithm>

#include <vector>
#include <bitset>
#include <string>
#include <utility>
#include <string.h>
#include <iostream>
#include <cstdlib>
#include <queue>

#include <boost/config.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/properties.hpp>

#include <boost/property_map/property_map.hpp>
#include <typeinfo>

/**
 *
 * Author : Christophe Riedinger
 * a class for computing the shortest non contractible cycle
 * on a Combinatorial Map
 * Input: the Generalized Map
 * Ouput: the shortest cycle non contractible
 */

namespace CGAL {



struct Average_functor
{
  template<class CellAttribute>
  void operator()(CellAttribute& ca1, CellAttribute& ca2)
  { ca1.info()=(ca1.info()+ ca2.info())/2; }
};

/*
 * The attributes used for the map
 *
 */

  struct Item //: public  //CGAL::Generic_map_min_items
  {
    template<class K>
    struct Dart_wrapper
    {
      typedef CGAL::My_cell_attribute_with_point< K, int, CGAL::Tag_true,Null_functor,Null_functor>
      //Average_functor >
      Vertex_attribute;
      typedef CGAL::Cell_attribute<K,int> Edge_attribute;
      typedef CGAL::Cell_attribute<K,int> Face_attribute;
      typedef CGAL::cpp11::tuple<Vertex_attribute,Edge_attribute,Vertex_attribute> Attributes;
    };
  };

  typedef CGAL::Linear_cell_complex_traits<3, CGAL::Exact_predicates_inexact_constructions_kernel> Traits;

  typedef CGAL::Linear_cell_complex_for_generalized_map<2,3,Traits,Item> LCC_3;

  typedef CGAL::Generalized_map<2> GMap2I;

  typedef LCC_3::Point                      Point;
  typedef LCC_3::FT   FT;

  template<class K=LCC_3>
  struct Edge_and_next_edge {
  public:
    typedef typename K::Dart_handle Dart_handle;
    Dart_handle de1;
    Dart_handle de2;
    Edge_and_next_edge<K>* dne1=0;
    Edge_and_next_edge<K>* dpe1;
    int degree=0;
  };

  /*
   *
   * The algorihtm class
   *
   */

  template<class K=LCC_3>
  class GM_algorithms
  {
  public :
    K map;
    std::map<int,int> boundary;
    typename K::size_type mark0;
    int edgesnumber=0;


    typedef typename K::Dart_handle              Dart_handle;

    typedef Dart_handle Dart_h;

    std::vector<int> V;
    typedef typename K::template Attribute_handle<0>::type Attrib_h0;
    typedef typename K::template Attribute_handle<2>::type Attrib_h20;
    std::vector<Attrib_h0> verticesvector;
    std::vector<Attrib_h20> facesvector;

    typedef CGAL::My_cell_attribute_with_point< K, int, CGAL::Tag_true,
                                                Average_functor > vattrib;

    /*
     *
     * constructor with a map
     *
     */

    GM_algorithms(K& g) {

      map=g;
    }

    /*
     *
     * default constructor
     *
     */

    GM_algorithms() {

    }

    void set_attributes() {

      typedef typename K ::size_type size_type;
      int i=0;
      int i2=0;
      int i1=0;
      typedef typename K::Vertex_attribute vertex;
      typedef typename K::template Attribute_handle<0>::type Attrib_h;
      typedef typename K::template Attribute_handle<2>::type Attrib_h2;

      for (typename K::template Attribute_range<0>::type::iterator it=map.template attributes<0>().begin();
           it!=map.template attributes<0>().end();it++){
        map.template erase_attribute<0>(it);
      }
      for (typename K::template Attribute_range<1>::type::iterator it=map.template attributes<1>().begin();
           it!=map.template attributes<1>().end();it++){
        map.template erase_attribute<1>(it);
      }
      for (typename K::template Attribute_range<2>::type::iterator it=map.template attributes<2>().begin();
           it!=map.template attributes<2>().end();it++){
        map.template erase_attribute<2>(it);
      }
      std::cout << "tototototo\n";

      verticesvector.clear();
      facesvector.clear();
      int first=1;

      for (typename K::template One_dart_per_cell_range<0>::
             iterator it=map.template one_dart_per_cell<0>().begin();
           it!=map.template one_dart_per_cell<0>().end();++it)
      {

        int index=0;

        first=0;

        typename K::template Attribute_handle<0>::type ah;
        ah=map.template create_vertex_attribute(map.template point_of_vertex_attribute(map.template vertex_attribute(it)));//map.template vertex_attribute(it);
        map.template set_vertex_attribute(it,ah );
        map.template info<0>(it)=i;

        verticesvector.push_back(ah);
        ah->template prev=-1;
        ah->template distance=0;
        i++;

      }

      first=1;
      for (typename K::template One_dart_per_cell_range<1>::
             iterator it=map.template one_dart_per_cell<1>().begin();
           it!=map.template one_dart_per_cell<1>().end();++it)
      {

        int index=0;

        first=0;

        typename K::template Attribute_handle<1>::type ah;
        ah=map.template create_attribute<1>();//map.template vertex_attribute(it);
        map.template set_attribute<1>(it,ah );
        map.template info<1>(it)=i1;
        i1++;
      }

      first=1;
      for (typename K:: template One_dart_per_cell_range<2>::
             iterator it=map.template one_dart_per_cell<2>().begin();
           it!=map.template one_dart_per_cell<2>().end();++it)
      {

        first=0;
        typename K::template Attribute_handle<2>::type ah=//map.template attribute<2>(it);
          map.template create_vertex_attribute(map.template point_of_vertex_attribute(map.template vertex_attribute(it)));//map.template vertex_attribute(it);

        //std::cout << (map.template darts_of_cell<2>(it).size()) << " facenumber\n";
        map.template set_attribute<2>(it,ah);
        map.template info<2>(it)=i2;
        ah->boundary=0;
        facesvector.push_back(ah);
        i2++;

      }


      for (std::map<int,int>::iterator it=boundary.begin();it!=boundary.end();it++) {

        facesvector[it->second]->boundary=1;

      }

      std::cout << "tototototo\n";


    }


    /*
     *
     * shortest cycle computation
     *
     */

    void compute_shortest_cycle_non_contractible(int geometric,std::map<int,int> boundary0,
                                                 std::vector<int>& W,int set_attribs,
                                                 int savefile) {//std::map<int,std::map<int,double> >& W

      //verticesvector=std::vector<Attrib_h0>();
      //facesvector=std::vector<Attrib_h20>();
      boundary=boundary0;
      std::cout << boundary.size() << " boundary \n";
      typedef typename K ::size_type size_type;
      int i=0;
      int i2=0;
      typedef typename K::Vertex_attribute vertex;
      typedef typename K::template Attribute_handle<0>::type Attrib_h;
      typedef typename K::template Attribute_handle<2>::type Attrib_h2;

      int first=1;

      if (set_attribs) {
        for (typename K::template One_dart_per_cell_range<0>::
               iterator it=map.template one_dart_per_cell<0>().begin();
             it!=map.template one_dart_per_cell<0>().end();++it)
        {

          int index=0;

          first=0;
          std::cout << i <<"\n";
          typename K::template Attribute_handle<0>::type ah;
          //ah=map.template create_vertex_attribute(//map.template vertex_attribute(it);
          ah=map.template create_vertex_attribute(map.template point_of_vertex_attribute(map.template vertex_attribute(it)));//map.template vertex_attribute(it);
          map.template set_vertex_attribute(it,ah );
          map.template info<0>(it)=i;

          verticesvector.push_back(ah);
          ah->template prev=-1;
          ah->template distance=0;
          i++;

        }

        int i1=0;
        for (typename K:: template One_dart_per_cell_range<1>::
               iterator it=map.template one_dart_per_cell<1>().begin();
             it!=map.template one_dart_per_cell<1>().end();++it)
        {

          first=0;
          typename K::template Attribute_handle<1>::type ah=//map.template attribute<2>(it);
            map.template create_attribute<1>(i1);

          std::cout << (map.template darts_of_cell<1>(it).size()) << " facenumber\n";
          map.template set_attribute<1>(it,ah);
          map.template info<1>(it)=i1;
          i1++;

        }

        first=1;
        for (typename K:: template One_dart_per_cell_range<2>::
               iterator it=map.template one_dart_per_cell<2>().begin();
             it!=map.template one_dart_per_cell<2>().end();++it)
        {

          first=0;
          typename K::template Attribute_handle<2>::type ah=//map.template attribute<2>(it);
            map.template create_attribute<2>(map.template point_of_vertex_attribute(map.template vertex_attribute(it)));//map.template vertex_attribute(it);

          std::cout << (map.template darts_of_cell<2>(it).size()) << " facenumber\n";
          map.template set_attribute<2>(it,ah);
          map.template info<2>(it)=i2;
          ah->boundary=0;
          facesvector.push_back(ah);
          i2++;

        }


        for (std::map<int,int>::iterator it=boundary.begin();it!=boundary.end();it++) {

          facesvector[it->second]->boundary=1;

        }
      }

      map.template display_characteristics(std::cout);
      std::cout<<", valid="<< map.template is_valid();



      std::cout << "toto\n";


      Graph g;
      std::vector<Vertex> vertices=create_CGAL_Graph2(g,geometric,W);



      boost::unordered_map<int,int> border;
      typedef boost::graph_traits < Graph >::vertex_descriptor Vertex;
      std::vector<Vertex> vertexvector;
      double mindist=10000;
      int min1=-1;

      int min2=-1;

      size_type amark2= map.get_new_mark();
      Dart_handle source;
      int sourceindex=-1;

      time_t t10 = time(0);
      std::cout << verticesvector.size() << " toto\n";
      for (int k=0;k<verticesvector.size();k++) {
        Dart_handle dh=verticesvector[k]->template dart();
        source=dh;
        for (int k2=0;k2<verticesvector.size();k2++) {
          verticesvector[k2]->template prev=-1;
          verticesvector[k2]->template distance=10000;
          verticesvector[k2]->template boundary=false;
        }
        map.template free_mark(amark2);
        amark2= map.get_new_mark();
        CGAL_Dijkstra(dh,amark2,g,vertices);


        boost::unordered_map<int,int> valids=remove_face_graph_sommetsF(border,amark2);
        std::cout << k << "K\n";
        double* pmin=compute_shortest_cycle_with_dual_graph3(source,valids,amark2,geometric);

        if (pmin[0]<mindist) {
          mindist=pmin[0];
          min1=(int)pmin[1];
          min2=(int)pmin[2];
          sourceindex=k;//map.template info<0>(source);
          V.clear();
          for (int k2=0;k2<verticesvector.size();k2++) {
            int i=verticesvector[k2]->prev;
            V.push_back(i);
          }

        }
      }


      if (0!=mindist&&savefile) {
        std::ofstream ifile("model2.off");
        CGAL::write_off<CGAL::LCC_3>(map,ifile );
        if(min1!=-1&&min2!=-1) {
          std::cout << sourceindex <<" 01\n";
          std::list<int> path2=compute_path((int)min1,
                                            sourceindex,
                                            V,1);
          std::cout << sourceindex <<" 02\n";
          std::list<int> path3=compute_path((int)min2,
                                            sourceindex,
                                            V,0);
          std::cout << sourceindex <<" 03\n";
          for (std::list<int>::iterator itp=path3.begin();itp!=path3.end();itp++) {
            path2.push_back(*itp);
          }
          std::cout << sourceindex <<" 04\n";
          std::vector<int> path0;
          for (std::list<int>::iterator itp=path2.begin();itp!=path2.end();itp++) {
            path0.push_back(*itp);
            std::cout << *itp << "\n";
          }
          path0.erase(path0.begin());

          FILE* fout2=fopen("minimalpathtest.vect","w");
          fprintf(fout2,"VECT\n");
          fprintf(fout2,"%d %d %d\n",1,(int)(path0.size()),1);
          fprintf(fout2,"%d\n",(int)(path0.size()));
          fprintf(fout2,"%d\n",1);
          for (int k=0;k<path0.size();k++) {
            float x=(float)map.template point_of_vertex_attribute(verticesvector[path0[k]]).x();
            float y=(float)map.template point_of_vertex_attribute(verticesvector[path0[k]]).y();
            float z=(float)map.template point_of_vertex_attribute(verticesvector[path0[k]]).z();
            fprintf(fout2,"%6.2f %6.2f %6.2f\n",x,
                    y,
                    z);
          }
          fprintf(fout2,"%d %d %d %d\n",0,1,0,1);
          fclose(fout2);
        }
      }
      std::cout << mindist << " mindist\n";
      time_t t20 = time(0);
      std::cout << ((t20-t10)) << "time ellapsed\n";
    }

    /**
     *
     * fills V1 the list of Sommet of the shortest non contractibel cycle originating from a source
     * source : a dart around the source vertex
     * V1 : a 2 element vertex containg the edge extremity of the shortest cycle
     * graphdual : the faces graph
     *
     */
    typedef typename K::template Attribute_handle<0>::type Attrib_h;
    typedef typename K::template Attribute_handle<2>::type Attrib_h2;
    double* compute_shortest_cycle_with_dual_graph3(Dart_h& source,boost::unordered_map<int,int>& valids,
                                                    typename K::size_type amark2,int geometric) {
      std::vector<Attrib_h>& mapS=verticesvector;
      int k=0;
      int mindist[2]={0,0};


      double pmin=10000;
      double* minpath=(double*)malloc(3*sizeof(double));
      minpath[0]=10000;
      minpath[1]=0;
      minpath[2]=0;

      bool valid=false;
      double p1size=5000;
      double p2size=5000;

      for (typename K::Dart_range::iterator it=map.template darts().begin();
           it!=map.template darts().end();++it)
      {

        //boost::unordered_map<int,int> neigh=find_neighbors_vertex(it);
        std::vector<int> neigh=find_neighbors_vertex(it);

        //for (std::vector<int>::iterator itn=neigh.begin();
        //		itn!=neigh.end();
        //		itn++) {
        int index= map.template info<0>(it);

        Dart_handle d2=map.template alpha<0>(it);
        int index2=map.template info<0>(d2);
        //if (index!=map.template info<0>(source)){
        //if (index2!=
        //	map.template info<0>(source)){
        if ((valids.find(map.template info<2>(d2))->second==1
             )&&
            valids.find(map.template info<2>(map.template alpha<2>(d2)))->second==1&&
            (!(verticesvector[index]->prev==index2)&&
             !(verticesvector[index2]->prev==index))){//&&
          //index!=index2) {

          p1size=verticesvector[index]->template distance;//;//p1.size();
          p2size=verticesvector[index2]
            ->template distance;
          double Wa=1.0;

          if (geometric) {

            Point p1=map.template point_of_vertex_attribute(verticesvector[index]);
            Point p2=map.template point_of_vertex_attribute(verticesvector[index2]);
            Wa=(double)(std::sqrt(std::pow(p2.x()-p1.x(),2)+
                                  std::pow(p2.y()-p1.y(),2)+
                                  std::pow(p2.z()-p1.z(),2)));
            //weight0=W.at(itn).at(*neighit);
          }

          //std::cout << "boucle\n";
          if (p1size+p2size+Wa<pmin){//&&p1.size()+p2.size()+1>4){//&&p1size+p2size!=2&&p1size+p2size>4)
            pmin=p1size+p2size+Wa;
            minpath[0]=pmin;minpath[1]=(double)index;
            minpath[2]=(double)index2;

          }
        }
        //}

        //}
        //}



      }
      std::cout << minpath[0] << " pmin\n";

      return minpath;
    }


    std::list<int> compute_path(int beginindex,int sourceindex,
                                std::vector<int>& V,int forward){

      std::list<int> path;
      int index= beginindex;
      path.push_back(index);
      std::cout << V[index] << "\n";

      while (true) {
        index=V[index];
        if (index==-1)
          break;
        std::cout << V[index] << "\n";
        if (forward)
          path.push_back(index);
        else
          path.push_front(index);


      }
      return path;
    }


    void display(LCC_3::Dart_handle d)
    {
      std::cout << "display";
      for (typename K::template One_dart_per_incident_cell_range<0,2>::
             iterator it=map.template one_dart_per_incident_cell<0,2>
             (d).begin(),
             itend=map.template one_dart_per_incident_cell<0,2>
             (d).end();
           it!=itend; ++it)
      {
        std::cout << map.point(it) << "; ";
      }
      std::cout<<std::endl;
    }

    typedef double Weight;
    typedef boost::property<boost::edge_weight_t, Weight> WeightProperty;
    typedef boost::property<boost::vertex_name_t, std::string> NameProperty;

    typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS,
                                    NameProperty, WeightProperty > Graph;

    typedef boost::graph_traits < Graph >::vertex_descriptor Vertex;

    typedef boost::property_map < Graph, boost::vertex_index_t >::type IndexMap;
    typedef boost::property_map < Graph, boost::vertex_name_t >::type NameMap;

    typedef boost::iterator_property_map < Vertex*, IndexMap, Vertex, Vertex& > PredecessorMap;
    typedef boost::iterator_property_map < Weight*, IndexMap, Weight, Weight& > DistanceMap;


    /*boost::unordered_map<int,int> find_neighbors_vertex(typename K::Dart_handle& d) {
      boost::unordered_map<int,int> neighbors;
      Dart_handle d0=map.template alpha<0>(d);
      int i0=map.template info<0>(d0);
      neighbors.insert(std::make_pair(i0,i0));
      bool first=true;
      while (true) {
      d=map.template alpha<2>(d);
      d=map.template alpha<1>(d);
      Dart_handle d2=map.template alpha<0>(d);
      int i=map.template info<0>(d2);
      neighbors.insert(std::make_pair(i,i));
      if (i==i0&&!first)
      return neighbors;
      first=false;
      }
      return neighbors;
      }*/

    std::vector<int> find_neighbors_vertex(typename K::Dart_handle& d) {
      //boost::unordered_map<int,int> neighbors;
      std::vector<int> neighbors;
      std::vector<Dart_handle> essai;
      Dart_handle d0=map.template alpha<0>(d);
      int i0=map.template info<0>(d0);
      neighbors.push_back(i0);
      bool first=true;
      while (true) {
        d=map.template alpha<2>(d);
        d=map.template alpha<1>(d);
        Dart_handle d2=map.template alpha<0>(d);
        essai.push_back(d);
        int i=map.template info<0>(d2);
        //if (i==i0&&!first) {
        if (d0==d2) {
          /*for (int k=0;k<essai.size();k++) {
            printf("%p\n",&(essai[k]));
            }
            std::cout << "essai\n";*/
          return neighbors;
        }
        neighbors.push_back(i);
        first=false;
      }

      return neighbors;
    }

    boost::unordered_map<int,int> find_neighbors_face(Dart_handle& d,
                                                      typename K::size_type& amark2,std::vector<int>& valids) {
      boost::unordered_map<int,int> neighbors;
      int info0=map.template info<0>(d);
      Dart_handle d0=map.template alpha<2>(d);
      int i0=map.template info<2>(d0);
      if(!map.template is_marked(d,amark2)&&valids[i0]==1)
        neighbors.insert(std::make_pair(i0,i0));
      bool first=true;
      while (true) {
        d=map.template alpha<0>(d);
        d=map.template alpha<1>(d);
        int info1=map.template info<0>(d);
        Dart_handle d2=map.template alpha<2>(d);
        int i=map.template info<2>(d2);

        //if (info0==info1)
        if (d2==d0)
          return neighbors;
        if(!map.template is_marked(d,amark2)&&valids[i]==1)
          neighbors.insert(std::make_pair(i,i));
        first=false;
      }
      return neighbors;
    }


    /**
     * Given the faces graph, removes recursively the appropriate faces in this graph
     * S2 : the face graph
     * lremoved : a list of the of the removed faces
     * mapremoved : unused
     * source : unused
     */

    boost::unordered_map<int,int> remove_face_graph_sommetsF(
                                                             boost::unordered_map<int,int>& border,
                                                             typename K::size_type amark2) {
      std::map<int,int> l1;
      std::vector<int>  valids;
      int bc=0;
      for (int k=0;k< facesvector.size();k++) {
        valids.push_back(1);

        l1.insert(std::make_pair(k,k));

      }
      int removed0=0;

      while (true) {
        int removed=0;


        int lsize=l1.size();
        std::map<int,int> l2;
        for (std::map<int,int>::iterator it0=l1.begin();it0!=l1.end();it0++) {

          int it=it0->second;

          if (valids[it]==1) {
            int sommetsvalid=0;
            Dart_handle d1=facesvector[it]->template dart();

            bool border=is_boundary(d1);

            boost::unordered_map<int,int>  neighbors=
              find_neighbors_face(d1,amark2,valids);

            if (neighbors.size()==1&&!border) {
              int facenumber = map.template info<2>(d1);
              valids[facenumber]=-1;
              removed++;
              Dart_handle dart=facesvector[neighbors.begin()->second]->dart();

              l2.insert(std::make_pair(neighbors.begin()->second,
                                       neighbors.begin()->second));
            } else if (neighbors.size()==0&&!border){
              int facenumber = it;
              valids[facenumber]=-1;
              removed++;
            } else {

              l2.insert(std::make_pair(it,it));

            }


          }
        }


        l1=l2;

        if (removed==0)
          break;
      }

      boost::unordered_map<int,int> valids2;
      for (int k=0;k<valids.size();k++) {

        valids2.insert(std::make_pair(k,valids[k]));
      }
      return valids2;

    }

    bool is_boundary(Dart_handle& d) {

      if(facesvector[map.template info<2>(d)]->boundary==1) {
        return true;
      }

      return false;
    }


    /* wraps the boost implementation of the Dijkstra algorithm
     * source : a dart of the GMap around the source vertex
     * amark2 : the marker for free edges
     * g: an empty boost graph
     * vertices : an empty vector of boost Vertex
     *
     */

    std::vector<Vertex> create_CGAL_Graph(Graph& g,int geometric,std::map<int,std::map<int,double> >& W) {

      std::vector<Vertex> vertices;//(map.template count_all_cells()[0]);
      std::vector<Vertex> vertices2;
      boost::unordered_map<int,boost::unordered_map<int,int> > edgesdone;

      boost::unordered_map<int,int> vertexdone;
      for (typename K::template One_dart_per_cell_range<0>::
             iterator it=map.template one_dart_per_cell<0>().begin();
           it!=map.template one_dart_per_cell<0>().end();++it)
      {

        char buffer [6];
        typedef typename K::Dart_handle itd;

        int a=(map.template info<0>(it));
        if (vertexdone.find(a)==vertexdone.end()) {
          sprintf (buffer, "%d", a);
          Vertex v0 = boost::add_vertex(std::string(buffer), g);
          boost::unordered_map<int,int> V;
          edgesdone.insert(std::make_pair(a,V));
          vertexdone.insert(std::make_pair(a,a));
          vertices.push_back(v0);
          std::cout << buffer << "vertex number\n";
        }
      }
      for (int k=0;k<vertices.size();k++){
        vertices2.push_back(vertices[vertices.size()-1-k]);
      }
      for (typename K::template One_dart_per_cell_range<0>::
             iterator it=map.template one_dart_per_cell<0>().begin();
           it!=map.template one_dart_per_cell<0>().end();++it)
      {

        int itn=(int)(map.template info<0>(it));
        ///boost::unordered_map<int,int> neigh=find_neighbors_vertex(it);
        std::vector<int> neigh=find_neighbors_vertex(it);

        for (std::vector<int>::iterator neighit=neigh.begin();
             neighit!=neigh.end();neighit++){

          bool dontdo=false;
          Weight weight0 = 1.0;
          if (geometric) {

            Point p1=map.template point_of_vertex_attribute(verticesvector[itn]);
            Point p2=map.template point_of_vertex_attribute(verticesvector[*neighit]);
            weight0=(double)(std::sqrt(std::pow(p2.x()-p1.x(),2)+
                                       std::pow(p2.y()-p1.y(),2)+
                                       std::pow(p2.z()-p1.z(),2)));
            //weight0=W.at(itn).at(*neighit);

          }
          for (int k=0;k<edgesdone.at(itn).size();k++) {
            if (edgesdone.at(itn).find(*neighit)!=edgesdone.at(itn).end()||
                edgesdone.at(*neighit).find(itn)!=edgesdone.at(*neighit).end()) {
              dontdo=true;
            }
          }
          if (dontdo==false) {

            boost::add_edge(vertices[itn], vertices[*neighit], weight0, g);
            edgesdone.at(itn).insert(std::make_pair(*neighit,*neighit));
            edgesdone.at(*neighit).insert(std::make_pair(*neighit,*neighit));

          }
        }
      }
      std::cout << vertices.size() << " " << verticesvector.size() << " vertex vertex";
      return vertices;
    }

    /* wraps the boost implementation of the Dijkstra algorithm
     * source : a dart of the GMap around the source vertex
     * amark2 : the marker for free edges
     * g: an empty boost graph
     * vertices : an empty vector of boost Vertex
     *
     */

    std::vector<Vertex> create_CGAL_Graph2(Graph& g,int geometric,std::vector<int>& W) {

      std::vector<Vertex> vertices;//(map.template count_all_cells()[0]);
      std::vector<Vertex> vertices2;
      boost::unordered_map<int,boost::unordered_map<int,int> > edgesdone;

      boost::unordered_map<int,int> vertexdone;
      for (typename K::template One_dart_per_cell_range<0>::
             iterator it=map.template one_dart_per_cell<0>().begin();
           it!=map.template one_dart_per_cell<0>().end();++it)
      {

        char buffer [6];
        typedef typename K::Dart_handle itd;

        int a=(map.template info<0>(it));
        if (vertexdone.find(a)==vertexdone.end()) {
          sprintf (buffer, "%d", a);
          Vertex v0 = boost::add_vertex(std::string(buffer), g);
          boost::unordered_map<int,int> V;
          edgesdone.insert(std::make_pair(a,V));
          vertexdone.insert(std::make_pair(a,a));
          vertices.push_back(v0);
          std::cout << buffer << "vertex number\n";
        }
      }
      for (int k=0;k<vertices.size();k++){
        vertices2.push_back(vertices[vertices.size()-1-k]);
      }
      for (typename K::template One_dart_per_cell_range<0>::
             iterator it=map.template one_dart_per_cell<0>().begin();
           it!=map.template one_dart_per_cell<0>().end();++it)
      {

        int itn=(int)(map.template info<0>(it));
        ///boost::unordered_map<int,int> neigh=find_neighbors_vertex(it);
        std::vector<int> neigh=find_neighbors_vertex(it);

        for (std::vector<int>::iterator neighit=neigh.begin();
             neighit!=neigh.end();neighit++){

          bool dontdo=false;
          Weight weight0 = 1.0;
          if (geometric) {

            Point p1=map.template point_of_vertex_attribute(verticesvector[itn]);
            Point p2=map.template point_of_vertex_attribute(verticesvector[*neighit]);
            weight0=(double)(std::sqrt(std::pow(p2.x()-p1.x(),2)+
                                       std::pow(p2.y()-p1.y(),2)+
                                       std::pow(p2.z()-p1.z(),2)));
            //weight0=W.at(itn).at(*neighit);

          }

          if (W.size()==map.template count_all_cells()[1]) {
            Dart_handle d1=verticesvector[itn]->dart();
            Dart_handle d2=verticesvector[*neighit]->dart();
            int edgeindice=map.template info<1>(d1);
            weight0=W[edgeindice];
          }

          for (int k=0;k<edgesdone.at(itn).size();k++) {
            if (edgesdone.at(itn).find(*neighit)!=edgesdone.at(itn).end()||
                edgesdone.at(*neighit).find(itn)!=edgesdone.at(*neighit).end()) {
              dontdo=true;
            }
          }
          if (dontdo==false) {

            boost::add_edge(vertices[itn], vertices[*neighit], weight0, g);

          }
        }

        for (std::vector<int>::iterator neighit=neigh.begin();
             neighit!=neigh.end();neighit++){
          edgesdone.at(itn).insert(std::make_pair(*neighit,*neighit));
          edgesdone.at(*neighit).insert(std::make_pair(itn,itn));
        }
      }
      std::pair<boost::graph_traits<Graph>::edge_iterator, boost::graph_traits<Graph>::edge_iterator> edgeIteratorRange = boost::edges(g);
      int cc=0;
      for(boost::graph_traits<Graph>::edge_iterator edgeIterator = edgeIteratorRange.first; edgeIterator != edgeIteratorRange.second; ++edgeIterator)
      {
        std::cout << (*edgeIterator).m_source << " " <<
          (*edgeIterator).m_target << " " << (cc++) <<" edges\n";
      }

      std::cout << vertices.size() << " " << verticesvector.size() << " vertex vertex";
      return vertices;
    }

    /* wraps the boost implementation of the Dijkstra algorithm
     * source : a dart of the GMap around the source vertex
     * amark2 : the marker for free edges
     * g: an empty boost graph
     * vertices : an empty vector of boost Vertex
     *
     */
    void CGAL_Dijkstra(Dart_handle& source,typename K::size_type& amark2,Graph& g,
                       std::vector<Vertex>& vertices
                       ) {
      typedef typename K::Vertex_attribute vertex;
      std::vector<Vertex> predecessors(boost::num_vertices(g)); // To store parents

      std::vector<Weight> distances(boost::num_vertices(g)); // To store distances

      IndexMap indexMap = boost::get(boost::vertex_index, g);

      PredecessorMap predecessorMap(&predecessors[0], indexMap);
      DistanceMap distanceMap(&distances[0], indexMap);
      boost::dijkstra_shortest_paths(g, vertices[map.template info<0>(source)]
                                     , boost::distance_map(distanceMap).predecessor_map(predecessorMap));

      NameMap nameMap = boost::get(boost::vertex_name, g);
      int k=0;

      for (int k=0;k<vertices.size();k++) {

        Vertex v=vertices[k];
        if (k!=map.template info<0>(source))
          verticesvector[k]->prev=atoi(nameMap[predecessorMap[vertices[v]]].c_str());
        else
          verticesvector[k]->prev=-1;
        verticesvector[k]->distance= (double)distanceMap[vertices[v]];

      }
      int count=0;

      for (typename K::Dart_range::
             iterator it=map.template darts().begin();
           it!=map.template darts().end();++it)
      {
        Dart_handle d1=it;
        count ++;
        int prev1=verticesvector[map.template info<0>(it)]->prev;

        int prev2=verticesvector[map.template info<0>(map.template alpha<0>(it))]->prev;
        int ancestors=0;
        if((prev1==map.template info<0>(map.template alpha<0>(it))||
            prev2==map.template info<0>(it))
           ){

          map.template mark(it,amark2);
          map.template mark(map.template alpha<2>(it),amark2);
          map.template mark(map.template alpha<0>(it),amark2);
          map.template mark(map.template alpha<2>(map.template alpha<0>(it)),amark2);
        }

      }

    }

    void CGAL_breadth_first_search(Dart_handle& source,typename K::size_type& amark2,Graph& g,
                                   std::vector<Vertex>& vertices
                                   ) {
      typedef typename K::Vertex_attribute vertex;
      std::vector<Vertex> predecessors(boost::num_vertices(g)); // To store parents

      std::vector<Weight> distances(boost::num_vertices(g)); // To store distances

      IndexMap indexMap = boost::get(boost::vertex_index, g);

      PredecessorMap predecessorMap(&predecessors[0], indexMap);
      DistanceMap distanceMap(&distances[0], indexMap);
      boost::breadth_first_search(g, vertices[map.template info<0>(source)]
                                  , boost::distance_map(distanceMap).predecessor_map(predecessorMap));

      NameMap nameMap = boost::get(boost::vertex_name, g);
      int k=0;

      for (int k=0;k<vertices.size();k++) {

        Vertex v=vertices[k];
        if (k!=map.template info<0>(source))
          verticesvector[k]->prev=atoi(nameMap[predecessorMap[vertices[v]]].c_str());
        else
          verticesvector[k]->prev=-1;
        verticesvector[k]->distance= (double)distanceMap[vertices[v]];

      }
      int count=0;

      for (typename K::Dart_range::
             iterator it=map.template darts().begin();
           it!=map.template darts().end();++it)
      {
        Dart_handle d1=it;
        count ++;
        int prev1=verticesvector[map.template info<0>(it)]->prev;

        int prev2=verticesvector[map.template info<0>(map.template alpha<0>(it))]->prev;
        int ancestors=0;
        if((prev1==map.template info<0>(map.template alpha<0>(it))||
            prev2==map.template info<0>(it))
           ){

          map.template mark(it,amark2);
          map.template mark(map.template alpha<2>(it),amark2);
          map.template mark(map.template alpha<0>(it),amark2);
          map.template mark(map.template alpha<2>(map.template alpha<0>(it)),amark2);
        }

      }

    }



    /**
     *
     * A fast implementation of Dijkstra algorihtm
     *
     * source : source of the paths
     * amark2 : the darts along dijkstra edges are marked
     * mapS2 : internal face graph structure
     *
     */

    void dijkstra_path3FV(typename K::Dart_handle& source,typename K ::size_type amark2
                          ,bool geometric) {
      int indexsommets=0;
      std::list<int> Q;


      std::vector<int> indexed;
      for (int k=0;k<verticesvector.size();k++) {
        Attrib_h0 it=verticesvector[k];
        it->distance=1000000;
        it->used=0;
        it->prev=-1;
        it->exist=0;

      }
      int j=0;


      int dart=0;

      int i=0;

      int debut=map.template info<0>(source);

      Q.push_back(map.template info<0>(verticesvector[debut]->dart()));
      verticesvector[debut]->distance=0;
      verticesvector[debut]->exist=0;


      int distinfinity=1000000;
      while (Q.size()>0) {
        double dist=-1;int Qr=-1;
        double distancemin=1000000;Attrib_h0 u;
        int Qr2=-1;
        for (std::list<int>::iterator it=Q.begin();it!=Q.end();it++) {
          if (verticesvector[*it]->distance<distancemin) {
            distancemin=verticesvector[*it]->distance;

            u=verticesvector[*it];

            Qr2=*it;
          }

        }
        for (std::list<int>::iterator it=Q.begin();it!=Q.end();it++) {
          if (Qr2==*it) {
            verticesvector[*it]->used=1;
            Q.erase(it);

            break;
          }
        }

        typename K::Dart_handle dh0=verticesvector[Qr2]->template dart();
        //boost::unordered_map<int,int> neigh=find_neighbors_vertex(dh0);
        std::vector<int> neigh=find_neighbors_vertex(dh0);
        for (std::vector<int>::iterator it2=neigh.begin();it2!=neigh.end();it2++){
          if (verticesvector[*it2]->exist==0&&verticesvector[*it2]->used==0) {
            Q.push_back(map.template info<0>(verticesvector[*it2]->dart()));
            verticesvector[*it2]->exist=1;
          }
        }

        for (std::vector<int>::iterator it2=neigh.begin();it2!=neigh.end();it2++){
          double alt=0;
          if (geometric==false)
            alt=distancemin+1;

          if (verticesvector[*it2]->exist==0&&verticesvector[*it2]->used==0) {
            Q.push_back(map.template info<0>(verticesvector[*it2]->dart()));
            verticesvector[*it2]->exist=1;
          }

          if (verticesvector[*it2]->used==0) {

            if (verticesvector[*it2]->distance>alt) {
              verticesvector[*it2]->distance=alt;
              int n=*it2;
              verticesvector[*it2]->prev=Qr2;
              for (std::list<int>::iterator it=Q.begin();it!=Q.end();it++) {
                if ((int)(*it2)==(int)(*it)) {
                  Q.erase(it);
                  break;
                }
              }
              Q.push_front(n);

            }
          }

        }
      }


      for (int k=0;k<verticesvector.size();k++) {
        Attrib_h0 it=verticesvector[k];

        std::cout << it->distance <<" distance\n";
        std::cout << it->used <<" used\n";
        std::cout << it->prev <<" prev\n";
      }

      for (typename K ::Dart_range::iterator it=map.template darts().begin();it!=map.template darts().end();it++) {

        Dart_handle d1=it;

        int prev1=verticesvector[map.template info<0>(it)]->prev;

        int prev2=verticesvector[map.template info<0>(map.template alpha<0>(it))]->prev;


        int ancestors=0;
        if(prev1==map.template info<0>(map.template alpha<0>(it))||
           prev2==map.template info<0>(it)
           ){
          std::cout << "dart marked\n";
          map.template mark(it,amark2);
          map.template mark(map.template alpha<2>(it),amark2);
          map.template mark(map.template alpha<0>(it),amark2);
          map.template mark(map.template alpha<2>(map.template alpha<0>(it)),amark2);
        }

      }
    }


    std::vector<int> bread_first(int i0,int i1,int i2,std::vector<Attrib_h0>& tree) {

      /*for (int k=0;k<verticesvector.size();k++) {
        verticesvector[k]->prev=-1;
        verticesvector[k]->index0=-1;
        }*/
      int cc0=map.template count_all_cells()[0];
      std::map<int,int> l1;
      std::vector<int> leafs;
      l1.insert(std::make_pair(0,0));
      tree.push_back(verticesvector[0]);
      std::cout << map.template info<0>(verticesvector[0]->dart())<< "toto\n";
      std::map<int,int> done;
      done.insert(std::make_pair(0,0));
      int Dist=0;
      while (true) {
        bool neighbors=false;
        std::map<int,int> l2;

        for (std::map<int,int>::iterator it=l1.begin();it!=l1.end();it++) {
          Attrib_h0 a=verticesvector[it->second];
          std::cout << it->second << " toto\n";
          done.insert(std::make_pair(it->second,it->second));
          int dontdo=0;
          Dart_handle dh=verticesvector[it->second]->dart();
          std::vector<int> neigh=this->find_neighbors_vertex(dh);
          for (int k=0;k<neigh.size();k++) {

            Attrib_h0 a2=verticesvector[neigh[k]];

            if ((done.find(neigh[k])==done.end()||done.find(it->second)==done.end())) {
              l2.insert(std::make_pair(neigh[k],neigh[k]));
              a2->prev=it->second;
              a2->distance=Dist;
              tree.push_back(a2);
              a2->index0=tree.size()-1;
              done.insert(std::make_pair(neigh[k],neigh[k]));
              //done.insert(std::make_pair(map.template info<0>(d),map.template info<0>(d)));
            } else {
              leafs.push_back(neigh[k]);
            }
          }
        }
        Dist++;
        l1=l2;
        if (l1.size()==0)
          break;
      }

      std::map<int,int> leafsm;
      for (int k=0;k<leafs.size();k++) {
        leafsm.insert(std::make_pair(leafs[k],leafs[k]));
      }
      leafs.clear();
      for (std::map<int,int>::iterator it=leafsm.begin();it!=leafsm.end();it++) {
        leafs.push_back(it->second);
      }

      for (int k=0;k<leafs.size();k++) {

        int start=leafs[k];
        while (true) {
          start=verticesvector[start]->prev;
          if (start==-1)
            break;
        }

      }
      std::cout << tree.size()<< " no cycle in tree\n";

      return leafs;
    }

    typedef Edge_and_next_edge<K> Edge_and_next_edge0;
    std::vector<Edge_and_next_edge0> make_darts_tree(std::vector<Attrib_h0>& tree
                                                     ,std::vector<Edge_and_next_edge0> leafs0,std::vector<int>& leafs
                                                     ) {
      mark0=map.template get_new_mark();
      std::vector<Edge_and_next_edge0> darts_tree;
      //std::cout << "toto\n";
      for (int k=0;k<tree.size();k++) {
        int D=map.template info<0>(tree[k]->dart());
        int D2=tree[k]->prev;
        //std::cout << D << " " <<(tree[k]->distance)<< " "<<
        //		D2 <<"\n";
        if (tree[k]->prev!=-1) {
          Attrib_h0 a1=tree[k];
          for (typename K::template Dart_of_orbit_range<1,2>::iterator it2=map.template darts_of_orbit<1,2>(a1->dart()).begin();
               it2!=map.template darts_of_orbit<1,2>(a1->dart()).end();it2++) {

            if (map.template info<0>(map.template alpha<0>(it2))==D2) {
              Dart_handle d1=it2;//debut arete
              Dart_handle d2=map.template alpha<0>(verticesvector[D2]->dart());//fin arete
              Edge_and_next_edge0 e;
              e.de1=d1;
              e.de2=d2;
              map.template mark(it2,mark0);
              map.template mark(d2,mark0);
              darts_tree.push_back(e);
              std::cout << "marked\n";
              break;
            }
          }
        } else {
          Dart_handle root=tree[k]->dart();
          Edge_and_next_edge0 e;
          e.de1=root;
          e.de2=root;
          map.template mark(root,mark0);
          //map.template mark(e.de2,marker);
          darts_tree.push_back(e);
          e.degree=-1;
        }
      }
      //std::cout << "toto\n";
      for (int k2=0;k2<darts_tree.size();k2++) {
        for (int k=k2+1;k<darts_tree.size();k++) {
          if (k!=k2) {
            if(map.template info<0>(darts_tree[k].de2)==
               map.template info<0>(darts_tree[k2].de1)) {
              darts_tree[k].dne1=&(darts_tree[k2]);
            }
          }
        }
        //std::cout << (darts_tree.size()-k2) << " \n";
      }
      /*for (int k=0;k<darts_tree.size();k++) {
        if darts_tree[k].degree!=-1
        darts_tree[k].degree=0;
        //std::cout << map.template info<0>((darts_tree[k].dne1)->de1) << "***\n";
        }*/
      for (int k=0;k<darts_tree.size();k++) {
        if(darts_tree[k].dne1&&(darts_tree[k].dne1)->degree!=-1)
          (darts_tree[k].dne1)->degree=(darts_tree[k].dne1)->degree+1;
      }

      //std::cout << "toto\n";
      typename std::vector<Edge_and_next_edge0>::iterator it=darts_tree.begin();
      //for (int k=0;k<leafs.size();k++) {
      for (int k2=0;k2<darts_tree.size();k2++) {

        //if (map.template info<0>((darts_tree[k2].de1))==leafs[k]) {
        if (darts_tree[k2].degree==0) {
          leafs0.push_back(*it);
          std::cout << leafs0.size() << " " << darts_tree.size() << " degree\n";
        }
        it++;
      }
      //}
      std::cout << "end_darts_tree\n";
      edgesnumber=darts_tree.size();
      return leafs0;
    }

    void contract_tree_primal2() {



      int K0=0;
      int size0=248;//map.template count_all_cells()[0]-1;

      for (typename K::Dart_range::iterator it=
             map.template darts().begin();
           it!=map.template darts().end();it++){
        //if (map.template count_all_cells()[1]==0)
        //	break;
        //std::cout << "try contraction1\n";
        if (map.template is_marked(it,mark0)) {
          //std::cout << "try contraction2\n";
          if (map.template is_contractible<1>(it)) {
            std::cout << K0 <<" is_contractible0\n";
            map.template contract_cell<1>(it,true);
            /*Dart_handle it2=map.template alpha<0>(it);
              Edge_and_next_edge<K> E1;
              E1.de1=it;
              E1.de2=it2;
              contract_edge(E1,mark0);*/
            K0++;

          }
        }

      }

      std::cout <<map.template count_all_cells()[0] <<"finished\n";
      map.template display_characteristics(std::cout);

    }

    void set_attributes(bool erase_attributes) {

      typedef typename K ::size_type size_type;
      int i=0;
      int i2=0;
      int i1=0;
      typedef typename K::Vertex_attribute vertex;
      typedef typename K::template Attribute_handle<0>::type Attrib_h;
      typedef typename K::template Attribute_handle<2>::type Attrib_h2;

      if (erase_attributes) {

        for (typename K::Vertex_attribute_range::iterator it=map.vertex_attributes().begin();
             it!=map.vertex_attributes().end();it++){
          map.template erase_vertex_attribute(it);
        }
        for (typename K::template Attribute_range<1>::type::iterator it=map.template attributes<1>().begin();
             it!=map.template attributes<1>().end();it++){
          map.template erase_attribute<1>(it);
        }
        for (typename K::template Attribute_range<2>::type::iterator it=map.template attributes<2>().begin();
             it!=map.template attributes<2>().end();it++){
          map.template erase_attribute<2>(it);
        }
      }

      std::cout << "tototototo\n";

      verticesvector.clear();
      facesvector.clear();
      int first=1;

      for (typename K::template One_dart_per_cell_range<0>::
             iterator it=map.template one_dart_per_cell<0>().begin();
           it!=map.template one_dart_per_cell<0>().end();++it)
      {

        int index=0;

        first=0;

        typename K::template Attribute_handle<0>::type ah;
        ah=map.template create_vertex_attribute(map.template point_of_vertex_attribute(map.template vertex_attribute(it)));//map.template vertex_attribute(it);
        map.template set_vertex_attribute(it,ah );
        map.template info<0>(it)=i;

        verticesvector.push_back(ah);
        ah->template prev=-1;
        ah->template distance=0;
        i++;

      }
      std::cout << "tototototo\n";
      first=1;
      if (false) {
        for (typename K::template One_dart_per_cell_range<1>::
               iterator it=map.template one_dart_per_cell<1>().begin();
             it!=map.template one_dart_per_cell<1>().end();++it)
        {
          typename K::template Attribute_handle<1>::type ah;
          ah=map.template create_attribute<1>();//map.template vertex_attribute(it);
          map.template set_attribute<1>(it,ah );
          map.template info<1>(it)=i1;
          i1++;
        }
      }
      std::cout << "tototototo\n";
      first=1;
      for (typename K:: template One_dart_per_cell_range<2>::
             iterator it=map.template one_dart_per_cell<2>().begin();
           it!=map.template one_dart_per_cell<2>().end();++it)
      {

        first=0;
        typename K::template Attribute_handle<2>::type ah=//map.template attribute<2>(it);
          map.template create_vertex_attribute(map.template point_of_vertex_attribute(map.template vertex_attribute(it)));//map.template vertex_attribute(it);

        //std::cout << (map.template darts_of_cell<2>(it).size()) << " facenumber\n";
        map.template set_attribute<2>(it,ah);
        map.template info<2>(it)=i2;
        ah->boundary=0;
        facesvector.push_back(ah);
        i2++;

      }


      for (std::map<int,int>::iterator it=boundary.begin();it!=boundary.end();it++) {

        facesvector[it->second]->boundary=1;

      }

      //map.template onsplit_functor<1>()=CGAL::Split_functor<K>(map);
      //map.template onmerge_functor<1>()=CGAL::Merge_functor();


      std::cout << "tototototo\n";


    }

    void make_A_Map() {

      Dart_handle d1=map.template make_combinatorial_polygon(4);
      Dart_handle d0=d1;
      Dart_handle d2=map.template alpha<1>(d0);
      d2=map.template alpha<0>(d2);
      d2=map.template alpha<1>(d2);
      //map.template link_alpha<2>(d0,d2);
      Dart_handle d3=map.template alpha<0>(d0);
      Dart_handle d4=map.template alpha<0>(d2);
      //map.template link_alpha<2>(d3,d4);
      Dart_handle d12=map.template alpha<1>(d1);
      d12=map.template alpha<0>(d12);
      Dart_handle d32=map.template alpha<1>(d3);
      Dart_handle d33=map.template alpha<0>(d32);

      map.template link_alpha<2>(d2,d0,true);
      map.template link_alpha<2>(d3,d4,true);
      map.template link_alpha<2>(d1,d32,true);
      map.template link_alpha<2>(d12,d33,true);
      //Dart_handle d2=_map.template make_combinatorial_polygon(7);
      int i=0;
      /*for (typename K::template One_dart_per_cell_range<0>::
        iterator it=map.template one_dart_per_cell<0>().begin();
        it!=map.template one_dart_per_cell<0>().end();++it)
        {

        int index=0;

        typename K::template Attribute_handle<0>::type ah;
        ah=map.template create_vertex_attribute(Point(0,0,0));//map.template vertex_attribute(it);
        map.template set_vertex_attribute(it,ah );
        map.template info<0>(it)=i;
        std::cout << "1\n";
        verticesvector.push_back(ah);
        ah->template prev=-1;
        ah->template distance=0;
        i++;

        }

        int i2=0;
        for (typename K:: template One_dart_per_cell_range<2>::
        iterator it=map.template one_dart_per_cell<2>().begin();
        it!=map.template one_dart_per_cell<2>().end();++it)
        {

        typename K::template Attribute_handle<2>::type ah=//map.template attribute<2>(it);
        map.template create_vertex_attribute(Point(0,0,0));

        map.template set_attribute<2>(it,ah);
        map.template info<2>(it)=i2;
        ah->boundary=0;
        facesvector.push_back(ah);
        i2++;

        }*/
      std::cout << "1\n";
      std::cout << map.template is_valid() << "\n";


    }

    void contract_tree_primal(std::vector<Edge_and_next_edge0> leafs0,int i0,int i1,int i2) {

      typename K::size_type marker;
      int K0=0;
      int size0=248;//map.template count_all_cells()[0]-1;
      while (true) {
        std::cout <<leafs0.size() << " before is_contractible\n";
        int count=0;
        std::map<Edge_and_next_edge0*,Edge_and_next_edge0*> leafs1;
        for (typename std::vector<Edge_and_next_edge0>::iterator it=leafs0.begin();it!=leafs0.end();it++) {

          typename std::vector<Edge_and_next_edge0>::iterator l=it;
          std::cout << "is_contractible0\n";

          //edgesnumber--;
          map.template contract_cell<1>(l->de1,false);
          //contract_edge(*l,marker);
          //map.template contract_cell<1>(l->de1,false);

          if (l->dne1->degree>0)
            l->dne1->degree=l->dne1->degree-1;

          std::cout << count++ << "contracted\n";
          if (l->dne1->degree==0)
            leafs1.insert(std::make_pair(l->dne1,l->dne1));

          K0++;


        }

        std::cout << map.template count_all_cells()[1] << " " << edgesnumber<< " en cours\n";
        if (leafs1.size()==0)
          break;
        std::cout << " en cours2\n";
        //std::deque<Edge_and_next_edge0*> leafs10;
        leafs0.clear();
        for (typename std::map<Edge_and_next_edge0*,Edge_and_next_edge0*>::iterator it0=leafs1.begin();
             it0!=leafs1.end();it0++) {
          leafs0.push_back(*(it0->second));
        }


      }
      std::cout << "finished\n";
      for (typename K::Dart_range::iterator it=map.template darts().begin();
           it!=map.template darts().end();it++) {
        if (map.template is_marked(it,marker))
          map.template mdarts.erase(it);
      }


      set_attributes(true);
      map.template display_characteristics(std::cout);
      //if (leafs.size()==size||leafs.size()==1)
      //	break;
      //verticesvector.clear();
      //facesvector.clear();
    }

    void contract_edge(Edge_and_next_edge0& E,typename K::size_type marker) {
      Dart_handle d1=E.de1;
      Dart_handle d2=E.de2;
      Dart_handle d3=map.template alpha<2>(d1);
      Dart_handle d4=map.template alpha<2>(d2);
      Dart_handle d12=map.template alpha<1>(d1);
      Dart_handle d22=map.template alpha<1>(d2);
      Dart_handle d32=map.template alpha<1>(d3);
      Dart_handle d42=map.template alpha<1>(d4);

      if (!map.template is_free<1>(d1))
        map.template unlink_alpha<1>(d1);
      if (!map.template is_free<1>(d2))
        map.template unlink_alpha<1>(d2);
      if (!map.template is_free<1>(d3))
        map.template unlink_alpha<1>(d3);
      if (!map.template is_free<1>(d4))
        map.template unlink_alpha<1>(d4);
      if (!map.template is_free<2>(d1))
        map.template unlink_alpha<2>(d1);
      if (!map.template is_free<2>(d2))
        map.template unlink_alpha<2>(d2);
      if (!map.template is_free<2>(d3))
        map.template unlink_alpha<2>(d3);
      if (!map.template is_free<2>(d4))
        map.template unlink_alpha<2>(d4);
      /*if (!map.template is_free<0>(d1))
        map.template unlink_alpha<0>(d1);
        if (!map.template is_free<0>(d2))
        map.template unlink_alpha<0>(d2);
        if (!map.template is_free<0>(d3))
        map.template unlink_alpha<0>(d3);
        if (!map.template is_free<0>(d4))
        map.template unlink_alpha<0>(d4);*/
      /*map.template mdarts.erase(d1);
        map.template mdarts.erase(d2);
        map.template mdarts.erase(d3);
        map.template mdarts.erase(d4);*/
      map.template link_alpha<1>(d12,d32);
      map.template link_alpha<1>(d22,d42);

      map.template mark(d1,marker);
      map.template mark(d2,marker);
      map.template mark(d3,marker);
      map.template mark(d4,marker);

    }


  };

}

#endif /* SRC_GM_algorithms_H_ */
