/*
 * Test_algo1.cpp
 *
 *  Created on: 29 août 2017
 *      Author: riedingc
 */
//#include "Generalized_map_algorithmtest2.h"
#include <CGAL/GM_algorithms.h>
//#include "Generate_weights_filetest.cpp"
#include <CGAL/Generalized_map.h>
#include <CGAL/Cell_attribute.h>
#include <boost/unordered_map.hpp>
#include <algorithm>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Timer.h>

//typedef CGAL::Linear_cell_complex_for_generalized_map<2,3,CGAL::Traits,CGAL::MyitemPt> LCC_3;

class Test_algo1test {

public:

	void test(char* name,int geometric,int boundary) {

		//CGAL::Generalized_map<2,CGAL::My_item> map1;

		//CGAL::Generalized_map_algorithms4 gma(map2);
		CGAL::LCC_3 lcc;
		load_and_simplify_off(lcc, name,
				true, 100);
		std::cout << "toto\n";
		//CGAL::Generate_weights_filetest<CGAL::LCC_3> gwft;
		std::vector<int> W;
		/*if (geometric)
			W=gwft.read_weights_file("weightsft");*/

		CGAL::GM_algorithms<CGAL::LCC_3 > gma(lcc);
		std::map<int,int> boundaries;
		std::vector<unsigned int> size=lcc.count_all_cells();
		if (boundary!=0) {
			for (int k=0;k<boundary;k++)
				boundaries.insert(std::make_pair(size[2]-k-1,size[2]-k-1));
		}
		gma.boundary=boundaries;

		std::cout << "toto \n";

		//std::map<int,std::map<int,double> > W;
				///gma.display(gma.map.darts().begin());
				gma.compute_shortest_cycle_non_contractible(geometric,boundaries,W,1,1);



				CGAL::LCC_3 lcc0;
				load_and_simplify_off(lcc0, name,
								true, 100);


				CGAL::GM_algorithms<CGAL::LCC_3 > gma0(lcc0);
				gma0.boundary=boundaries;
				//gma0.set_attributes(false);
				//typedef typename CGAL::LCC_3::Attribute_handle<0>::type Attrib_h0;
				//std::vector<Attrib_h0> tree;
				//std::vector<int> leafs=gma0.bread_first(0,1,2,tree);
				//std::cout << "tutu\n";
				//typedef CGAL::Edge_and_next_edge<CGAL::LCC_3> Edge_and_next_edge0;
				//std::vector<Edge_and_next_edge0> leafs0;

				//leafs0=gma0.make_darts_tree(tree,leafs0,leafs);
				//std::cout << "toto\n";

				typedef double Weight;
					 typedef boost::property<boost::edge_weight_t, Weight> WeightProperty;
					 typedef boost::property<boost::vertex_name_t, std::string> NameProperty;
				typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS,
					 	    NameProperty, WeightProperty > Graph;
				typedef boost::graph_traits < Graph >::vertex_descriptor Vertex;
				Graph g;
				std::vector<int > W2;
				gma0.set_attributes(false);
				std::vector<Vertex> vertices = gma0.create_CGAL_Graph2(g,0,W2);
				CGAL::LCC_3::Dart_handle dh=gma0.map.darts().begin();
				gma0.mark0=gma0.map.get_new_mark();
				//gma0.CGAL_breadth_first_search(dh,gma0.mark0,g,vertices);
				gma0.CGAL_Dijkstra(dh,gma0.mark0,g,vertices);
				gma0.contract_tree_primal2();

				gma0.map.display_characteristics(std::cout);
				gma0.set_attributes(false);
				gma0.map.unmark_all(gma0.mark0);
				std::vector<int> Wdual;
				Graph g2;
				std::vector<Vertex> verticesdual =gma0.create_CGAL_dual_Graph2(g2,0,Wdual);
				dh=gma0.map.darts().begin();
				gma0.CGAL_dual_Dijkstra(dh,gma0.mark0,g2,verticesdual);
				gma0.remove_edges();
				gma0.map.display_characteristics(std::cout);
				///gma0.compute_shortest_cycle_non_contractible(geometric,boundaries,W2,1,0);
				gma0.map. display_characteristics(std::cout);
				std::cout << "\n";
				//CGAL::write_off(gma.map,"torecontracted.off");

				CGAL::Generalized_map<2,CGAL::Generic_map_min_items> lcc2;
				CGAL::GM_algorithms<CGAL::Generalized_map<2,CGAL::Generic_map_min_items> > gma2(lcc2);
				gma2.make_A_Map();

				std::cout << gma2.map.is_valid() << "\n";
				gma2.boundary=boundaries;
				std::cout << "AMAP\n";
				std::vector<int > W3;
				///gma2.compute_shortest_cycle_non_contractible(geometric,boundaries,W3,0,0);
	}



		void load_and_simplify_off(CGAL::LCC_3& lcc, const std::string& filename,
		bool updateattribs, int percent)
		{
		std::ifstream ifile(filename.c_str());
		if (ifile)
		{
		CGAL::load_off<CGAL::LCC_3>(lcc, ifile);

		lcc.display_characteristics(std::cout);


		CGAL::Timer timer;
		CGAL::LCC_3::Dart_handle dh;
		std::size_t nb=lcc.number_of_darts();//(*percent)/200;
		timer.start();
		//if (!updateattribs) lcc.set_automatic_attributes_management(false);
		/*for (CGAL::LCC_3::Dart_range::iterator it=lcc.darts().begin(),
		itend=lcc.darts().end(); it!=itend && nb>0; )
		{
		dh=it++;
		if ( it!=itend && it==lcc.alpha<2>(dh) ) ++it;
		lcc.remove_cell<1>(dh);
		--nb;
		}*/
		//if ( updateattribs )
			//lcc.set_automatic_attributes_management(true);
		timer.stop();
		lcc.display_characteristics(std::cout);
		std::cout<<", valid="<< lcc.is_valid()
		<<" time: "<<timer.time()<<" seconds." << std::endl;



		}
		}


};

