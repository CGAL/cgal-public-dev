#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_traits_d.h>
#include <CGAL/Extreme_points_options_d.h>

typedef CGAL::Cartesian_d<double>               Kernel_d;
typedef Kernel_d::Point_d                       Point_d;
typedef CGAL::Extreme_points_traits_d<Point_d>  EP_Traits_d;
typedef CGAL::Extreme_points_d<EP_Traits_d>     EP_d;

template <class InputStream, class OutputIterator>
void read_input(InputStream &in, OutputIterator out, int &n, int &d) {
	// read *.con file
	std::string s;
	do {
		getline(in,s);
	} while (s[0]=='%'); // skipping comments
	std::stringstream ss(s);
	ss>>d>>n;
	
	std::vector<Point_d> points(n);
	for (int i=0;i<n;++i) {
		std::vector<double> p(d);
		for (int j=0;j<d;++j)
			in>>p[j];
		*out++=Point_d(d,p.begin(),p.end());
	}
}

int n,d,p=1;
std::vector<Point_d> points;

void test(EP_d ep, EP_d ep2) {
    int k=n/20; // n==2500, k==125
    
    const int ROUNDS = 4; // maximal 4 rounds..
    
    // add parts of the point set and compare with static implementation
    for (int i=0;i<ROUNDS;++i) {
        std::vector<Point_d> extreme_points;
        ep.insert(points.begin() + (i*k), points.begin() + ((i+1)*k));
        ep.extreme_points(std::back_inserter(extreme_points));
        
        std::vector<Point_d> extreme_points_ref;
        extreme_points_d_dula_helgason(points.begin(), points.begin() + ((i+1)*k), std::back_inserter(extreme_points_ref));
        
        std::cout<<"First "<<(i+1)*k<<" points:"<<std::endl;
        std::cout<<"Dynamic class found "<<extreme_points.size()<<" extreme points"<<std::endl;
        std::cout<<"Static function found "<<extreme_points_ref.size()<<" extreme points"<<std::endl;
        
        assert(extreme_points.size()==extreme_points_ref.size());
        
        std::sort(extreme_points.begin(),extreme_points.end(), Kernel_d::Less_lexicographically_d());
        std::sort(extreme_points_ref.begin(),extreme_points_ref.end(), Kernel_d::Less_lexicographically_d());
        
        // check that the different implementations produce the same output
        assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points_ref.begin()));
        
        // test classification function..
        std::cout<<"Testing classification function"<<std::endl;
        
        std::set<Point_d, Kernel_d::Less_lexicographically_d > xp(extreme_points.begin(), extreme_points.end());
        
        for (std::vector<Point_d>::iterator it=points.begin();it!=points.begin()+(i+1)*k;it++) {
            if (xp.find(*it)!=xp.end()) { // extreme point
                assert(ep.classify(*it,true)==CGAL::ON_BOUNDARY);
                assert(ep.classify(*it)==CGAL::ON_BOUNDARY);
            } else {
                // internal point
                assert(ep.classify(*it,true)==CGAL::ON_UNBOUNDED_SIDE);
                assert(ep.classify(*it)==CGAL::ON_UNBOUNDED_SIDE);
            }
        }
    }
    
    std::cout<<"Duplicated points tests:"<<std::endl;
    
    // same thing but insert points twice for checking handling of duplicates
    for (int i=0;i<ROUNDS;++i) {
        std::cout<<"First "<<(i+1)*k<<" points:"<<std::endl;
        
        std::vector<Point_d> extreme_points;
        // insert everything twice
        ep2.insert(points.begin() + (i*k), points.begin() + ((i+1)*k));
        ep2.insert(points.begin() + (i*k), points.begin() + ((i+1)*k));
        ep2.extreme_points(std::back_inserter(extreme_points));
        
        std::vector<Point_d> extreme_points_ref;
        extreme_points_d_dula_helgason(points.begin(), points.begin() + ((i+1)*k),std::back_inserter(extreme_points_ref));
        
        std::cout<<"Dynamic class found "<<extreme_points.size()<<" extreme points"<<std::endl;
        std::cout<<"Static function found "<<extreme_points_ref.size()<<" extreme points"<<std::endl;
        
        assert(extreme_points.size()==extreme_points_ref.size());
        
        std::sort(extreme_points.begin(),extreme_points.end(), Kernel_d::Less_lexicographically_d());
        std::sort(extreme_points_ref.begin(),extreme_points_ref.end(), Kernel_d::Less_lexicographically_d());
        
        // check that the different implementations produce the same output
        assert(std::equal(extreme_points.begin(),extreme_points.end(),extreme_points_ref.begin()));
        
        // test classification function..
        std::cout<<"Testing classification function"<<std::endl;
        
        std::set<Point_d, Kernel_d::Less_lexicographically_d> xp(extreme_points.begin(), extreme_points.end());
        
        for (std::vector<Point_d>::iterator it=points.begin();it!=points.begin()+(i+1)*k;it++) {
            if (xp.find(*it)!=xp.end()) { // extreme point
                assert(ep2.classify(*it,true)==CGAL::ON_BOUNDARY);
                assert(ep2.classify(*it)==CGAL::ON_BOUNDARY);
            } else {
                // internal point
                assert(ep2.classify(*it,true)==CGAL::ON_UNBOUNDED_SIDE);
                assert(ep2.classify(*it)==CGAL::ON_UNBOUNDED_SIDE);
            }
        }
    }

    std::cout<<"ok"<<std::endl<<std::endl;
}

int main(int argc, char **argv) {
    std::ifstream fin("./data/05by02500at01.con");
    read_input(fin,std::back_inserter(points),n,d);

    std::cout;
    std::cout<<"Testing with default options"<<std::endl;
    std::cout<<"----------------------------"<<std::endl;
    EP_d ep(d);
    EP_d ep2(d);
    test(ep,ep2);
    
    // testing options
    CGAL::Extreme_points_options_d options;
    options.set_anti_cycling(true);    
    assert(CGAL::QP_BLAND==options.get_qp_options().get_pricing_strategy());

    std::cout;
    std::cout<<"Testing with EP_CHOOSE_APPROPRIATE"<<std::endl;
    std::cout<<"----------------------------------"<<std::endl;
    options.set_algorithm(CGAL::EP_CHOOSE_APPROPRIATE);
    ep=EP_d(d,options);
    ep2=EP_d(d,options);
    test(ep,ep2);

    std::cout;
    std::cout<<"Testing with EP_SIMPLE"<<std::endl;
    std::cout<<"--------------------- "<<std::endl;
    options.set_algorithm(CGAL::EP_SIMPLE);
    ep=EP_d(d,options);
    ep2=EP_d(d,options);
    test(ep,ep2);

    std::cout;
    std::cout<<"Testing with EP_DULA_HELGASON"<<std::endl;
    std::cout<<"-----------------------------"<<std::endl;
    options.set_algorithm(CGAL::EP_DULA_HELGASON);
    ep=EP_d(d,options);
    ep2=EP_d(d,options);
    test(ep,ep2);

    ep.clear();
    ep2.clear();    

    std::cout<<"Finished successfully!"<<std::endl;
    return 0;
}
