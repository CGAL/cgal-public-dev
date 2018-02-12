#ifndef ISLAND_TRIANGULATE_HOLE_POLYLINE_H
#define ISLAND_TRIANGULATE_HOLE_POLYLINE_H

#include <vector>
#include <tuple>
#include <stack>
#include <CGAL/Combination_enumerator.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>


namespace CGAL {
namespace internal {


// Domain structure //
// ---------------- //

template <typename PointRange>
struct Domain
{

  Domain() {}

  // boundary should always be given without the stupid last(first) one.
  // not used
  Domain(PointRange& boundary) : boundary(boundary)
  {
    std::size_t n_p = boundary.size();
    if(boundary.size() > 0)
      assert(boundary[n_p - 1] != boundary[0]);
  }

  // constructor with indices
  Domain(std::vector<int> ids) : b_ids(ids) {}

  ~Domain() {}

  void add_hole(std::vector<int>& ids)
  {
    holes_list.push_back(ids);

    for(int i=0; i<ids.size(); ++i)
      h_ids.push_back(ids[i]);
  }

  void clear_holes()
  {
    holes_list.clear();
    h_ids.clear();
  }

  bool is_empty()
  {
    return b_ids.size() == 2 ? true : false;
  }

  bool has_islands()
  {
    return holes_list.empty() ? false : true;
  }

  std::pair<int, int> get_access_edge()
  {
    std::size_t number_of_points = b_ids.size();
    assert(number_of_points >= 2);

    int i = b_ids[0];
    int k = b_ids[number_of_points - 1];

    //assert(i != k);

    return std::make_pair(i, k);
  }

  void print_boundary()
  {
    print(boundary);
  }

  void print_size()
  {
    std::cout << boundary.size() << std::endl;
  }

  //data
  PointRange boundary; // not used in the main algorithm
  std::vector<PointRange> holes; // not used

  std::vector<int> b_ids;
  std::vector<int> h_ids;
  std::vector<std::vector<int>> holes_list;

};


template <typename T>
void print(T &v)
{
  for(int i=0; i<v.size(); ++i)
  {
    std::cout << v[i] << " ";//<< std::endl;
  }
  std::cout << std::endl;
}

template <typename T, typename PointRange>
void print(T &v, std::ofstream& out, int& i, PointRange points)
{
  if(v.size() <= 2)
    return;

  std::string filename("data/seq"+std::to_string(i)+".polylines.txt");
  out.open(filename); //std::ofstream::app

  out << v.size() + 1 << " ";
  for(int i=0; i<v.size(); ++i)
  {
    out << points[v[i]] << " ";//<< std::endl;
  }
  out << points[v[0]] << std::endl;
  out << std::endl;

  out.close();
  i++;

}



template <typename T>
void print_append(T &v, std::ofstream& out)
{
  if(v.size() <= 2)
    return;

  std::string filename("data/domains.dat");
  out.open(filename , std::ofstream::app); //

  for(int i=0; i<v.size(); ++i)
  {
    out << v[i] << " ";
  }
  out << std::endl;

  out.close();

}


template <typename PointRange>
void print_triangle(int i, int m, int k, std::ofstream& out, int& ii, PointRange points)
{

  std::string filename("data/tr"+std::to_string(ii)+".polylines.txt");
  out.open(filename); //std::ofstream::app

  out << 4 << " ";
  out << points[i] << " " << points[m] << " " << points[k] << " " << points[i];
  out << std::endl;

  out.close();
  ii++;

}


// partition permutations // 
// ---------------------- //

struct Phi
{
  typedef std::pair<std::vector<int>, std::vector<int>> Sub_domains_pair;
  void put(std::vector<int>& left, std::vector<int>& right)
  {
    sub_domains_list.push_back(std::make_pair(left, right)); // preallocate list
  }

  std::size_t size()
  {
    return sub_domains_list.size();
  }

  bool empty()
  {
    return size() == 0 ? true : false;
  }

  std::vector<int> lsubset(const int i)
  {
    assert(i >= 0);
    assert(i < sub_domains_list.size());
    return sub_domains_list[i].first;
  }

  std::vector<int> rsubset(const int i)
  {
    assert(i >= 0);
    assert(i < sub_domains_list.size());
    return sub_domains_list[i].second;
  }

  std::vector<Sub_domains_pair> sub_domains_list;
};


void do_permutations(std::vector<std::vector<int>>& hole_list, Phi& subsets)
{
  if(hole_list.empty())
    return;

  std::vector<int> hs;

  for(int n = 0; n < hole_list.size(); ++n)
    hs.push_back(n);

  const int first = hs.front();
  const int last = hs.back();
  //std::sort(hs.begin(), hs.end()); // already sorted

  for(int s = 0; s <= hs.size(); ++s) // s = number of holes on one (left) side
  {
    std::vector<int> p1(s);
    std::vector<int> p2(hole_list.size() - s);

    if(s == 0)
    {
      subsets.put(p1, hs);

      print(p1); std::cout << "-- "; print(hs); std::cout << std::endl;
      continue;
    }

    CGAL::Combination_enumerator<int> permutations(s, first, last+1);

    int p = 0;
    while(!permutations.finished())
    {
      for(int i=0; i<s; ++i)
      {
        p1[i] = permutations[i];
      }

      ++permutations;

      std::sort(p1.begin(), p1.end());
      std::set_symmetric_difference(p1.begin(), p1.end(), hs.begin(), hs.end(),
                                    p2.begin());

      assert(p1.size() == s);
      assert(p2.size() == hs.size() - s);

      print(p1); std::cout << "-- "; print(p2); std::cout << std::endl;

      subsets.put(p1, p2);
    }

  }
}

// split //
// ----- //

template <typename PointRange>
void split_domain(const Domain<PointRange>& init_domain,
                    Domain<PointRange>& left_dom, Domain<PointRange>& right_dom,
                    const int& i, const int& pid, const int& k)
{

  typedef std::vector<int> Ids;
  Ids ids = init_domain.b_ids;
  Ids left;
  Ids right;

  // i, k indices of access edge = first and last

  // find position of pid
  Ids::iterator it;
  it = std::find(ids.begin(), ids.end(), pid);
  assert(it != ids.end());

  left.insert(left.begin(), ids.begin(), it + 1);
  right.insert(right.begin(), it, ids.end());


  assert(left.front() == i);
  assert(left.back() == pid);
  assert(right.front() == pid);
  assert(right.back() == k);

  left_dom.b_ids = left; // todo: avoid copying
  right_dom.b_ids = right;

}


void reorder_island(std::vector<int>& h_ids, const int& v)
{

  // 1) find v's position
  std::vector<int>::iterator it = find(h_ids.begin(), h_ids.end(), v);
  assert(it != h_ids.end());

  // 2) rotate by the third vertex of t
  std::size_t dist = std::distance(h_ids.begin(), it); // std::size_t?
  std::rotate(h_ids.begin(), h_ids.begin() + dist, h_ids.end());

  // 3) add the first removed element
  h_ids.push_back(h_ids[0]);

  // 4) reverse order
  std::reverse(h_ids.begin(), h_ids.end());

}


void do_not_reorder_island(std::vector<int>& h_ids, const int& v)
{

  // 1) find v's position
  std::vector<int>::iterator it = find(h_ids.begin(), h_ids.end(), v);
  assert(it != h_ids.end());

  // 2) rotate by the third vertex of t
  std::size_t dist = std::distance(h_ids.begin(), it); // std::size_t?
  std::rotate(h_ids.begin(), h_ids.begin() + dist, h_ids.end());

  // 3) add the first removed element
  h_ids.push_back(h_ids[0]);

}


void merge_id_sets(std::vector<int>& b_ids,
                   const int& i, const int& v, const int& k,
                   std::vector<int>& hole_ids, bool reorder)
{
  std::size_t initial_b_size = b_ids.size();

  if(reorder)
    reorder_island(hole_ids, v);
  else
    do_not_reorder_island(hole_ids, v);

  // insertion position = just after k
  // k is at position n - 1 = last element.
  // just append at the end - i is the first point on b_ids
  // and k is the last. t triangle is (i, v, k)
  typename std::vector<int>::iterator insertion_point = b_ids.end();
  b_ids.insert(insertion_point, hole_ids.begin(), hole_ids.end());

  //assert(*(b_ids.begin() + i) == b_ids[i] );

  assert(b_ids[initial_b_size - 1] == k);
  assert(b_ids[0] == i);
  assert(b_ids[initial_b_size] == v);
  assert(b_ids[b_ids.size() - 1] == v);
  assert(b_ids.size() == initial_b_size + hole_ids.size());

}

template<typename PointRange>
void join_domain(const Domain<PointRange>& domain, Domain<PointRange>& D1, Domain<PointRange>& D2,
                 const int& i, const int& v, const int& k)
{
  typedef std::vector<int> Ids;
  Ids id_set1 = domain.b_ids;
  Ids id_set2 = domain.b_ids;
  Ids hole_ids1 = domain.h_ids; // for now assume just one hole.
  Ids hole_ids2 = domain.h_ids; // for now assume just one hole.

  // merge once without reordering
  merge_id_sets(id_set1, i, v, k, hole_ids1, false);
  D1.b_ids = id_set1;

  // merge again with reordering island
  merge_id_sets(id_set2, i, v, k, hole_ids2, true);
  D2.b_ids = id_set2;
}


struct Tracer
{

  template <class LookupTable>
  void operator()(const LookupTable& lambda, int v0, int v1)
  {
    CGAL_assertion_code( const int n = lambda.n; )
    std::stack<std::pair<int, int> > ranges;
    ranges.push(std::make_pair(v0, v1));

    while(!ranges.empty()) {
      std::pair<std::size_t, std::size_t> r = ranges.top();
      ranges.pop();
      CGAL_assertion(r.first >= 0 && r.first < n);
      CGAL_assertion(r.second >= 0 && r.second < n);

      // if on border
      if(r.first + 1 == r.second) {
        continue; }
      if(r.first + 1 == n) {
        continue;
      }

      std::size_t la = lambda.get(r.first, r.second);
      if(la == -1) {
          std::cerr << "out hole" << std::endl;
          //*out_hole++ = std::make_pair(r.first, r.second);
        continue;
      }

      CGAL_assertion(la >= 0 && la < n);
     // CGAL_assertion(r.first < la && r.second > la); not with islands
      //auto triangle = std::make_tuple(r.first, la, r.second);

      std::vector<std::size_t> triangle = {r.first, la, r.second};
      collection.push_back(triangle);

      ranges.push(std::make_pair(r.first, la));
      ranges.push(std::make_pair(la, r.second));
    }
  }


  //std::vector<std::tuple<int, int, int>> collection;
  std::vector<std::vector<std::size_t>> collection;
};






template<typename PointRange, typename WeightCalculator,
         typename WeightTable, typename LambdaTable>
class Triangulate
{
  typedef typename WeightCalculator::Weight Weight;


public:

  Triangulate(Domain<PointRange> domain,
              PointRange allpoints,
              WeightTable& W,
              LambdaTable& l,
              const WeightCalculator & WC) :
              points(allpoints),
              W(W),
              lambda(l),
              domain(domain),
              WC(WC){}

  std::size_t do_triangulation(int& i, int& k, std::size_t& count)
  {

    print_i = 1;

    init_triangulation();

    processDomain(domain, i, k, count);

    //ambda.print("data/lambda-rec.dat");
    //W.print("data/weight-rec.dat");


  }

  void collect_triangles(std::vector<std::vector<std::size_t>>& triplets,
                         int& i, int& k)
  {
    Tracer tracer;
    tracer(lambda, i, k);
    triplets = tracer.collection;

    assert(triplets.size() > 0);
  }

  template <typename PolygonMesh>
  void visualize(PointRange& points, std::vector<std::vector<std::size_t>>& polygon_soup,
                 PolygonMesh& mesh)
  {
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygon_soup, mesh);
  }




private:


  void init_triangulation()
  {
    for(int i=0; i < domain.b_ids.size(); ++i)
    {
      init_b.push_back(domain.b_ids[i]);
    }

    for(int i=0; i < domain.h_ids.size(); ++i)
    {
      init_island.push_back(domain.h_ids[i]);
    }

  }


  // main loop //
  // --------- //
  void processDomain(Domain<PointRange> domain, int& i, int& k, std::size_t& count)
  {
    // (i, k) = acccess edge

    if (i == k)
    {
      std::cout << "on domain: ";
      for(int j=0; j<domain.b_ids.size(); ++j)
      {
        std::cout << domain.b_ids[j] << " ";
      }
      std::cout << std::endl;
      std::cout <<"i == k: " << i << "=" << k << " returning..." <<std::endl;
      return;
    }


    //std::cout << "count: " << count << std::endl;
    //print_append(domain.b_ids, out_domain);

    if(!domain.b_ids.size() == 2)
    {
      std::cout << "domain = ";
      for(int j=0; j<domain.b_ids.size(); ++j)
      {
        std::cout << domain.b_ids[j] << " ";
      }
      std::cout << std::endl;
    }


    // empty domain
    if(domain.b_ids.size() == 2)
    return;


    // base case
    if(domain.b_ids.size() == 3 && domain.holes_list.empty())
    {

      assert(domain.b_ids[0] == i); // access edge source
      assert(domain.b_ids[2] == k); // access edge target


      int m = domain.b_ids[1]; //third vertex

      std::cout<<"BASE CASE - evaluating triangle = ("<<i<<","<<m<<","<<k<<")"<<std::endl;
      calculate_weight(i, m, k);

      count++;

      //std::cin.get();

      return;
    }
    assert(domain.b_ids.size() >= 3);




    // CASE I - if there are islands, join until there are no islands.
    for(int pid : domain.h_ids)
    {
      std::cout << "------------- JOIN DOMAIN ------------" << std::endl;

      // avoid source & target of e_D
      if(pid == i || pid == k)
        continue;

      Domain<PointRange> D1;
      Domain<PointRange> D2;
      // D1, D2 correspond to the two different hole orientations
      join_domain(domain, D1, D2, i, pid, k);
      // get a new e_D - todo: const reference
      std::pair<int, int> e_D1 = D1.get_access_edge();
      std::pair<int, int> e_D2 = D2.get_access_edge();


      std::cout << "new domain D1 after join = ";
      for(int j=0; j<D1.b_ids.size(); ++j)
      {
        std::cout << D1.b_ids[j] << " ";
      }
      std::cout << std::endl;

      std::cout << "new domain D2 after join = ";
      for(int j=0; j<D2.b_ids.size(); ++j)
      {
        std::cout << D2.b_ids[j] << " ";
      }
      std::cout << std::endl;

      // second ordering
      processDomain(D2, e_D2.first, e_D2.second, count);
      // after the subdomains left and right have been processed
      assert(domain.has_islands());
      assert(domain.b_ids[0] == i);
      assert(domain.b_ids[domain.b_ids.size() - 1] == k);

      std::cout << "After CASE I";
      std::cout<<"triangle t= ("<<i<<","<<pid<<","<<k<<")"<<std::endl;
      calculate_weight(i, pid, k);
      count++;

      // first ordering
      processDomain(D1, e_D1.first, e_D1.second, count);
      // after the subdomains left and right have been processed
      assert(domain.has_islands());
      assert(domain.b_ids[0] == i);
      assert(domain.b_ids[domain.b_ids.size() - 1] == k);

      std::cout << "After CASE I";
      std::cout<<"triangle t= ("<<i<<","<<pid<<","<<k<<")"<<std::endl;
      calculate_weight(i, pid, k);      
      count++;

      std::cout << "--FINISHED with first ordering--, onto the SECOND" << std::endl;

      std::cin.get();


    }



    // create a new vector on which pid will run
    std::vector<int> third_verts;
    third_verts.reserve(domain.b_ids.size() - 2);
    // without the first and the last (source, target of access edge)
    third_verts.insert(third_verts.begin(), domain.b_ids.begin() +1, domain.b_ids.end() - 1);
    assert(third_verts.size() == domain.b_ids.size() - 2);


    // CASE II
    for(int pid : third_verts)
    {

      // avoid source & target of e_D
      //if(pid == i || pid == k)
      //  continue;

      // return if boundary is only 3 v. and has holes inside
      if(domain.b_ids.size() == 3 && domain.has_islands())
        return;

      // print triangle t
      //print_triangle(i, pid, k, out_domain, print_i, points);

      // split to two sub-domains
      Domain<PointRange> D1;
      Domain<PointRange> D2;
      // essentially splitting boundaries
      split_domain(domain, D1, D2, i, pid, k);
      // D1, D2 have just new boundaries - no hole information.

      assert(D1.b_ids[0] == i);
      assert(D2.b_ids[D2.b_ids.size() - 1] == k);


      // get new access edges for each
      std::pair<int, int> e_D1 = D1.get_access_edge();
      std::pair<int, int> e_D2 = D2.get_access_edge();


      std::cout << "splitting domain = ";
      for(int j=0; j<domain.b_ids.size(); ++j)
      {
        std::cout << domain.b_ids[j] << " ";
      }
      std::cout << " with pid= " << pid << std::endl;

      std::cout << "D1 = ";
      for(int j=0; j<D1.b_ids.size(); ++j)
      {
        std::cout << D1.b_ids[j] << " ";
      }
      std::cout << std::endl;

      std::cout << "D2 = ";
      for(int j=0; j<D2.b_ids.size(); ++j)
      {
        std::cout << D2.b_ids[j] << " ";
      }
      std::cout << std::endl;



      // assign all combination of holes to subdomains and process each pair
      Phi partition_space; // todo : pre-calculate this once.
      do_permutations(domain.holes_list, partition_space);
      if(partition_space.empty())
      {
        // when the domain has been merged so that there is no holes inside
        processDomain(D1, e_D1.first, e_D1.second, count);
        processDomain(D2, e_D2.first, e_D2.second, count);
      }
      else
      {

        std::cout << "--- entering PARTITION SPACE ---" << std::endl;
        // when t is formed with a vertex on the boundary of a domain with holes
        for(std::size_t p = 0; p < partition_space.size(); ++p)
        {
          std::vector<int> lholes = partition_space.lsubset(p);
          std::vector<int> rholes = partition_space.rsubset(p);

          D1.clear_holes();
          D2.clear_holes();

          for(int lh : lholes)
            D1.add_hole(domain.holes_list[lh]);

          for(int rh : rholes)
            D2.add_hole(domain.holes_list[rh]);


          if(D1.is_empty() && !D2.has_islands())
          {
            continue;
          }

          if(D2.is_empty() && !D1.has_islands())
          {
            continue;
          }


          processDomain(D1, e_D1.first, e_D1.second, count);
          processDomain(D2, e_D2.first, e_D2.second, count);
        }

      }

      // calculate weight of triangle t - after the subdomains left and right have been checked

      std::cout << "back to domain = ";
      for(int j=0; j<domain.b_ids.size(); ++j)
      {
        std::cout << domain.b_ids[j] << " ";
      }
      std::cout << " with pid= " << pid << std::endl;


      if(i == pid || k == pid)
      {
        std::cout << "aborting ("<<i<<","<<pid<<","<<k<<")"<<", manifold edge" << std::endl;
        continue;
      }

      std::cout<<"triangle t= ("<<i<<","<<pid<<","<<k<<")"<<std::endl;

      calculate_weight(i, pid, k);
      count++;

    }
  }


  bool are_vertices_in_island(const int& i, const int& m, const int& k)
  {

    std::vector<int>::iterator it1, it2, it3;

    it1 = std::find(init_island.begin(), init_island.end(), i);
    it2 = std::find(init_island.begin(), init_island.end(), m);
    it3 = std::find(init_island.begin(), init_island.end(), k);

    return (it1 != init_island.end()) && (it2 != init_island.end()) && (it3 != init_island.end()) ?  true : false;

  }


  bool are_vertices_on_boundary(const int& i, const int& m, const int& k)
  {

    std::vector<int>::iterator it1, it2, it3;

    it1 = std::find(init_b.begin(), init_b.end(), i);
    it2 = std::find(init_b.begin(), init_b.end(), m);
    it3 = std::find(init_b.begin(), init_b.end(), k);

    return (it1 != init_b.end()) && (it2 != init_b.end()) && (it3 != init_b.end()) ?  true : false;

  }


  void calculate_weight(int& i, int& m, int& k)
  {


    if(are_vertices_in_island(i, m, k))
    {
      std::cout << "vertices are all in island! no weight caclulated" << std::endl;
      return;
    }

    /* just for testing - should not be used
    if(are_vertices_on_boundary(i, m, k))
    {
      std::cout << "vertices are all on boundary! no weight caclulated" << std::endl;
      return;
    }
    */


    // i, m, k are global indices
    assert(m != i);
    assert(m != k);

   // if(i > k)
   // {
   //   std::swap(i, k); // needed to store the edge (i,k) sorted. Maybe move this in the Look_up_map.
   // }
    assert(i < k);


    if(i == 0 && k == 2 && m == 3)
    {
      std::cout << "stop" << std::endl;
    }

    PointRange Q;
    const Weight& w_imk = WC(points, Q, i, m, k, lambda);


    if(w_imk == Weight::NOT_VALID())
    {
      std::cerr << "non-manifold edge"  << std::endl;
      return;
    }

    auto lw = W.get(i,m);
    auto rw = W.get(m,k); // should swap! (3,2) is not in the table but (2,3) is. - or change the get in the table
    const Weight& w = W.get(i,m) + W.get(m,k) + w_imk;
    //const Weight& w = w_imk;

    if(lambda.get(i, k) == -1 || w < W.get(i, k)) {
      W.put(i, k, w);
      lambda.put(i, k, m);

      W.print("data/weight.dat");
      lambda.print("data/lambda.dat");
      std::cout << " value updated for ("<< i << " " << k << ") with: "
                << lambda.get(i, k)<< std::endl;
      std::cin.get();

    }
  }



  // data
  PointRange points;

  //std::set<std::vector<std::size_t>> memoized;

  WeightTable W;
  LambdaTable lambda;

  Domain<PointRange> domain;
  const WeightCalculator& WC;


  // initial boundary & island
  std::vector<int> init_b;
  std::vector<int> init_island;


  std::ofstream out_domain;

  int print_i;


};





} // namespace internal
} // namsepace CGAL





#endif // ISLAND_TRIANGULATE_HOLE_POLYLINE_H
