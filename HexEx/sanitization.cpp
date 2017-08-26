#include "hexextr.h"

void HexExtr::sanitize(){/**
* Sanitizes the parametrization so as to prevent numerical inaccuracies.
*/
  check_singularity();
  truncate_precision();
}

void HexExtr::check_singularity(){ /**
* if an edge is singular, all the darts incident on the edge are set to indicate that the edge is singular. The singular_edges keep a count of the number of darts incident on the singular edge
*/
  for(LCC_3::One_dart_per_cell_range<1>::iterator e = input_tet_mesh.one_dart_per_cell<1>().begin(), eend = input_tet_mesh.one_dart_per_cell<1>().end(); e != eend; e++){
    if(is_singular(e)){
      int sing_edge = 1;
      for(LCC_3::Dart_of_cell_range<1>::iterator d = input_tet_mesh.darts_of_cell<1>(e).begin(), dend = input_tet_mesh.darts_of_cell<1>(e).end(); d != dend; d++){
        if((input_tet_mesh.info(d)).singular == true) sing_edge++;
      }
      for(LCC_3::Dart_of_cell_range<1>::iterator d = input_tet_mesh.darts_of_cell<1>(e).begin(), dend = input_tet_mesh.darts_of_cell<1>(e).end(); d != dend; d++){
        (input_tet_mesh.info(d)).singular = true; (input_tet_mesh.info(d)).singular_edges = sing_edge; 
      }
    } 
  }
}

void HexExtr::truncate_precision(){/**
* Takes the maximum coordinates of the parametrization, finds delta by taking ceil of raising the 2 to log(max), thus getting the upper limit of max in binary system.
*/
  
  for(LCC_3::One_dart_per_cell_range<0>::iterator v=input_tet_mesh.one_dart_per_cell<0>().begin(), vend=input_tet_mesh.one_dart_per_cell<0>().end(); v!=vend; ++v){
    double max = std::numeric_limits<double>::min();
    Point p; 
    for(LCC_3::Dart_of_cell_range<0>::iterator dh = input_tet_mesh.darts_of_cell<0>(v).begin(), dhend = input_tet_mesh.darts_of_cell<0>(v).end(); dh != dhend; dh++){
      p = (input_tet_mesh.info(dh)).parameters; 
      max = std::max(max, std::abs(p[0])); max = std::max(max, std::abs(p[1])); max = std::max(max, std::abs(p[2]));
    }
    p =  (input_tet_mesh.info(v)).parameters;
    double delta = std::pow(2, std::ceil(std::log(max)/std::log(2)));

    std::vector<double> temp(3);
    for(unsigned int i = 0; i < 3; ++i){
      int sign = std::signbit(p[i]) ? -1 : 1;
      temp[i] = p[i] + sign*delta;
      temp[i] = temp[i] - sign*delta;
    }
//Now store tmp into parameters. This "sanitizes" the parameters. 
    ((input_tet_mesh.info(v)).parameters) = Point(temp[0], temp[1], temp[2]);
    
//singularity check
    if((input_tet_mesh.info(v)).singular == true){
      fix_singularity(v);
    }
    propagate_parameters(v);  
  }

}

void HexExtr::propagate_parameters(Dart_handle& dh){/**
* iterates through all the vertices incident on a face and updates the new parameters of the vertices on opposite face using the transformatuions in g, thus propagating the updates parameters. 
*/
  Point p = (input_tet_mesh.info(dh)).parameters;
  for(LCC_3::One_dart_per_incident_cell_range<2,0>::iterator f = input_tet_mesh.one_dart_per_incident_cell<2,0>(dh).begin(), fend = input_tet_mesh.one_dart_per_incident_cell<2,0>(dh).end(); f!=fend; f++){// for the vertex in each cell
    (input_tet_mesh.info(f)).parameters = p;
    if(input_tet_mesh.is_free<2>(f)) continue;
    else{
      Dart_handle fopp = input_tet_mesh.opposite(input_tet_mesh.beta(f, 2));

      //update the parameter 
      (input_tet_mesh.info(fopp)).parameters = p.transform(g[(input_tet_mesh.info(f)).cell_no][(input_tet_mesh.info(fopp)).cell_no]);
    }    
  }
}

Aff_transformation HexExtr::get_transition(Dart_handle& dh, Dart_handle& singular_edge){/**
* if two darts belong to the same tet, return identity. If not, look up the transition matrix from g iteratively for all faces around an edge, multiply them and return.
*/

// If dh and singular_edge belong to the same tet, transition function is identity
  if((input_tet_mesh.info(dh)).cell_no == (input_tet_mesh.info(singular_edge).cell_no))
    return Aff_transformation(1,0,0,0,1,0,0,0,1,1); //return identity

  else{
    Aff_transformation temp(1,0,0,0,1,0,0,0,1,1);
    for(LCC_3::One_dart_per_incident_cell_range<2,1>::iterator f = input_tet_mesh.one_dart_per_incident_cell<2,1>(dh).begin(),
fend = input_tet_mesh.one_dart_per_incident_cell<2,1>(dh).end(); f != fend; f++){
      if(!input_tet_mesh.is_free<3>(f)){// iterate through all faces incident on an edge
        if((input_tet_mesh.info(dh)).cell_no == (input_tet_mesh.info(singular_edge).cell_no)) return temp;
        else temp = temp * g[(input_tet_mesh.info(f)).cell_no][(input_tet_mesh.info(input_tet_mesh.beta(f, 3))).cell_no]; 
      }
    }
    return temp;   
  }  
}

bool HexExtr::is_singular(Dart_handle& e){// is the edge singular?

//intializing identity and temp to identity transformation:
  Aff_transformation temp(1, 0, 0, 0, 1, 0, 0, 0, 1, 1), identity(1, 0, 0, 0, 1, 0, 0, 0, 1, 1);

//iterate through all the faces incident on the edge, and multiply the transition matrices between two adjacent tets to the get the trasition matrix around the edge:
  for(LCC_3::One_dart_per_incident_cell_range<2,1>::iterator f = input_tet_mesh.one_dart_per_incident_cell<2,1>(e).begin(),
fend = input_tet_mesh.one_dart_per_incident_cell<2,1>(e).end(); f != fend; f++){
      if(!input_tet_mesh.is_free<3>(f)){
        temp = temp * g[(input_tet_mesh.info(f)).cell_no][(input_tet_mesh.info(input_tet_mesh.beta(f, 3))).cell_no]; 
      }
    }

// if the transition matrix is identity, the edge is not singular. If not, it is singular:
  return !(temp.m(0,0) == 1 && temp.m(0,1) == 0 && temp.m(0,2) == 0 && temp.m(0,3) == 0
&&  temp.m(1,0) == 0 &&temp.m(1,1) == 1 &&temp.m(1,2) == 0 && temp.m(1,3) == 0 && temp.m(2,0) == 0 
&& temp.m(2,1) == 0 && temp.m(2,2) == 1 && temp.m(2,3) == 0); 
}


Dart_handle HexExtr::get_singular_edge(Dart_handle& v){//get the edge incident to given dart which is singular.
  for(LCC_3::One_dart_per_incident_cell_range<1,0>::iterator e = input_tet_mesh.one_dart_per_incident_cell<1,0>(v).begin(), eend = input_tet_mesh.one_dart_per_incident_cell<1,0>(v).end(); e != eend; e++){ //iterate through all the edges incident on the vertex to get the singlar edge
    if(is_singular(e)) return e;
  }
  return Dart_handle();
}

Point_3 HexExtr::get_projected_parameters(Dart_handle& dh, Dart_handle& singular_edge){/**
* given an edge defined by dart handle dh, we find the direction to which it points in dir. We get the transition matrix between dh and singular edge and apply the transformation */

//finding from and to points of dh 
  Point_3 from = (input_tet_mesh.info(dh)).parameters, to = (input_tet_mesh.info(input_tet_mesh.opposite(dh))).parameters;

//this gives the direction:
  Vector_3 dir(from, to);

  Aff_transformation at= get_transition(dh, singular_edge);
  dir = dir.transform(at);
  Point_3 p(from[0]+dir[0], from[1]+dir[1], from[2]+dir[2]);
  return p;
}

void HexExtr::fix_singularity(Dart_handle& dh){ //fixes singular edges by modifying parameters given the deviation is less than e-4.
  Point_3 param = (input_tet_mesh.info(dh)).parameters; Point_3 p;

  if(input_tet_mesh.info(dh).singular_edges == 2){
    Dart_handle singular_edge = get_singular_edge(dh);
    p = get_projected_parameters(dh, singular_edge);
 }
  else{
     p = Point_3(std::round(param[0]), std::round(param[1]), std::round(param[2]));
  }
  Vector_3 v(param, p);
  double eps = 1e-4;
  if(v.squared_length() < eps)
    (input_tet_mesh.info(dh)).parameters = p;
}


