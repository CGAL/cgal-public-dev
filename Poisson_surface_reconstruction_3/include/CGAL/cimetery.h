//Poisson_reconstruction_function.h
  //grad calculation method
  bool m_gradfit;

  const bool gradfit() const
  {
    return m_gradfit;
  }

  bool& gradfit()
  {
    return m_gradfit;
  }

  void compute_grads()
   {
     if(m_gradfit){
       m_tr->compute_grad_fit();
       return;
     }
     m_tr->compute_grad_per_cell();
     m_tr->compute_grad_per_vertex();
   }

//Reconstruction_triangulation_3.h

  void compute_grad_fit()
  {
    for(auto it = this->finite_vertices_begin(); it != this->finite_vertices_end(); it++){
      Vertex_handle v = it;
      std::vector<Cell_handle> cells;
		  this->incident_cells(v, std::back_inserter(cells));

		  std::vector<Vertex_handle> init_vertices;
		  init_vertices.clear();
		  this->incident_vertices(v, std::back_inserter(init_vertices));//, init_vertices.end()));

		  std::set<Vertex_handle> vertices;
		  vertices.clear();
		  for(auto it = init_vertices.begin();
			  it != init_vertices.end();
			  it++)
		  {
			  vertices.insert(*it);
			  this->incident_vertices(*it, std::inserter(vertices, vertices.end()));
		  }

		  v->df() = grad_fit(this, vertices, v);
	  }
  }

  //convolution code

  void grad_convolution(){ //unweighted, repeats 100 times
    for(int i = 0; i < 20; i++)
    {
      for(auto it = this->finite_vertices_begin();
        it != this->finite_vertices_end(); it++)
      {
        std::vector<Vertex_handle> vertices;
        this->incident_vertices(it, std::back_inserter(vertices));
        Vector grad(0.0, 0.0, 0.0);
        for(auto v = vertices.begin(); v != vertices.end(); v++)
        {
          grad += (*v)->df();
        }
        grad /= vertices.size();
        it->df() = grad;
      }
    }
  }

  void grad_weighted_convolution(){ //weighted, repeats 20 times
    FT alpha = 10.0;
    for(int i = 0; i < 20; i++)
    {
      for(auto it = this->finite_vertices_begin();
      it != this->finite_vertices_end(); it++)
      {
        std::vector<Vertex_handle> vertices;
        this->incident_vertices(it, std::back_inserter(vertices));
        Vector grad(0.0, 0.0, 0.0);
        FT weights = 0.0;
        for(auto v = vertices.begin(); v != vertices.end(); v++)
        {
          Vector vec(it->point(), (*v)->point());
          FT dist = std::sqrt(vec * vec);
          FT weight = std::exp(-alpha * dist);
          grad += weight * (*v)->df();
          weights += weight;
        }
        grad /= weights;
        it->df() = grad;
      }
    }
  }

  //compute_grads function in CGAL_POISSON_RECONSTRUCTION_FUNCTION_H

  void compute_grads()
   {
     switch(m_smooth){
       case 1:
         m_tr->compute_grad_per_vertex();
         break;
       case 2:
         m_tr->compute_grad_bounding_sphere();
         break;
       default:
         break;
     }
   }

   //gradient calculation using sphere averaging for each vertex;

     void compute_grad_bounding_sphere()
     {
       for(auto it = this->finite_vertices_begin(); it != this->finite_vertices_end(); it++)
       {
         compute_grad_bounding_sphere(it);
       }
     }

  //gradient calculation using cell volume weighted averaging for each vertex:  
    void compute_grad_per_vertex()
    {
      for(auto it = this->finite_vertices_begin(); it != this->finite_vertices_end(); it++){
        compute_df(it);
      }
    }
