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
