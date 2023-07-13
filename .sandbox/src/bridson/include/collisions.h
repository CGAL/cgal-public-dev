

template <class SimulationMeshRange,
          class OutputIterator,
          class FT>
OutputIterator do_collide(
    const SimulationMeshRange& range,
    OutputIterator out,
    FT time_step 
){

  typedef typename SimulationMeshRange::const_iterator SimulationMeshIterator;

  // Build AABB tree whose boxes bound the swept primitives of each mesh


  bool report_overlap =  choose_parameter(get_parameter(np, internal_np::overlap_test),false);

  typedef CGAL::Box_intersection_d::ID_FROM_BOX_ADDRESS Box_policy;
  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, SimulationMeshIterator, Box_policy> Mesh_box;

  std::vector<Mesh_box> boxes;
  boxes.reserve(std::distance(range.begin(), range.end()));

  for(SimulationMeshIterator it = range.begin(); it != range.end(); ++it)
  {
    boxes.push_back( Mesh_box(Polygon_mesh_processing::bbox(*it), it) );
  }

  std::vector<Mesh_box*> boxes_ptr(boost::make_counting_iterator(&boxes[0]),
                                   boost::make_counting_iterator(&boxes[0]+boxes.size()));

  typedef typename boost::range_value<NamedParametersRange>::type NP_rng;
  typedef typename boost::range_value<TriangleMeshRange>::type TriangleMesh;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters, NP_rng>::type GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  //get all the pairs of meshes intersecting (no strict inclusion test)
  std::ptrdiff_t cutoff = 2000;
  internal::Mesh_callback<TriangleMeshRange, GT, OutputIterator, NamedParametersRange> callback(range, out, report_overlap, gt, nps);
  CGAL::box_self_intersection_d(boxes_ptr.begin(), boxes_ptr.end(), callback, cutoff);
  return callback.m_iterator;
}