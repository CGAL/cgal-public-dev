
namespace CGAL {

/*!
\ingroup PkgBasicViewerClasses

The class `Graphics_scene_options` is used to tune the way that the cells of a given data structure of \cgal are considered.
The different `std::function` can be modified to change for example the behavior of the drawing.
`VolumeDescriptor` can be `void` for data structures that do not represent volumes.

This class is a model of `GraphicsSceneOptions` when `VolumeDescriptor` is `void`, or a model of `GraphicsSceneOptionsWithVolumes` otherwise (`VolumeDescriptor` non `void`).

\tparam DS a data structure of \cgal.
\tparam VertexDescriptor a descriptor of vertices of `DS`.
\tparam EdgeDescriptor a descriptor of edges of `DS`.
\tparam FaceDescriptor a descriptor of faces of `DS`.
\tparam VolumeDescriptor a descriptor of volumes of `DS`. `void` by default.

\cgalModels{GraphicsSceneOptions or GraphicsSceneOptionsWithVolumes}
*/

template <typename DS,
          typename VertexDescriptor,
          typename EdgeDescriptor,
          typename FaceDescriptor,
          typename VolumeDescriptor=void>
struct Graphics_scene_options
{
public:
  /// constructs default graphics scene options.
  ///
  /// Vertices, edges, and faces are enabled. The `draw_vertex`,
  /// `draw_edge`, and `draw_face` callbacks return `true` by default,
  /// while the `colored_vertex`, `colored_edge`, `colored_face`, and
  /// `face_wireframe` callbacks return `false` by default.
  Graphics_scene_options();

  typedef VertexDescriptor vertex_descriptor;
  typedef EdgeDescriptor edge_descriptor;
  typedef FaceDescriptor face_descriptor;
  typedef VolumeDescriptor volume_descriptor;

};

} // End namespace CGAL
