#include "config.h"
#include <CGAL/Three/Scene_group_item.h>
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/Real_timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/IO/PLY_writer.h>

#include <CGAL/Levels_of_detail.h>
#include <CGAL/Random.h>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <boost/function_output_iterator.hpp>

#include "ui_Level_of_detail_plugin.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QMessageBox>

using namespace CGAL::Three;

typedef std::map<int, CGAL::Level_of_detail::Semantic_label> Map_l2sl;
typedef boost::shared_ptr<Map_l2sl> Map_l2sl_ptr;

struct Semantic_map_from_labels
{
  typedef Point_set::Index key_type;
  typedef CGAL::Level_of_detail::Semantic_label value_type;
  typedef CGAL::Level_of_detail::Semantic_label reference;
  typedef boost::readable_property_map_tag category;

  Point_set* points;
  Point_set::Property_map<int> label_map;
  Map_l2sl_ptr map_l2sl;

  Semantic_map_from_labels () { }
  Semantic_map_from_labels (Point_set* points) : points (points)
                                               , map_l2sl (new Map_l2sl())
  {
    label_map = points->property_map<int>("label").first;
  }

  friend value_type get (const Semantic_map_from_labels& map, const key_type& key)
  {
    int l = map.label_map[key];

    typename Map_l2sl::const_iterator
      found = map.map_l2sl->find(l);
    if (found == map.map_l2sl->end())
      return CGAL::Level_of_detail::Semantic_label::UNASSIGNED;

    return found->second;
  }
};

struct Insert_point_colored_by_index
{
  typedef std::pair<Kernel::Point_3, int> argument_type;
  typedef void result_type;
          
  Point_set& points;
  Point_set::Property_map<unsigned char> red, green, blue;

  Insert_point_colored_by_index (Point_set& points) : points (points)
  {
    red = points.template add_property_map<unsigned char>("r", 0).first;
    green = points.template add_property_map<unsigned char>("g", 0).first;
    blue = points.template add_property_map<unsigned char>("b", 0).first;
    points.check_colors();
  }

  void operator() (const argument_type& a)
  {
    Point_set::iterator pt = points.insert (a.first);
    if (a.second == -1)
      return;
    
    CGAL::Random rand(a.second);
    red[*pt] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
    green[*pt] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
    blue[*pt] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
  }
};

struct Add_polyline_from_segment
{
  typedef Kernel::Segment_3 argument_type;
  typedef void result_type;

  Scene_polylines_item::Polylines_container& polylines;

  Add_polyline_from_segment (Scene_polylines_item::Polylines_container& polylines) : polylines (polylines) { }

  void operator() (const argument_type& a)
  {
    polylines.push_back (Scene_polylines_item::Polyline());
    polylines.back().push_back (a.source());
    polylines.back().push_back (a.target());
  }
};

struct Add_polygon_with_visibility_color
{
  typedef std::pair<std::vector<std::size_t>, CGAL::Level_of_detail::Visibility_label> argument_type;
  typedef void result_type;

  std::vector<std::vector<std::size_t> >& polygons;
  std::vector<CGAL::Color>& fcolors;

  Add_polygon_with_visibility_color (std::vector<std::vector<std::size_t> >& polygons,
                                     std::vector<CGAL::Color>& fcolors)
    : polygons (polygons), fcolors (fcolors) { }

  void operator() (const argument_type& a)
  {
    polygons.push_back (a.first);
    
    unsigned char r, g, b;
    if (a.second == CGAL::Level_of_detail::Visibility_label::OUTSIDE)
    {
      r = 186; g = 189; b = 182;
    }
    else if (a.second == CGAL::Level_of_detail::Visibility_label::INSIDE)
    {
      r = 245; g = 121; b = 0;
    }
    else
    {
      r = 77; g = 131; b = 186;
    }

    fcolors.push_back (CGAL::Color(r,g,b));
  }

};
        
struct Add_triangle_with_building_color
{
  typedef std::pair<CGAL::cpp11::array<std::size_t, 3>, int> argument_type;
  typedef void result_type;

  std::vector<std::vector<std::size_t> >& polygons;
  std::vector<CGAL::Color>& fcolors;

  Add_triangle_with_building_color (std::vector<std::vector<std::size_t> >& polygons,
                                    std::vector<CGAL::Color>& fcolors)
    : polygons (polygons), fcolors (fcolors) { }

  void operator() (const argument_type& a)
  {
    polygons.push_back (std::vector<std::size_t>(3));
    for (std::size_t i = 0; i < 3; ++ i)
      polygons.back()[i] = a.first[i];
    
    unsigned char r, g, b;

    if (a.second < 0)
    {
      r = 128; g = 128; b = 128;
    }
    else
    {
      CGAL::Random rand(a.second);
      r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
    }
    
    fcolors.push_back (CGAL::Color(r,g,b));
  }

};
        
struct Add_triangle_with_metaface_id
{
  typedef std::pair<CGAL::cpp11::array<std::size_t, 3>, int> argument_type;
  typedef void result_type;

  std::vector<Kernel::Point_3>& vertices;
  std::vector<std::vector<std::size_t> >& polygons;
  Scene_polylines_item::Polylines_container& polylines;

  typedef std::map<std::pair<std::size_t, std::size_t>, int> Map_edges;
  Map_edges map_edges;

  Add_triangle_with_metaface_id (std::vector<Kernel::Point_3>& vertices,
                                 std::vector<std::vector<std::size_t> >& polygons,
                                 Scene_polylines_item::Polylines_container& polylines)
    : vertices (vertices), polygons (polygons), polylines (polylines) { }

  void operator() (const argument_type& a)
  {
    polygons.push_back (std::vector<std::size_t>(3));
    for (std::size_t i = 0; i < 3; ++ i)
    {
      polygons.back()[i] = a.first[i];


      // handle edge
      std::size_t idx_a = a.first[i];
      std::size_t idx_b = a.first[(i+1)%3];
      if (idx_a > idx_b)
        std::swap (idx_a, idx_b);
      
      typename Map_edges::iterator iter;
      bool inserted = false;
      boost::tie (iter, inserted) = map_edges.insert (std::make_pair (std::make_pair (idx_a, idx_b), a.second));
      if (!inserted && iter->second != a.second) // edge between two metafaces
      {
        polylines.push_back (Scene_polylines_item::Polyline());
        polylines.back().push_back (vertices[idx_a]);
        polylines.back().push_back (vertices[idx_b]);
      }
    }
  }
};
        
typedef CGAL::Level_of_detail::Level_of_detail<Kernel, Point_set, Point_set::Point_map,
                                               Semantic_map_from_labels,
                                               CGAL::Level_of_detail::Visibility_from_semantic_map<Semantic_map_from_labels>,
                                               CGAL::Tag_true> LOD;


class Polyhedron_demo_level_of_detail_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  QAction* actionLOD;
  
public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    mw = mainWindow;
    actionLOD = new QAction(tr("Level Of Detail"), mainWindow);
    actionLOD->setObjectName("actionLOD");
    autoConnectActions();
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionLOD;
  }
  
  //! Applicable if the currently selected item is a
  //! points_with_normal_item.
  bool applicable(QAction*) const {
    Scene_points_with_normal_item* item =
      qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if (!item)
      return false;
    if (!item->point_set()->property_map<int>("label").second)
      return false;
    return true;
  }

public Q_SLOTS:
  void on_actionLOD_triggered();
}; // end Polyhedron_demo_level_of_detail_plugin

class Polyhedron_demo_lod_dialog : public QDialog, private Ui::LevelOfDetailDialog
{
  Q_OBJECT

  Point_set* points;
public:
  
  Polyhedron_demo_lod_dialog(Point_set* points, QWidget * /*parent*/ = 0)
    : points (points)
  {
    setupUi(this);
    compute_detailed_parameters();

    connect(scaleDoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(compute_detailed_parameters()));
    connect(minimumWallLengthDoubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(compute_detailed_parameters()));
  }

  // Building Detection
  double alpha_shape_size() const { return alphaShapeSizeDoubleSpinBox->value(); }
  double cluster_radius() const { return clusterRadiusDoubleSpinBox->value(); }
  double building_grid_cell_width() const { return gridCellWidthDoubleSpinBox->value(); }
  double max_regularized_angle() const
  {
    return std::cos (CGAL_PI * maximumRegularizedAngleDoubleSpinBox->value() / 180.);
  }
  double max_regularized_gap() const { return maximumRegularizedGapDoubleSpinBox->value(); }
  double noise_level() const { return noiseLevelDoubleSpinBox->value(); }
  double normal_threshold() const
  {
    return std::cos (CGAL_PI * normalThresholdDoubleSpinBox->value() / 180.);
  }
  double min_wall_length() const { return minimumWallLengthDoubleSpinBox_2->value(); }

  // Tree Detection
  double tree_grid_cell_width() const { return gridCellWidthDoubleSpinBox_2->value(); }
  double min_tree_height() const { return minimumTreeHeightDoubleSpinBox->value(); }
  double min_tree_radius() const { return minimumTreeRadiusDoubleSpinBox->value(); }

  // Partitionning
  double min_face_width() const { return minimumFaceWidthDoubleSpinBox->value(); }
  std::size_t num_intersections() const { return numberOfIntersectionsPerLineSpinBox->value(); }
  bool make_consistent_visibility() const { return makeLocallyConsistentVisibilityCheckBox->isChecked(); }

  // Footprints and 3D Models
  std::size_t num_faces() const { return minimumNumberOfFacesPerBuildingSpinBox->value(); }
  double segment_constraints_threshold() const { return segmentConstraintsThresholdDoubleSpinBox->value(); }
  double smooth_ground_precision() const { return smoothGroundPrecisionDoubleSpinBox->value(); }
  double tree_precision() const { return treePrecisionDoubleSpinBox->value(); }

  bool lod0() const { return lod0CheckBox->isChecked(); }
  bool lod1() const { return lod1CheckBox->isChecked(); }
  bool lod2() const { return lod2CheckBox->isChecked(); }
  bool detailed() const { return detailedOutput->isChecked(); }

  QString ground_indices() const { return groundLineEdit->text(); }
  QString boundary_indices() const { return buildingBoundariesLineEdit->text(); }
  QString interior_indices() const { return buildingInteriorLineEdit->text(); }
  QString vegetation_indices() const { return vegetationLineEdit->text(); }

  void set_ground_index (const int& i) const
  {
    return groundLineEdit->setText(tr("%1 %2").arg(groundLineEdit->text()).arg(i));
  }
  void set_boundary_index (const int& i) const { return buildingBoundariesLineEdit->setText(tr("%1").arg(i)); }
  void set_interior_index (const int& i) const { return buildingInteriorLineEdit->setText(tr("%1").arg(i)); }
  void set_vegetation_index (const int& i) const
  {
    return vegetationLineEdit->setText(tr("%1 %2").arg(vegetationLineEdit->text()).arg(i));
  }

  QTextEdit* comment_section() { return textEdit; }

public Q_SLOTS:
  
  void compute_detailed_parameters()
  {
    double scale = scaleDoubleSpinBox->value();
    double wall_length = minimumWallLengthDoubleSpinBox->value();
    
    alphaShapeSizeDoubleSpinBox->setValue(0.5 * scale);
    noiseLevelDoubleSpinBox->setValue(scale);
    clusterRadiusDoubleSpinBox->setValue(scale);
    minimumWallLengthDoubleSpinBox_2->setValue(wall_length);
    gridCellWidthDoubleSpinBox->setValue(0.25 * scale);
    
    gridCellWidthDoubleSpinBox_2->setValue(scale);
    minimumTreeHeightDoubleSpinBox->setValue(1.5 * wall_length);
    minimumTreeRadiusDoubleSpinBox->setValue(wall_length);

    minimumFaceWidthDoubleSpinBox->setValue(0.5 * scale);

    segmentConstraintsThresholdDoubleSpinBox->setValue(0.5 * scale);
    smoothGroundPrecisionDoubleSpinBox->setValue(scale);
    treePrecisionDoubleSpinBox->setValue(scale);
  }
};


struct array_to_vector
{
  std::vector<std::vector<std::size_t> >& vectors;

  array_to_vector (std::vector<std::vector<std::size_t> >& vectors) : vectors (vectors) { }

  void operator() (const CGAL::cpp11::array<std::size_t, 3>& ar)
  {
    vectors.push_back (std::vector<std::size_t>(3));
    vectors.back()[0] = ar[0];
    vectors.back()[1] = ar[1];
    vectors.back()[2] = ar[2];
  }
};

void Polyhedron_demo_level_of_detail_plugin::on_actionLOD_triggered()
{
  CGAL::Random rand(time(0));
  
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    // Gets point set
    Point_set* points = item->point_set();
    if(points == NULL)
        return;
    Polyhedron_demo_lod_dialog dialog (points);

    // prefill dialog if label info found
    const std::string& comments = item->comments();
    std::istringstream iss (comments);
    std::string word;
    while (iss >> word)
    {
      if (word == "label")
      {
        int number;
        if (iss >> number)
        {
          std::string label;
          if (iss >> label)
          {
            if (label == "ground" || label == "bridge_deck" || label == "water")
              dialog.set_ground_index(number);
            else if (label == "facade" || label == "wall")
              dialog.set_boundary_index(number);
            else if (label == "roof" || label == "building")
              dialog.set_interior_index(number);
            else if (label == "vegetation" ||
//                     label == "low_veget" ||
                     label == "med_veget" ||
                     label == "high_veget" ||
                     label == "tree")
              dialog.set_vegetation_index(number);
          }
        }
      }
    }
    
    dialog.comment_section()->setText(comments.c_str());
    if(!dialog.exec())
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    QApplication::processEvents();
    
    CGAL::Real_timer task_timer; task_timer.start();

    Semantic_map_from_labels semantic_map (points);

    std::istringstream gi (dialog.ground_indices().toStdString());
    std::istringstream bi (dialog.boundary_indices().toStdString());
    std::istringstream ii (dialog.interior_indices().toStdString());
    std::istringstream vi (dialog.vegetation_indices().toStdString());
    int idx;
    while (gi >> idx)
    {
      std::cerr << idx << " is ground" << std::endl;
      semantic_map.map_l2sl->insert
        (std::make_pair (idx, CGAL::Level_of_detail::Semantic_label::GROUND));
    }
    while (bi >> idx)
    {
      std::cerr << idx << " is building boundary" << std::endl;
      semantic_map.map_l2sl->insert
        (std::make_pair (idx, CGAL::Level_of_detail::Semantic_label::BUILDING_BOUNDARY));
    }
    while (ii >> idx)
    {
      std::cerr << idx << " is building interior" << std::endl;
      semantic_map.map_l2sl->insert
        (std::make_pair (idx, CGAL::Level_of_detail::Semantic_label::BUILDING_INTERIOR));
    }
    while (vi >> idx)
    {
      std::cerr << idx << " is vegetation" << std::endl;
      semantic_map.map_l2sl->insert
        (std::make_pair (idx, CGAL::Level_of_detail::Semantic_label::VEGETATION));
    }

    LOD lod (*points, points->point_map(), semantic_map);

    Scene_group_item* group = NULL;
    if (dialog.detailed())
    {
      group = new Scene_group_item(tr("%1 (LOD detailed output)").arg(item->name()));
      group->setVisible(false);
      scene->addItem(group);
    }
      
    lod.compute_planar_ground();
				
    lod.detect_building_boundaries(dialog.alpha_shape_size(), // alpha shape size
                                   dialog.noise_level(), // region growing epsilon
                                   dialog.cluster_radius(), // region growing cluster epsilon
                                   dialog.normal_threshold(), // region growing normal threshold
                                   dialog.min_wall_length(), // region growing min wall length
                                   dialog.max_regularized_angle(),
                                   dialog.max_regularized_gap(),
                                   dialog.building_grid_cell_width()); // grid cell width

    if (dialog.detailed())
    {
      Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item;
      new_item->setName("Boundary points");

      lod.output_filtered_boundary_points (new_item->point_set()->point_back_inserter());

      new_item->setVisible(false);
      new_item->setColor (Qt::red);
      scene->addItem(new_item);
      scene->changeGroup(new_item, group);
    }

    if (dialog.detailed())
    {
      Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item;
      new_item->setName("Detected 2D regions (points)");

      Insert_point_colored_by_index inserter (*(new_item->point_set()));

      lod.output_segmented_boundary_points (boost::make_function_output_iterator(inserter));

      new_item->setVisible(false);
      new_item->setColor (Qt::black);
      scene->addItem(new_item);
      scene->changeGroup(new_item, group);
    }

    if (dialog.detailed())
    {
      Scene_polylines_item* new_item = new Scene_polylines_item;
      new_item->setName("Detected 2D regions (segments)");

      Add_polyline_from_segment adder (new_item->polylines);
      lod.output_boundary_edges (boost::make_function_output_iterator(adder));
        
      new_item->setVisible(false);
      new_item->setColor (Qt::black);
      scene->addItem(new_item);
      scene->changeGroup(new_item, group);
    }
      
    lod.detect_trees(dialog.tree_grid_cell_width(),
                     dialog.min_tree_height(),
                     dialog.min_tree_radius());
    
    if (dialog.detailed())
    {
      Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item;
      new_item->setName("Segmented trees (points)");

      Insert_point_colored_by_index inserter (*(new_item->point_set()));

      lod.output_tree_points (boost::make_function_output_iterator(inserter));

      new_item->setVisible(false);
      new_item->setColor (Qt::black);
      scene->addItem(new_item);
      scene->changeGroup(new_item, group);
    }

    lod.partition(dialog.min_face_width(),
                  dialog.num_intersections(),
                  dialog.make_consistent_visibility());

    if (dialog.detailed())
    {
      Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item;
      new_item->setName("Partitioning");

      std::vector<Kernel::Point_3> vertices;
      std::vector<std::vector<std::size_t> > polygons;
      std::vector<CGAL::Color> fcolors;
      std::vector<CGAL::Color> vcolors;

      lod.output_partition_to_polygon_soup (std::back_inserter (vertices),
                                            std::back_inserter (polygons));
      for (std::size_t i = 0; i < polygons.size(); ++ i)
      {
        unsigned char r, g, b;
        r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        fcolors.push_back (CGAL::Color(r,g,b));
      }

      new_item->load (vertices, polygons, fcolors, vcolors);
      new_item->setVisible(false);
      scene->addItem(new_item);
      scene->changeGroup(new_item, group);

      std::ofstream out("debug.ply", std::ios::binary);
      out.precision (std::numeric_limits<double>::digits10 + 2);
      CGAL::write_PLY (out, new_item->points(), new_item->polygons());
    }

    if (dialog.detailed())
    {
      Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item;
      new_item->setName("Visibility");

      std::vector<Kernel::Point_3> vertices;
      std::vector<std::vector<std::size_t> > polygons;
      std::vector<CGAL::Color> fcolors;
      std::vector<CGAL::Color> vcolors;

      Add_polygon_with_visibility_color adder (polygons, fcolors);

      lod.output_partition_with_visibility_to_polygon_soup
        (std::back_inserter (vertices),
         boost::make_function_output_iterator(adder));

      new_item->load (vertices, polygons, fcolors, vcolors);
      new_item->setVisible(false);
      scene->addItem(new_item);
      scene->changeGroup(new_item, group);
    }
      
    lod.compute_footprints(dialog.segment_constraints_threshold(),
                           dialog.num_faces());

    if (dialog.detailed())
    {
      Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item;
      new_item->setName("Buildings");

      std::vector<Kernel::Point_3> vertices;
      std::vector<std::vector<std::size_t> > polygons;
      std::vector<CGAL::Color> fcolors;
      std::vector<CGAL::Color> vcolors;

      Add_triangle_with_building_color adder (polygons, fcolors);

      lod.output_building_footprints_to_triangle_soup
        (std::back_inserter (vertices),
         boost::make_function_output_iterator(adder));
        
      new_item->load (vertices, polygons, fcolors, vcolors);
      new_item->setRenderingMode(Flat);
      new_item->setVisible(false);
      scene->addItem(new_item);
      scene->changeGroup(new_item, group);
    }

    if (dialog.detailed())
    {
      Scene_polylines_item* new_item = new Scene_polylines_item;
      new_item->setName("Building walls");

      Add_polyline_from_segment adder (new_item->polylines);
      lod.output_building_footprints_to_segment_soup (boost::make_function_output_iterator(adder));
        
      new_item->setVisible(false);
      new_item->setColor (Qt::black);
      scene->addItem(new_item);
      scene->changeGroup(new_item, group);
    }

    std::vector<Kernel::Point_3> vertices;
    std::vector<std::vector<std::size_t> > polygons;
    std::vector<CGAL::Color> fcolors;
    std::vector<CGAL::Color> vcolors;
    
    std::size_t first_building_facet;
    std::size_t first_vegetation_facet;
    std::size_t first_wall_facet;
    
    if (dialog.lod0())
    {
      Scene_polygon_soup_item* lod0_item = new Scene_polygon_soup_item;
      lod0_item->setVisible(false);
      
      boost::tie (first_building_facet, first_vegetation_facet)
        = lod.output_lod0_to_triangle_soup
        (std::back_inserter (vertices),
         boost::make_function_output_iterator (array_to_vector(polygons)));
    
      // Fill colors according to facet type
      for (std::size_t i = 0; i < first_building_facet; ++ i)
        fcolors.push_back (CGAL::Color(186, 189, 182));
      for (std::size_t i = first_building_facet; i < first_vegetation_facet; ++ i)
        fcolors.push_back (CGAL::Color(245, 121, 0));
      for (std::size_t i = first_vegetation_facet; i < polygons.size(); ++ i)
        fcolors.push_back (CGAL::Color(138, 226, 52));
    
      lod0_item->load (vertices, polygons, fcolors, vcolors);
      lod0_item->setName(tr("%1 (LOD0)").arg(item->name()));
      lod0_item->setRenderingMode(Flat);
      lod0_item->setVisible (false);

      scene->addItem(lod0_item);
    }
    
    lod.extrude_footprints();
      
    lod.compute_smooth_ground(dialog.smooth_ground_precision());

    if (dialog.detailed())
    {
      Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item;
      new_item->setName("LOD1 metafaces");

      std::vector<Kernel::Point_3> vertices;
      std::vector<std::vector<std::size_t> > polygons;
      std::vector<CGAL::Color> fcolors;
      std::vector<CGAL::Color> vcolors;

      Add_triangle_with_building_color adder (polygons, fcolors);

      lod.output_lod1_to_triangle_soup
        (std::back_inserter (vertices),
         boost::make_function_output_iterator(adder));
        
      new_item->load (vertices, polygons, fcolors, vcolors);
      new_item->setRenderingMode(Flat);
      new_item->setVisible(false);
      scene->addItem(new_item);
      scene->changeGroup(new_item, group);
    }

    if (dialog.lod1())
    {
      Scene_group_item* lod1_item = new Scene_group_item(tr("%1 (LOD1)").arg(item->name()));
      lod1_item->setVisible(false);
      scene->addItem(lod1_item);
    
      Scene_polygon_soup_item* lod1_faces = new Scene_polygon_soup_item;
      lod1_faces->setName("Faces");
    
      Scene_polylines_item* lod1_edges = new Scene_polylines_item;
      lod1_edges->setName("Metaedges");
      lod1_edges->setColor (Qt::black);

      vertices.clear();
      polygons.clear();
      fcolors.clear();
      vcolors.clear();

      Add_triangle_with_metaface_id adder (vertices, polygons, lod1_edges->polylines);

      std::tie (first_building_facet, first_wall_facet, first_vegetation_facet)
        = lod.output_lod1_to_triangle_soup
        (std::back_inserter (vertices),
         boost::make_function_output_iterator (adder));

      std::cerr << vertices.size() << " ; " << polygons.size() << " ; " << lod1_edges->polylines.size() << std::endl;
    
      // Fill colors according to facet type
      for (std::size_t i = 0; i < first_building_facet; ++ i)
        fcolors.push_back (CGAL::Color(186, 189, 182));
      for (std::size_t i = first_building_facet; i < first_wall_facet; ++ i)
        fcolors.push_back (CGAL::Color(245, 121, 0));
      for (std::size_t i = first_wall_facet; i < first_vegetation_facet; ++ i)
        fcolors.push_back (CGAL::Color(77, 131, 186));
      for (std::size_t i = first_vegetation_facet; i < polygons.size(); ++ i)
        fcolors.push_back (CGAL::Color(138, 226, 52));

      lod1_faces->load (vertices, polygons, fcolors, vcolors);
      lod1_faces->setRenderingMode(Flat);
      scene->addItem(lod1_faces);
      scene->changeGroup(lod1_faces, lod1_item);
      scene->addItem(lod1_edges);
      scene->changeGroup(lod1_edges, lod1_item);
    }
    
    lod.fit_tree_models(dialog.tree_precision());

    if (dialog.lod2())
    {
      Scene_group_item* lod2_item = new Scene_group_item(tr("%1 (LOD2)").arg(item->name()));
      lod2_item->setVisible(false);
      scene->addItem(lod2_item);
    
      Scene_polygon_soup_item* lod2_faces = new Scene_polygon_soup_item;
      lod2_faces->setName("Faces");
    
      Scene_polylines_item* lod2_edges = new Scene_polylines_item;
      lod2_edges->setName("Metaedges");
      lod2_edges->setColor (Qt::black);

      vertices.clear();
      polygons.clear();
      fcolors.clear();
      vcolors.clear();

      Add_triangle_with_metaface_id adder2 (vertices, polygons, lod2_edges->polylines);
    
      std::tie (first_building_facet, first_wall_facet, first_vegetation_facet)
        = lod.output_lod2_to_triangle_soup
        (std::back_inserter (vertices),
         boost::make_function_output_iterator (adder2));

      std::cerr << vertices.size() << " ; " << polygons.size() << " ; " << lod2_edges->polylines.size() << std::endl;
    
      // Fill colors according to facet type
      for (std::size_t i = 0; i < first_building_facet; ++ i)
        fcolors.push_back (CGAL::Color(186, 189, 182));
      for (std::size_t i = first_building_facet; i < first_wall_facet; ++ i)
        fcolors.push_back (CGAL::Color(245, 121, 0));
      for (std::size_t i = first_wall_facet; i < first_vegetation_facet; ++ i)
        fcolors.push_back (CGAL::Color(77, 131, 186));
      for (std::size_t i = first_vegetation_facet; i < polygons.size(); ++ i)
        fcolors.push_back (CGAL::Color(138, 226, 52));
    
      lod2_faces->load (vertices, polygons, fcolors, vcolors);
      lod2_faces->setRenderingMode(Flat);
      scene->addItem(lod2_faces);
      scene->changeGroup(lod2_faces, lod2_item);
      scene->addItem(lod2_edges);
      scene->changeGroup(lod2_edges, lod2_item);
    }
    
    item->setVisible(false);

    task_timer.stop();

    std::cerr << "Done in " << task_timer.time() << " second(s)" << std::endl;
    
    QApplication::restoreOverrideCursor();

  }
}


#include "Level_of_detail_plugin.moc"
