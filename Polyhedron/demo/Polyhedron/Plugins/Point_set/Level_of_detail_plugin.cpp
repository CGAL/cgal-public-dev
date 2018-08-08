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

#include <CGAL/Level_of_detail.h>
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
  }

  double scale() const { return scaleDoubleSpinBox->value(); }
  double noise_level() const { return noiseLevelDoubleSpinBox->value(); }
  std::size_t min_points_per_wall() const { return std::size_t(minimumPointsPerWallSpinBox->value()); }
  double normal_threshold() const { return normalThresholdDoubleSpinBox->value(); }
  bool regularize() const { return regularizeSegmentsCheckBox->isChecked(); }
  double maximum_regularized_angle() const { return maximumRegularizedAngleDoubleSpinBox->value(); }
  
  bool detailed() const { return detailedOutput->isChecked(); }

  QString ground_indices() const { return groundLineEdit->text(); }
  QString boundary_indices() const { return buildingBoundariesLineEdit->text(); }
  QString interior_indices() const { return buildingInteriorLineEdit->text(); }
  QString vegetation_indices() const { return vegetationLineEdit->text(); }

  void set_ground_index (const int& i) const { return groundLineEdit->setText(tr("%1").arg(i)); }
  void set_boundary_index (const int& i) const { return buildingBoundariesLineEdit->setText(tr("%1").arg(i)); }
  void set_interior_index (const int& i) const { return buildingInteriorLineEdit->setText(tr("%1").arg(i)); }
  void set_vegetation_index (const int& i) const
  {
    return vegetationLineEdit->setText(tr("%1 %2").arg(vegetationLineEdit->text()).arg(i));
  }

  QTextEdit* comment_section() { return textEdit; }
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
            if (label == "ground")
              dialog.set_ground_index(number);
            else if (label == "facade" || label == "wall")
              dialog.set_boundary_index(number);
            else if (label == "roof" || label == "building")
              dialog.set_interior_index(number);
            else if (label == "vegetation" ||
                     label == "low_veget" ||
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

    double scale = dialog.scale();
    double noise_tolerance = dialog.noise_level();

    Scene_group_item* group = NULL;
    if (dialog.detailed())
    {
      group = new Scene_group_item(tr("%1 (LOD detailed output)").arg(item->name()));
      scene->addItem(group);
    }
      
    lod.compute_planar_ground();
				
    lod.detect_building_boundaries(scale / 2., // alpha shape size
                                   noise_tolerance, // region growing epsilon
                                   scale, // region growing cluster epsilon
                                   dialog.normal_threshold(), // region growing normal threshold
                                   dialog.min_points_per_wall(), // region growing min points
                                   0, 0, // no regularization
                                   scale / 4.); // grid cell width

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
      
    lod.partition(scale / 2.);

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

      CGAL::Level_of_detail::internal::Indexer<Kernel::Point_2> indexer;
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
      
    lod.compute_footprints(scale / 2.);

    if (dialog.detailed())
    {
      Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item;
      new_item->setName("Buildings");

      CGAL::Level_of_detail::internal::Indexer<Kernel::Point_2> indexer;
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
      
    lod.extrude_footprints();
      
    lod.compute_smooth_ground(noise_tolerance);

    Scene_polygon_soup_item* lod0_item = new Scene_polygon_soup_item;

    std::vector<Kernel::Point_3> vertices;
    std::vector<std::vector<std::size_t> > polygons;
    
    std::size_t first_building_facet
      = lod.output_lod0_to_triangle_soup
      (std::back_inserter (vertices),
       boost::make_function_output_iterator (array_to_vector(polygons)));
    
    std::vector<CGAL::Color> fcolors;
    std::vector<CGAL::Color> vcolors;

    // Fill colors according to facet type
    for (std::size_t i = 0; i < first_building_facet; ++ i)
      fcolors.push_back (CGAL::Color(186, 189, 182));
    for (std::size_t i = first_building_facet; i < polygons.size(); ++ i)
      fcolors.push_back (CGAL::Color(245, 121, 0));
    
    lod0_item->load (vertices, polygons, fcolors, vcolors);
    lod0_item->setName(tr("%1 (LOD0)").arg(item->name()));
    lod0_item->setRenderingMode(Flat);
    lod0_item->setVisible (false);
    scene->addItem(lod0_item);
    
    Scene_polygon_soup_item* lod1_item = new Scene_polygon_soup_item;

    vertices.clear();
    polygons.clear();
    fcolors.clear();
    vcolors.clear();

    std::size_t first_wall_facet;
    boost::tie (first_building_facet, first_wall_facet)
      = lod.output_lod1_to_triangle_soup
      (std::back_inserter (vertices),
       boost::make_function_output_iterator (array_to_vector(polygons)));

    // Fill colors according to facet type
    for (std::size_t i = 0; i < first_building_facet; ++ i)
      fcolors.push_back (CGAL::Color(186, 189, 182));
    for (std::size_t i = first_building_facet; i < first_wall_facet; ++ i)
      fcolors.push_back (CGAL::Color(245, 121, 0));
    for (std::size_t i = first_wall_facet; i < polygons.size(); ++ i)
      fcolors.push_back (CGAL::Color(77, 131, 186));
    
    lod1_item->load (vertices, polygons, fcolors, vcolors);
    lod1_item->setName(tr("%1 (LOD1)").arg(item->name()));
    lod1_item->setRenderingMode(Flat);
    scene->addItem(lod1_item);

    item->setVisible(false);

    task_timer.stop();

    std::cerr << "Done in " << task_timer.time() << " second(s)" << std::endl;
    
    QApplication::restoreOverrideCursor();

  }
}


#include "Level_of_detail_plugin.moc"
