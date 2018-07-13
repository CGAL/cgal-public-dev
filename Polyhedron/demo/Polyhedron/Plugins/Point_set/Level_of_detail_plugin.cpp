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

#include "ui_Level_of_detail_plugin.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QMessageBox>

using namespace CGAL::Three;

typedef CGAL::Level_of_detail::Level_of_detail<Kernel, Point_set, Point_set::Point_map> LOD;


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

  Kernel::Point_3 position_on_plane (const Kernel::Plane_3& plane, const Kernel::Point_2& point)
  {
    static Kernel::Vector_3 vertical (0., 0., 1.);
    Kernel::Line_3 line (Kernel::Point_3 (point.x(), point.y(), 0.), vertical);
    
    CGAL::cpp11::result_of<Kernel::Intersect_3(Kernel::Line_3, Kernel::Plane_3)>::type
      inter = CGAL::intersection (line, plane);
    
    if (inter)
      if (const Kernel::Point_3* p = boost::get<Kernel::Point_3>(&*inter))
        return *p;

    std::cerr << "Error: can't compute 3D position" << std::endl;
    return Kernel::Point_3 (0., 0., 0.);
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

struct Visibility_map_from_labels
{
  typedef Point_set::Index key_type;
  typedef double value_type;
  typedef double reference;
  typedef boost::readable_property_map_tag category;

  Point_set* points;
  Point_set::Property_map<int> label_map;
  Map_l2sl_ptr map_l2sl;

  Visibility_map_from_labels (Point_set* points, Map_l2sl_ptr map_l2sl)
    : points (points), map_l2sl (map_l2sl)
  {
    label_map = points->property_map<int>("label").first;
  }

  friend value_type get (const Visibility_map_from_labels& map, const key_type& key)
  {
    int l = map.label_map[key];

    typename Map_l2sl::const_iterator
      found = map.map_l2sl->find(l);
    if (found == map.map_l2sl->end())
      return 0.;
    if (found->second == CGAL::Level_of_detail::Semantic_label::BUILDING_INTERIOR)
      return 1.;
    if (found->second == CGAL::Level_of_detail::Semantic_label::BUILDING_BOUNDARY)
      return 0.5;

    return 0.; // ground, unassigned, vegetation
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

    LOD::Parameters parameters;
    parameters.verbose() = true;
    parameters.scale() = dialog.scale();
    parameters.epsilon() = dialog.noise_level();
    parameters.update_dependent();
    
    parameters.alpha_shape_size() = dialog.scale() / 2.;
    parameters.grid_cell_width() = dialog.scale() / 4.;
    
    parameters.region_growing_2_epsilon() = dialog.noise_level();
    parameters.region_growing_2_cluster_epsilon() = dialog.scale();
    parameters.region_growing_2_min_points() = dialog.min_points_per_wall();
    parameters.region_growing_2_normal_threshold() = dialog.normal_threshold();

    parameters.no_regularization() = !dialog.regularize();
    parameters.segment_regularizer_2_max_difference_in_meters() = dialog.noise_level();
    parameters.kinetic_partitioning_2_min_face_width() = dialog.scale() / 2.;
    parameters.segment_regularizer_2_max_angle_in_degrees() = dialog.maximum_regularized_angle();
    
    LOD lod (*points, points->point_map(), parameters);
    Semantic_map_from_labels semantic_map (points);
    Visibility_map_from_labels visibility_map (points, semantic_map.map_l2sl);

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

    if (dialog.detailed())
    {
      Scene_group_item* group = new Scene_group_item(tr("%1 (LOD detailed output)").arg(item->name()));
      scene->addItem(group);
      
      const typename LOD::Data_structure& data = lod.get_internal_data_structure();
      
      std::cerr << "Building LOD with detailed output" << std::endl;
      
      lod.split_semantic_data(semantic_map);
				
      lod.fit_ground_plane();

      const Kernel::Plane_3& ground_plane = data.ground_plane();
				
      lod.extract_building_boundaries();

      {
        Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item;
        new_item->setName("Boundary points");

        for (std::size_t i = 0; i < data.filtered_building_boundary_points().size(); ++ i)
          new_item->point_set()->insert (points->point (*(data.filtered_building_boundary_points()[i])));
        new_item->setVisible(false);
        new_item->setColor (Qt::red);
        scene->addItem(new_item);
        scene->changeGroup(new_item, group);
      }
      
      lod.simplify_building_boundaries();

      {
        Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item;
        new_item->setName("Simplified boundary points");

        for (std::size_t i = 0; i < data.simplified_building_boundary_points().size(); ++ i)
          new_item->point_set()->insert (points->point (*(data.simplified_building_boundary_points()[i])));
        new_item->setVisible(false);
        new_item->setColor (Qt::red);
        scene->addItem(new_item);
        scene->changeGroup(new_item, group);
      }
      
      lod.detect_lines();

      {
        Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item;
        new_item->setName("Detected 2D regions (points)");

        Point_set::Property_map<unsigned char>
          red = new_item->point_set()->template add_property_map<unsigned char>("r", 128).first;
        Point_set::Property_map<unsigned char>
          green = new_item->point_set()->template add_property_map<unsigned char>("g", 128).first;
        Point_set::Property_map<unsigned char>
          blue = new_item->point_set()->template add_property_map<unsigned char>("b", 128).first;
        new_item->point_set()->check_colors();

        for (typename LOD::Data_structure::Detected_regions::const_iterator
                 it = data.detected_2d_regions().begin();
               it !=  data.detected_2d_regions().end(); ++ it)
        {
          unsigned char r, g, b;

          r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
          g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
          b = static_cast<unsigned char>(64 + rand.get_int(0, 192));

          for (std::size_t i = 0; i < it->size(); ++ i)
          {
            Point_set::iterator pt
              = new_item->point_set()->insert (points->point (*((*it)[i])));
            red[*pt] = r; green[*pt] = g; blue[*pt] = b;
          }
        }
        new_item->setVisible(false);
        new_item->setColor (Qt::black);
        scene->addItem(new_item);
        scene->changeGroup(new_item, group);
      }

      {
        Scene_polylines_item* new_item = new Scene_polylines_item;
        new_item->setName("Detected 2D regions (segments)");

        const CGAL::Level_of_detail::Segment_from_region_property_map_2
          <Kernel,
           typename LOD::Data_structure::Point_identifiers,
           typename LOD::Point_map_2,
           typename LOD::Data_structure::Regularized_segments>
          segment_from_region_map_2(data.detected_2d_regions(), lod.point_map_2());
        
        for (typename LOD::Data_structure::Detected_regions::const_iterator
                 it = data.detected_2d_regions().begin();
               it !=  data.detected_2d_regions().end(); ++ it)
        {
          const Kernel::Segment_2 &segment = get (segment_from_region_map_2, *it);
          new_item->polylines.push_back (Scene_polylines_item::Polyline());
          new_item->polylines.back().push_back (position_on_plane(ground_plane, segment.source()));
          new_item->polylines.back().push_back (position_on_plane(ground_plane, segment.target()));
        }
        
        new_item->setVisible(false);
        new_item->setColor (Qt::black);
        scene->addItem(new_item);
        scene->changeGroup(new_item, group);
      }
      
      lod.regularize_segments();

      if (dialog.regularize())
      {
        Scene_polylines_item* new_item = new Scene_polylines_item;
        new_item->setName("Regularized segments");

        for (typename LOD::Data_structure::Regularized_segments::const_iterator
                 it = data.regularized_segments().begin();
               it !=  data.regularized_segments().end(); ++ it)
        {
          new_item->polylines.push_back (Scene_polylines_item::Polyline());
          new_item->polylines.back().push_back (position_on_plane(ground_plane, it->source()));
          new_item->polylines.back().push_back (position_on_plane(ground_plane, it->target()));
        }
        
        new_item->setVisible(false);
        new_item->setColor (Qt::black);
        scene->addItem(new_item);
        scene->changeGroup(new_item, group);
      }
      
      lod.create_partitioning();

      {
        Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item;
        new_item->setName("Partitioning");

        CGAL::Level_of_detail::internal::Indexer<Kernel::Point_2> indexer;
        std::vector<Kernel::Point_3> vertices;
        std::vector<std::vector<std::size_t> > polygons;
        std::vector<CGAL::Color> fcolors;
        std::vector<CGAL::Color> vcolors;

        for (typename LOD::Data_structure::Partition_faces_2::const_iterator
               it = data.partition_faces_2().begin();
             it != data.partition_faces_2().end(); ++ it)
        {
          unsigned char r, g, b;
          r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
          g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
          b = static_cast<unsigned char>(64 + rand.get_int(0, 192));

          polygons.push_back (std::vector<std::size_t>());
          fcolors.push_back (CGAL::Color(r,g,b));
          
          for (typename LOD::Data_structure::Partition_face_2::const_iterator
                 fit = it->begin(); fit != it->end(); ++ fit)
          {
            std::size_t idx = indexer(*fit);
            if (idx == vertices.size())
              vertices.push_back (position_on_plane (ground_plane, *fit));
            polygons.back().push_back (idx);
          }
        }

        new_item->load (vertices, polygons, fcolors, vcolors);
        new_item->setVisible(false);
        scene->addItem(new_item);
        scene->changeGroup(new_item, group);

        std::ofstream out("debug.ply", std::ios::binary);
        out.precision (std::numeric_limits<double>::digits10 + 2);
        CGAL::write_PLY (out, new_item->points(), new_item->polygons());
      }
      
      lod.compute_visibility(visibility_map);

      {
        Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item;
        new_item->setName("Visibility");

        CGAL::Level_of_detail::internal::Indexer<Kernel::Point_2> indexer;
        std::vector<Kernel::Point_3> vertices;
        std::vector<std::vector<std::size_t> > polygons;
        std::vector<CGAL::Color> fcolors;
        std::vector<CGAL::Color> vcolors;

        for (typename LOD::Data_structure::Partition_faces_2::const_iterator
               it = data.partition_faces_2().begin();
             it != data.partition_faces_2().end(); ++ it)
        {
          const CGAL::Level_of_detail::Visibility_label visibility_label
            = it->visibility_label();
          
          unsigned char r, g, b;
          if (visibility_label == CGAL::Level_of_detail::Visibility_label::OUTSIDE)
          {
            r = 186; g = 189; b = 182;
          }
          else if (visibility_label == CGAL::Level_of_detail::Visibility_label::INSIDE)
          {
            r = 245; g = 121; b = 0;
          }
          else
          {
            r = 77; g = 131; b = 186;
          }

          polygons.push_back (std::vector<std::size_t>());
          fcolors.push_back (CGAL::Color(r,g,b));
          
          for (typename LOD::Data_structure::Partition_face_2::const_iterator
                 fit = it->begin(); fit != it->end(); ++ fit)
          {
            std::size_t idx = indexer(*fit);
            if (idx == vertices.size())
              vertices.push_back (position_on_plane (ground_plane, *fit));
            polygons.back().push_back (idx);
          }
        }

        new_item->load (vertices, polygons, fcolors, vcolors);
        new_item->setVisible(false);
        scene->addItem(new_item);
        scene->changeGroup(new_item, group);
      }
      
      lod.create_triangulation();

      lod.find_buildings();

      {
        Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item;
        new_item->setName("Buildings");

        CGAL::Level_of_detail::internal::Indexer<Kernel::Point_2> indexer;
        std::vector<Kernel::Point_3> vertices;
        std::vector<std::vector<std::size_t> > polygons;
        std::vector<CGAL::Color> fcolors;
        std::vector<CGAL::Color> vcolors;

        for (typename LOD::Triangulation::Finite_faces_iterator
               it = data.triangulation().finite_faces_begin();
             it != data.triangulation().finite_faces_end(); ++ it)
        {
          int building_id = it->info().group_number();

          CGAL::Random rand (building_id);
          unsigned char r, g, b;
          if (building_id < 0)
          {
            r = 128; g = 128; b = 128;
          }
          else
          {
            r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
            g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
            b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
          }

          polygons.push_back (std::vector<std::size_t>());
          fcolors.push_back (CGAL::Color(r,g,b));
          
          for (std::size_t i = 0; i < 3; ++ i)
          {
            std::size_t idx = indexer(it->vertex(i)->point());
            if (idx == vertices.size())
              vertices.push_back (position_on_plane (ground_plane, it->vertex(i)->point()));
            polygons.back().push_back (idx);
          }
        }
        new_item->load (vertices, polygons, fcolors, vcolors);
        new_item->setRenderingMode(Flat);
        new_item->setVisible(false);
        scene->addItem(new_item);
        scene->changeGroup(new_item, group);
      }
      
      lod.find_building_walls();

      {
        Scene_polylines_item* new_item = new Scene_polylines_item;
        new_item->setName("Building walls");

        for (typename LOD::Data_structure::Buildings::const_iterator
                 it = data.buildings().begin();
               it !=  data.buildings().end(); ++ it)
        {
          for (typename LOD::Data_structure::Building::Floor_edges::const_iterator
                 bit = it->floor_edges().begin();
               bit != it->floor_edges().end(); ++ bit)
          {
            new_item->polylines.push_back (Scene_polylines_item::Polyline());
            new_item->polylines.back().push_back (position_on_plane(ground_plane, bit->source()));
            new_item->polylines.back().push_back (position_on_plane(ground_plane, bit->target()));
          }
        }
        
        new_item->setVisible(false);
        new_item->setColor (Qt::black);
        scene->addItem(new_item);
        scene->changeGroup(new_item, group);
      }
      
      lod.fit_flat_building_roofs();
      
      lod.compute_triangulation_vertices_heights();
    }
    else
      lod.build (semantic_map, visibility_map);

    Scene_polygon_soup_item* lod0_item = new Scene_polygon_soup_item;

    std::vector<Kernel::Point_3> vertices;
    std::vector<std::vector<std::size_t> > polygons;
    
    std::size_t first_building_facet
      = lod.output_lod0_to_polygon_soup (vertices, polygons);
    
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
      = lod.output_lod1_to_polygon_soup (vertices, polygons);

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
