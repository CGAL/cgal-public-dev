#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include "ui_Mesh_segmentation_widget.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include "Scene.h"
#include "Color_map.h"

#include <CGAL/mesh_segmentation.h>
#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QTime>
#include <QAction>
#include <QDebug>
#include <QObject>
#include <QDockWidget>
//#include <QtConcurrentRun>
#include <map>
#include <algorithm>
#include <vector>
#include <CGAL/property_map.h>


template<class PolyhedronWithId, class ValueType>
struct Polyhedron_with_id_to_vector_property_map
    : public boost::put_get_helper<ValueType&,
             Polyhedron_with_id_to_vector_property_map<PolyhedronWithId, ValueType> >
{
public:
    typedef typename PolyhedronWithId::Facet_const_handle key_type;
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;

    Polyhedron_with_id_to_vector_property_map() : internal_vector(NULL) { }
    Polyhedron_with_id_to_vector_property_map(std::vector<ValueType>* internal_vector)
         : internal_vector(internal_vector) { }
        
    reference operator[](key_type key) const { return (*internal_vector)[key->id()]; }
private:
    std::vector<ValueType>* internal_vector;
};
    
class Polyhedron_demo_mesh_segmentation_plugin : 
    public QObject,
    public Polyhedron_demo_plugin_helper
{
    Q_OBJECT
        Q_INTERFACES(Polyhedron_demo_plugin_interface)
private:
    typedef std::map<Scene_polyhedron_item*, std::vector<double> > Item_sdf_map;
public:

    QList<QAction*> actions() const {
        return QList<QAction*>() << actionSegmentation;
    }

    bool applicable(QAction*) const {
      return 
        qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
    }    
    
    void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
        this->scene = scene_interface;
        this->mw = mainWindow;
        actionSegmentation = new QAction("Mesh Segmentation", mw);
        connect(actionSegmentation, SIGNAL(triggered()),this, SLOT(on_actionSegmentation_triggered()));

        // adding slot for itemAboutToBeDestroyed signal, aim is removing item from item-functor map.
        
        if( Scene* scene = dynamic_cast<Scene*>(scene_interface) ) {
            connect(scene, SIGNAL(itemAboutToBeDestroyed(Scene_item*)), this, SLOT(itemAboutToBeDestroyed(Scene_item*)));
        }
        
        init_color_map_sdf();
        init_color_map_segmentation();

        dock_widget = new QDockWidget("Mesh segmentation parameters", mw);
        dock_widget->setVisible(false); // do not show at the beginning
        ui_widget.setupUi(dock_widget);
        mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);
    
        connect(ui_widget.Partition_button,  SIGNAL(clicked()), this, SLOT(on_Partition_button_clicked()));   
        connect(ui_widget.SDF_button,  SIGNAL(clicked()), this, SLOT(on_SDF_button_clicked()));   
    }
    
    template<class SDFPropertyMap>
    void colorize_sdf(Scene_polyhedron_item* item, SDFPropertyMap sdf_values, std::vector<QColor>& color_vector);
    template<class SegmentPropertyMap> 
    void colorize_segmentation(Scene_polyhedron_item* item, SegmentPropertyMap segment_ids, std::vector<QColor>& color_vector);
    void check_and_set_ids(Polyhedron* polyhedron);
    void init_color_map_sdf();
    void init_color_map_segmentation();
    
    public Q_SLOTS:
        void on_actionSegmentation_triggered();
        void on_Partition_button_clicked();
        void on_SDF_button_clicked();
        void itemAboutToBeDestroyed(Scene_item*);
private:
    QAction*                      actionSegmentation;
    QDockWidget*                  dock_widget;
    Ui::Mesh_segmentation         ui_widget;
    
    std::vector<QColor>  color_map_sdf;
    std::vector<QColor>  color_map_segmentation;
    Item_sdf_map         item_sdf_map;
};

void Polyhedron_demo_mesh_segmentation_plugin::init_color_map_sdf()
{
    color_map_sdf = std::vector<QColor>(256);
    int r = 0, g = 0, b = 255;
    for(int i = 0; i <= 255; ++i)
    {
        if(i > 128 && i <= 192) { r = static_cast<int>( ((i - 128) / (192.0 - 128)) * 255 ); }        
        if(i > 0 && i <= 98)    { g = static_cast<int>( ((i) / (98.0)) * 255 ); }
        if(i > 191 && i <=255)  { g = 255 - static_cast<int>( ((i - 191) / (255.0 - 191)) * 255 ); }
        if(i > 64 && i <= 127)  { b = 255 - static_cast<int>( ((i - 64) / (127.0 - 64)) * 255 ); }
        color_map_sdf[i] = QColor(r, g, b);        
    }
}

void Polyhedron_demo_mesh_segmentation_plugin::init_color_map_segmentation()
{
    /*
    color_map_segmentation.push_back(QColor("#fce94f"));
    color_map_segmentation.push_back(QColor("#edd400"));
    color_map_segmentation.push_back(QColor("#c4a000"));
    color_map_segmentation.push_back(QColor("#fcaf3e"));
    color_map_segmentation.push_back(QColor("#f57900"));
    color_map_segmentation.push_back(QColor("#ce5c00")); 
    color_map_segmentation.push_back(QColor("#e9b96e"));
    color_map_segmentation.push_back(QColor("#c17d11"));
    color_map_segmentation.push_back(QColor("#8f5902")); 
    color_map_segmentation.push_back(QColor("#729fcf"));
    color_map_segmentation.push_back(QColor("#3465a4"));
    color_map_segmentation.push_back(QColor("#204a87"));
    color_map_segmentation.push_back(QColor("#ad7fa8"));
    color_map_segmentation.push_back(QColor("#75507b"));
    color_map_segmentation.push_back(QColor("#5c3566"));
    */
    
    color_map_segmentation.push_back(QColor( 173, 35, 35)); 
    color_map_segmentation.push_back(QColor( 87, 87, 87));    
    color_map_segmentation.push_back(QColor( 42, 75, 215)); 
    color_map_segmentation.push_back(QColor( 29, 105, 20)); 
    color_map_segmentation.push_back(QColor( 129, 74, 25)); 
    color_map_segmentation.push_back(QColor( 129, 38, 192)); 
    color_map_segmentation.push_back(QColor( 160, 160, 160)); 
    color_map_segmentation.push_back(QColor( 129, 197, 122)); 
    color_map_segmentation.push_back(QColor( 157, 175, 255)); 
    color_map_segmentation.push_back(QColor( 41, 208, 208)); 
    color_map_segmentation.push_back(QColor( 255, 146, 51)); 
    color_map_segmentation.push_back(QColor( 255, 238, 51)); 
    color_map_segmentation.push_back(QColor( 233, 222, 187)); 
    color_map_segmentation.push_back(QColor( 255, 205, 243)); 
    
}

void Polyhedron_demo_mesh_segmentation_plugin::itemAboutToBeDestroyed(Scene_item* scene_item)
{
    if(Scene_polyhedron_item* item = qobject_cast<Scene_polyhedron_item*>(scene_item)) {
      item_sdf_map.erase(item);
    }
}

void Polyhedron_demo_mesh_segmentation_plugin::on_actionSegmentation_triggered()
{ dock_widget->show(); }

void Polyhedron_demo_mesh_segmentation_plugin::on_SDF_button_clicked()
{
    Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_polyhedron_item* item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(!item) { return; }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    
    std::size_t number_of_rays = ui_widget.Number_of_rays_spin_box->value();
    double cone_angle = (ui_widget.Cone_angle_spin_box->value()  / 180.0) * CGAL_PI;
    bool create_new_item = ui_widget.New_item_check_box->isChecked();
    
    Item_sdf_map::iterator pair;
    Scene_polyhedron_item* active_item = item;

    if(create_new_item) {
        active_item = new Scene_polyhedron_item(*item->polyhedron()); 
        active_item->setGouraudMode();
    }
    
    pair = item_sdf_map.insert(
            std::make_pair(active_item, std::vector<double>()) ).first; 
    
    check_and_set_ids(pair->first->polyhedron());
    pair->second.resize(item->polyhedron()->size_of_facets(), 0.0);
    Polyhedron_with_id_to_vector_property_map<Polyhedron, double>  sdf_pmap(&pair->second);
    QTime time;
    time.start();
    std::pair<double, double> min_max_sdf = sdf_values(*(pair->first->polyhedron()), sdf_pmap, cone_angle, number_of_rays);
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

    std::cout << "SDF computation is completed. Min-SDF : " << min_max_sdf.first << " " "Max-SDF : " << min_max_sdf.second << std::endl;

    pair->first->set_color_vector_read_only(true);
    colorize_sdf(pair->first, sdf_pmap, pair->first->color_vector());
       
    pair->first->setName(tr("(SDF-%1-%2)").arg(number_of_rays).arg(ui_widget.Cone_angle_spin_box->value()));
    
    if(create_new_item) {
        index = scene->addItem(pair->first);
        item->setVisible(false);
        scene->setSelectedItem(index);
    }
    else {
      scene->itemChanged(index);
    }

    QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mesh_segmentation_plugin::on_Partition_button_clicked()
{    
    Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_polyhedron_item* item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(!item) { return; }
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    
    std::size_t number_of_clusters = ui_widget.Number_of_clusters_spin_box->value();
    double smoothness = ui_widget.Smoothness_spin_box->value();
    std::size_t number_of_rays = ui_widget.Number_of_rays_spin_box->value();
    double cone_angle = (ui_widget.Cone_angle_spin_box->value()  / 180.0) * CGAL_PI;
    bool create_new_item = ui_widget.New_item_check_box->isChecked();
    bool extract_segments = ui_widget.Extract_segments_check_box->isChecked();

    Item_sdf_map::iterator pair;
    if(create_new_item)
    {
        // create new item
        Scene_polyhedron_item* new_item = new Scene_polyhedron_item(*item->polyhedron()); 
        new_item->setGouraudMode(); 
        
        // copy SDF values of existing poly to new poly
        Item_sdf_map::iterator it = item_sdf_map.find(item);
        const std::vector<double>& sdf_data = it == item_sdf_map.end() ?
                                              std::vector<double>() : it->second;
        pair = item_sdf_map.insert(std::make_pair(new_item, sdf_data) ).first;
    }
    else
    {
        std::pair<Item_sdf_map::iterator, bool> res = 
          item_sdf_map.insert(std::make_pair(item, std::vector<double>()) );
        pair = res.first;
    }

    check_and_set_ids(pair->first->polyhedron());
    QTime time;
    time.start();
    if(pair->second.empty()) { // SDF values are empty, calculate
      pair->second.resize(pair->first->polyhedron()->size_of_facets(), 0.0);
      Polyhedron_with_id_to_vector_property_map<Polyhedron, double> sdf_pmap(&pair->second);
      sdf_values(*(pair->first->polyhedron()), sdf_pmap, cone_angle, number_of_rays); 
    }

    std::vector<std::size_t> internal_segment_map(pair->first->polyhedron()->size_of_facets());
    Polyhedron_with_id_to_vector_property_map<Polyhedron, std::size_t> segment_pmap(&internal_segment_map);
    Polyhedron_with_id_to_vector_property_map<Polyhedron, double> sdf_pmap(&pair->second);

    std::size_t nb_segments = segmentation_from_sdf_values(*(pair->first->polyhedron())
        ,sdf_pmap, segment_pmap, number_of_clusters, smoothness, extract_segments); 
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
    std::cout << "Segmentation is completed. Number of segments : " << nb_segments << std::endl;  
    pair->first->set_color_vector_read_only(true);  
     
    colorize_segmentation(pair->first, segment_pmap, pair->first->color_vector());
    pair->first->setName(tr("(Segmentation-%1-%2)").arg(number_of_clusters).arg(smoothness));   

    if(create_new_item) {
        index = scene->addItem(pair->first);
        item->setVisible(false);
        scene->setSelectedItem(index);
    }
    else {
      scene->itemChanged(index);
    }

    QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mesh_segmentation_plugin::check_and_set_ids(Polyhedron* polyhedron)
{
    Polyhedron::Facet_iterator a_facet = polyhedron->facets_begin();
    Polyhedron::Facet_iterator another_facet = ++polyhedron->facets_begin();
    if(a_facet->id() != another_facet->id()) { return; } // ids are OK
    std::size_t facet_id = 0;
    for(Polyhedron::Facet_iterator facet_it = polyhedron->facets_begin();
        facet_it != polyhedron->facets_end(); ++facet_it, ++facet_id)
    {
        facet_it->id() = facet_id;
    }
}

template<class SDFPropertyMap>
void Polyhedron_demo_mesh_segmentation_plugin::colorize_sdf(
     Scene_polyhedron_item* item,
     SDFPropertyMap sdf_values,  
     std::vector<QColor>& color_vector)
{
    Polyhedron* polyhedron = item->polyhedron();
    color_vector.clear();
    std::size_t patch_id = 0;
    for(Polyhedron::Facet_iterator facet_it = polyhedron->facets_begin(); 
        facet_it != polyhedron->facets_end(); ++facet_it, ++patch_id)   
    {
        double sdf_value = sdf_values[facet_it]; 
        int gray_color = static_cast<int>(255 * sdf_value);
        if(gray_color < 0 || gray_color >= 256) {
          color_vector.push_back(QColor::fromRgb(0,0,0));
        }
        else {
          color_vector.push_back(color_map_sdf[gray_color]);
        }
        facet_it->set_patch_id(static_cast<int>(patch_id));
    }
}

template<class SegmentPropertyMap>
void Polyhedron_demo_mesh_segmentation_plugin::colorize_segmentation(
     Scene_polyhedron_item* item,
     SegmentPropertyMap segment_ids,
     std::vector<QColor>& color_vector)
{
    Polyhedron* polyhedron = item->polyhedron();
    color_vector.clear();
    std::size_t max_segment = 0;
    for(Polyhedron::Facet_iterator facet_it = polyhedron->facets_begin(); 
        facet_it != polyhedron->facets_end(); ++facet_it)   
    {
        std::size_t segment_id = segment_ids[facet_it];
        facet_it->set_patch_id(static_cast<int>(segment_id));
        max_segment = (std::max)(max_segment, segment_id);      
    }
    for(std::size_t i = 0; i <= max_segment; ++i)   
    {
        QColor aColor = color_map_segmentation[(max_segment - i) % color_map_segmentation.size()]; 
        color_vector.push_back(aColor);     
    }    
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_mesh_segmentation_plugin, Polyhedron_demo_mesh_segmentation_plugin)

#include "Polyhedron_demo_mesh_segmentation_plugin.moc"
