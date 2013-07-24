#include "Scene_swept_volume_3_item_config.h"
#include <CGAL_demo/Io_plugin_interface.h>
#include <fstream>

class Io_swept_volume_3_plugin :
  public QObject,
  public Io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Io_plugin_interface);

public:
  QStringList nameFilters() const;
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);
};

QStringList Io_swept_volume_3_plugin::nameFilters() const {
  return QStringList() << "SV3 files (*.sv3)";
};

bool Io_swept_volume_3_plugin::canLoad() const {
  return true;
}


Scene_item* 
Io_swept_volume_3_plugin::load(QFileInfo fileinfo) {

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }
    
  // Try to read .sv3 in a swept_volume_3
  Scene_swept_volume_3_item* item = new Scene_swept_volume_3_item();
  item->setName(fileinfo.baseName());
  if(!item->load(in)){
    delete item;
    return 0;
  }
  return item;
}

bool Io_swept_volume_3_plugin::canSave(const Scene_item* item)
{
  // This plugin supports swept_volume_3s
  return qobject_cast<const Scene_swept_volume_3_item*>(item);
}

bool Io_swept_volume_3_plugin::save(const Scene_item* item, QFileInfo fileinfo)
{
  // This plugin supports swept_volume_3s 
  const Scene_swept_volume_3_item* sv3_item = 
    qobject_cast<const Scene_swept_volume_3_item*>(item);
  
  if(!sv3_item) return false;
  
  std::ofstream out(fileinfo.filePath().toUtf8());
  
  return (sv3_item && sv3_item->save(out));
}

#include <QtPlugin>
Q_EXPORT_PLUGIN2(Io_swept_volume_3_plugin, Io_swept_volume_3_plugin);
#include "Io_swept_volume_3_plugin.moc"
