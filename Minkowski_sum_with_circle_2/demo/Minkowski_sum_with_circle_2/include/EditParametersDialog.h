#ifndef KGON_PARAMETERS_2_MAIN_WINDOW_H
#define KGON_PARAMETERS_2_MAIN_WINDOW_H

// Qt headers
#include <QtGui>
#include <QDialog>

// for viewportsBbox
#include <CGAL/Qt/utility.h>
// the two base classes
#include "ui_Parameters.h"

class EditParametersDialog :
  public QDialog, 
  private Ui::Parameters
{
  Q_OBJECT
  
  public:
  
  struct Parameters
  {
    int type;
    int size;
    
    std::string radius;
    std::string epsilon;
    std::string delta;
    
    std::string radiusRatio;
    std::string epsilonRatio;
    std::string deltaRatio;
    
    bool radiusKeepRatio;
    bool epsilonKeepRatio;
    bool deltaKeepRatio;

    int sizeP;
    std::string bboxP;

    Parameters(): type(0), size(8),
                  radius("1/1"), epsilon("1/10"), delta("1/40"),
                  radiusRatio("1/1"), epsilonRatio("1/10"), deltaRatio("1/4"),
                  radiusKeepRatio(false), epsilonKeepRatio(false), deltaKeepRatio(false),
                  sizeP(0), bboxP("1/1")
                  
    {
    }
    
    bool equal(const Parameters& p) { 
      return (type == p.type) && (size == p.size) && 
        (radius == p.radius) && (epsilon == p.epsilon) && (delta == p.delta); 
    }
    
  };
  
  EditParametersDialog(QWidget *parent = 0);

  static int saved_type() {
    return m_saved.type;
  }

  static int saved_size() {
    return m_saved.size;
  }

  static void saved_size(const int& saved_size) {
    m_saved.size = saved_size;
  }

  static std::string saved_radius() {
    return m_saved.radius;
  }

  static std::string saved_epsilon() {
    return m_saved.epsilon;
  }

  static void saved_epsilon(const std::string& saved_epsilon) {
    m_saved.epsilon = saved_epsilon;
  }

  static std::string saved_delta() {
    return m_saved.delta;
  }

  static void saved_delta(const std::string& saved_delta) {
    m_saved.delta = saved_delta;
  }

  static void saved_polygon_size(const int& size) {
    m_saved.sizeP = size;
  }
  static void saved_polygon_bbox(const std::string& bbox_size) {
    m_saved.bboxP = bbox_size;
  }
  
signals:
  
  void changedParameters(const Parameters& parameters);
                                                      
public slots:
  
  void accept();
  void reject();
  
  void editType(int type);
  void editSize(int size);
  
  void editRadius(QString radius);
  void editEpsilon(QString epsilon);
  void editDelta(QString delta);
  
  void editRadiusRatio(QString radiusRatio);
  void editEpsilonRatio(QString epsilonRatio);
  void editDeltaRatio(QString deltaRatio);

  void changedRadius();
  void changedEpsilon();
  void changedDelta();

  void changedRadiusRatio();
  void changedEpsilonRatio();
  void changedDeltaRatio();
  
  void keepRadiusRatio(bool radiusRatio);
  void keepEpsilonRatio(bool epsilonRatio);
  void keepDeltaRatio(bool deltaRatio);

public:
private:

  struct ValueAndRatio
  {
    std::string& value;
    std::string& ratio;
    QLineEdit* valueParam;
    QLineEdit* ratioParam;
    const bool editParam;

    ValueAndRatio(std::string& i_value, std::string& i_ratio,
      QLineEdit* i_valueParam, QLineEdit* i_ratioParam):
      value(i_value),
      ratio(i_ratio),
      valueParam(i_valueParam),
      ratioParam(i_ratioParam),
      editParam(true)
      {}
      
    ValueAndRatio(std::string& i_value, std::string& i_ratio):
      value(i_value),
      ratio(i_ratio),
      editParam(false)
      {}
      

    void update(const std::string& input, const bool& keepRatio)
    {
      if(keepRatio)
      {
        updateValueFromRatio(ratio, input, value);
        if(editParam) updateParamFromValue(value, valueParam);
      }
      else
      {
        updateRatioFromValue(value, input, ratio);
        if(editParam) updateParamFromValue(ratio, ratioParam);
      }
    }
  };

  static void updateValueFromRatio(const std::string& mult, const std::string& val, std::string& res);
  static void updateRatioFromValue(const std::string& num, const std::string& den, std::string& res);
  static void updateParamFromValue(const std::string& val, QLineEdit* param);

  static Parameters m_saved;  
  Parameters m_edited;
};

#endif // KGON_PARAMETERS_2_MAIN_WINDOW_H
