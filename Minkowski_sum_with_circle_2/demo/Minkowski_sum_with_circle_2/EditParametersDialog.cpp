// Edit Dialog
#include <QtGui> 
#include "include/EditParametersDialog.h"
#include "include/typedefs.h"

EditParametersDialog::Parameters EditParametersDialog::m_saved;

EditParametersDialog::EditParametersDialog(QWidget *parent)
{
  setupUi(this); // this sets up GUI

  // make data consistent (in case bbox has changed update radius ratio)
  ValueAndRatio radius(m_saved.radius, m_saved.radiusRatio);
  // can't update radius at this point, so updating ratio
  // even if it is against current "keep ratio" setting
  radius.update(m_saved.bboxP, false);
  
  // update default setup values according to our data
  m_edited = m_saved;

  // set combo and spin boxes values
  comboBoxType->setCurrentIndex(m_edited.type);
  spinBoxSize->setValue(m_edited.size);
  spinBoxSizeP->setValue(m_edited.sizeP);

  // set check boxes and enable text editing accordingly
  checkBoxRadiusRatio->setChecked(m_edited.radiusKeepRatio);
  checkBoxEpsilonRatio->setChecked(m_edited.epsilonKeepRatio);
  checkBoxDeltaRatio->setChecked(m_edited.deltaKeepRatio);

  lineEditRadius->setEnabled(!m_edited.radiusKeepRatio);
  lineEditRadiusRatio->setEnabled(m_edited.radiusKeepRatio);
  
  lineEditEpsilon->setEnabled(!m_edited.epsilonKeepRatio);
  lineEditEpsilonRatio->setEnabled(m_edited.epsilonKeepRatio);
  
  lineEditDelta->setEnabled(!m_edited.deltaKeepRatio);
  lineEditDeltaRatio->setEnabled(m_edited.deltaKeepRatio);

  // set text values
  updateParamFromValue(m_edited.bboxP, lineEditBbox);

  updateParamFromValue(m_edited.radius, lineEditRadius);
  updateParamFromValue(m_edited.epsilon, lineEditEpsilon);
  updateParamFromValue(m_edited.delta, lineEditDelta);
  
  updateParamFromValue(m_edited.radiusRatio, lineEditRadiusRatio);
  updateParamFromValue(m_edited.epsilonRatio, lineEditEpsilonRatio);
  updateParamFromValue(m_edited.deltaRatio, lineEditDeltaRatio);

  connect( this, SIGNAL( changedParameters(const Parameters&) ), 
           parent, SLOT( changedParameters(const Parameters&) ) );
           
}

void EditParametersDialog::accept()
{
//std::clog << "ok " << std::endl;
  bool changed = !(m_edited.equal(m_saved));
  
  m_saved = m_edited;
  QDialog::accept();
  
  if (changed) {
    emit changedParameters(m_saved);
  }
}

void EditParametersDialog::reject()
{
//std::clog << "cancel " << std::endl;
  m_edited = m_saved;
  QDialog::reject();
}

void EditParametersDialog::editType(int type)
{
//std::clog << "type " << type << std::endl;
  m_edited.type = type;
}

void EditParametersDialog::editSize(int size)
{
//std::clog << "size " << size << std::endl;
  m_edited.size = size;
  
  // TODO: update epsilon accordingly
}

void EditParametersDialog::editRadius(QString radius)
{
  std::string strRadius = radius.toStdString();
//std::clog << "radius " << strRadius << std::endl;
  m_edited.radius = strRadius;
}

void EditParametersDialog::editEpsilon(QString epsilon)
{
  std::string strEpsilon = epsilon.toStdString();
//std::clog << "epsilon " << strEpsilon << std::endl;
  m_edited.epsilon = strEpsilon;
}

void EditParametersDialog::editDelta(QString delta)
{
  std::string strDelta = delta.toStdString();
//std::clog << "delta " << strDelta << std::endl;
  m_edited.delta = strDelta;
}

void EditParametersDialog::editRadiusRatio(QString radiusRatio)
{
  std::string strRadiusRatio = radiusRatio.toStdString();
//std::clog << "radiusRatio " << strRadiusRatio << std::endl;
  m_edited.radiusRatio = strRadiusRatio;
}

void EditParametersDialog::editEpsilonRatio(QString epsilonRatio)
{
  std::string strEpsilonRatio = epsilonRatio.toStdString();
//std::clog << "epsilonRatio " << strEpsilonRatio << std::endl;
  m_edited.epsilonRatio = strEpsilonRatio;
}

void EditParametersDialog::editDeltaRatio(QString deltaRatio)
{
  std::string strDeltaRatio = deltaRatio.toStdString();
//std::clog << "deltaRatio " << strDeltaRatio << std::endl;
  m_edited.deltaRatio = strDeltaRatio;
}

void EditParametersDialog::updateValueFromRatio(
  const std::string& mult, 
  const std::string& val, 
  std::string& res)
{
  Minsum_with_circle_2::Input_rational rat_mult(mult);
  Minsum_with_circle_2::Input_rational rat_val(val);

  Minsum_with_circle_2::Input_rational result = rat_mult * rat_val;
  std::ostringstream rstream;
  rstream << result;
  res = rstream.str();  
}


void EditParametersDialog::updateRatioFromValue(
  const std::string& num, 
  const std::string& den,
  std::string& res)
{
  Minsum_with_circle_2::Input_rational rat_num(num);
  Minsum_with_circle_2::Input_rational rat_den(den);

  Minsum_with_circle_2::Input_rational ratio = rat_num/rat_den;
  std::ostringstream rstream;
  rstream << ratio;
  res = rstream.str();
}

void EditParametersDialog::updateParamFromValue(
  const std::string& val, 
  QLineEdit* param)
{
  param->setText(QApplication::translate("Parameters", val.c_str(), 
    0, QApplication::UnicodeUTF8));
}


void EditParametersDialog::changedRadius()
{
  // update radius ratio accordingly
  updateRatioFromValue(m_edited.radius, m_edited.bboxP, m_edited.radiusRatio);
  updateParamFromValue(m_edited.radiusRatio, lineEditRadiusRatio);

  // update epsilon or its ratio accordingly
  ValueAndRatio epsilon(m_edited.epsilon, m_edited.epsilonRatio, lineEditEpsilon, lineEditEpsilonRatio);
  epsilon.update(m_edited.radius, m_edited.epsilonKeepRatio);
  
  // update delta or its ratio accordingly (if needed)
  ValueAndRatio delta(m_edited.delta, m_edited.deltaRatio, lineEditDelta, lineEditDeltaRatio);
  delta.update(m_edited.epsilon, m_edited.deltaKeepRatio);
}

void EditParametersDialog::changedEpsilon()
{
  // update epsilon ratio accordingly
  updateRatioFromValue(m_edited.epsilon, m_edited.radius, m_edited.epsilonRatio);
  updateParamFromValue(m_edited.epsilonRatio, lineEditEpsilonRatio);

   // update delta or its ratio accordingly 
  ValueAndRatio delta(m_edited.delta, m_edited.deltaRatio, lineEditDelta, lineEditDeltaRatio);
  delta.update(m_edited.epsilon, m_edited.deltaKeepRatio);

  // TODO: update size accordingly

}

void EditParametersDialog::changedDelta()
{
  // update delta ratio accordingly
  updateRatioFromValue(m_edited.delta, m_edited.epsilon, m_edited.deltaRatio);
  updateParamFromValue(m_edited.deltaRatio, lineEditDeltaRatio);
}

void EditParametersDialog::changedRadiusRatio()
{
  // update radius accordingly
  updateValueFromRatio(m_edited.radiusRatio, m_edited.bboxP, m_edited.radius);
  updateParamFromValue(m_edited.radius, lineEditRadius);

  // update epsilon or its ratio accordingly
  ValueAndRatio epsilon(m_edited.epsilon, m_edited.epsilonRatio, lineEditEpsilon, lineEditEpsilonRatio);
  epsilon.update(m_edited.radius, m_edited.epsilonKeepRatio);

  // update delta or its ratio accordingly (if needed)
  ValueAndRatio delta(m_edited.delta, m_edited.deltaRatio, lineEditDelta, lineEditDeltaRatio);
  delta.update(m_edited.epsilon, m_edited.deltaKeepRatio);
  
}

void EditParametersDialog::changedEpsilonRatio()
{
  // update epsilon accordingly
  updateValueFromRatio(m_edited.epsilonRatio, m_edited.radius, m_edited.epsilon);
  updateParamFromValue(m_edited.epsilon, lineEditEpsilon);

   // update delta or its ratio accordingly
  ValueAndRatio delta(m_edited.delta, m_edited.deltaRatio, lineEditDelta, lineEditDeltaRatio);
  delta.update(m_edited.epsilon, m_edited.deltaKeepRatio);

  // TODO: update size accordingly  
}

void EditParametersDialog::changedDeltaRatio()
{
  // update delta accordingly
  updateValueFromRatio(m_edited.deltaRatio, m_edited.epsilon, m_edited.delta);
  updateParamFromValue(m_edited.delta, lineEditDelta);
}
  
void EditParametersDialog::keepRadiusRatio(bool radiusRatio)
{
  m_edited.radiusKeepRatio = radiusRatio;
  
  lineEditRadius->setEnabled(!m_edited.radiusKeepRatio);
  lineEditRadiusRatio->setEnabled(m_edited.radiusKeepRatio);
}

void EditParametersDialog::keepEpsilonRatio(bool epsilonRatio)
{
  m_edited.epsilonKeepRatio = epsilonRatio;

  lineEditEpsilon->setEnabled(!m_edited.epsilonKeepRatio);
  lineEditEpsilonRatio->setEnabled(m_edited.epsilonKeepRatio);
}

void EditParametersDialog::keepDeltaRatio(bool deltaRatio)
{
  m_edited.deltaKeepRatio = deltaRatio;

  lineEditDelta->setEnabled(!m_edited.deltaKeepRatio);
  lineEditDeltaRatio->setEnabled(m_edited.deltaKeepRatio);
}

#include "EditParametersDialog.moc"
