#ifndef DIAL_SETTINGS_H
#define DIAL_SETTINGS_H

#include <QtGui>
#include <QtWidgets>
#include "ui_dialSettings.h"

class SettingsDialog :
  public QDialog,
  public Ui_dialSettings
{
    Q_OBJECT

public:
    SettingsDialog(QWidget *parent = 0);
    void accept();

private slots:
    void loadDefaults();
    void loadFromSettings();
    void saveToSettings();

private:
};

#endif // DIAL_SETTINGS_H
