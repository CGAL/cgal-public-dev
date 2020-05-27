#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_clippingplane.h"

class ClippingPlane : public QMainWindow
{
    Q_OBJECT

public:
    ClippingPlane(QWidget *parent = Q_NULLPTR);

private:
    Ui::ClippingPlaneClass ui;
};
