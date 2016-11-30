#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QGraphicsScene>
#include <CGAL/Qt/GraphicsViewNavigation.h>

#include "EventFilterManagerGroup.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
	explicit MainWindow(QWidget* parent = 0);
	~MainWindow();

private:
	Ui::MainWindow* ui;

	EventFilterManagerGroup* manager;

	void setScene();
	void setEventFilterManagerGroup();
	void setViewControl();
	void setPolygonInputDevice();
	void setToolbarActionGroup();
	void addToolbarAction();
	void addToolbarActionWithActionGroup();


private slots:

};

#endif // MAINWINDOW_H
