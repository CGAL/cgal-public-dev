#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <CGAL/Qt/DemosMainWindow.h>

class QDragEnterEvent;
class QDropEvent;

class Viewer;
class Scene;

namespace Ui {
	class MainWindow;
}

class MainWindow : public CGAL::Qt::DemosMainWindow
{
	Q_OBJECT
private:
	Ui::MainWindow* ui;

	// parameters
	double m_angle;
	double m_sizing;
	double m_approximation;

	Viewer* m_pViewer;
	Scene* m_pScene;

public:
	MainWindow(QWidget* parent = 0);
	~MainWindow();

	public slots:
		void updateViewerBBox();

		protected slots:

			// settings
			void quit();
			void readSettings();
			void writeSettings();
			void on_actionLoad_oriented_points_triggered();

			// view menu
			void on_actionView_mesh_triggered();
			void on_actionView_edges_triggered();
			void on_actionView_vertices_triggered();

			// algorithms menu
			void on_actionParameters_triggered();
			void on_actionMesh_sphere_triggered();
			void on_actionMesh_torus_triggered();
			void on_actionMesh_ellipsoid_triggered();
			void on_actionMesh_with_implicit_triggered();

};

#endif // ifndef MAINWINDOW_H
