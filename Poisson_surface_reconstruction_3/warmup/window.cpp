#include "window.h"
#include "viewer.h"
#include "scene.h"
#include "options.h"

#include "ui_mesh.h"
#include "ui_options.h"

#include <CGAL/Qt/debug.h>

#include <QTextStream>
#include <QUrl>
#include <QFileDialog>
#include <QSettings>
#include <QHeaderView>

MainWindow::MainWindow(QWidget* parent)
	: CGAL::Qt::DemosMainWindow(parent)
{
	ui = new Ui::MainWindow;
	ui->setupUi(this);

	// default parameters
	m_angle = 30.0;
	m_sizing = 0.15;
	m_approximation = 0.015;

	// saves some pointers from ui, for latter use.
	m_pViewer = ui->viewer;
	assert(m_pViewer != NULL);

	// does not save the state of the viewer
	m_pViewer->setStateFileName(QString::null);

	// setups scene
	m_pScene = new Scene();
	m_pViewer->setScene(m_pScene);

	// connects actionQuit (Ctrl+Q) and qApp->quit()
	connect(ui->actionQuit, SIGNAL(triggered()), this, SLOT(quit()));

	this->addRecentFiles(ui->menuFile, ui->actionQuit);

	readSettings();
}

MainWindow::~MainWindow()
{
	delete m_pScene;
	delete ui;
}

void MainWindow::quit()
{
	writeSettings();
	close();
}

void MainWindow::readSettings()
{
	this->readState("MainWindow", Size|State);
}

void MainWindow::writeSettings()
{
	this->writeState("MainWindow");
	std::cerr << "Write setting...done.\n";
}

void MainWindow::on_actionLoad_oriented_points_triggered(){
	QString fileName = QFileDialog::getOpenFileName(this,tr("Open Image file"), "",tr("XY file (*.xyz);;All Files (*)"));
	//std::cout << read_xy(fileName) << std::endl;
	m_pScene->read_xyz(fileName);
	m_pViewer->update();
}

void MainWindow::updateViewerBBox()
{
	const Bbox bbox(-2.0, -2.0, -2.0,
		2.0,  2.0,  2.0);
	const double xmin = bbox.xmin();
	const double ymin = bbox.ymin();
	const double zmin = bbox.zmin();
	const double xmax = bbox.xmax();
	const double ymax = bbox.ymax();
	const double zmax = bbox.zmax();
	qglviewer::Vec vec_min(xmin, ymin, zmin);
	qglviewer::Vec vec_max(xmax, ymax, zmax);
	m_pViewer->setSceneBoundingBox(vec_min, vec_max);
	m_pViewer->camera()->showEntireScene();
}




void MainWindow::on_actionView_vertices_triggered()
{
	m_pScene->toggle_view_vertices();
	m_pViewer->update();
}

void MainWindow::on_actionView_edges_triggered()
{
	m_pScene->toggle_view_edges();
	m_pViewer->update();
}

void MainWindow::on_actionView_mesh_triggered()
{
	m_pScene->toggle_view_mesh();
	m_pViewer->update();
}

void MainWindow::on_actionParameters_triggered()
{
	Dialog_options dlg;

	dlg.set_angle(m_angle);
	dlg.set_sizing(m_sizing);
	dlg.set_approximation(m_approximation);

	if(dlg.exec() == QDialog::Accepted)
	{
		m_angle = dlg.get_angle();
		m_sizing = dlg.get_sizing();
		m_approximation = dlg.get_approximation();
	}
}

void MainWindow::on_actionMesh_sphere_triggered()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);
	m_pScene->mesh_sphere(m_angle, m_sizing, m_approximation);
	QApplication::restoreOverrideCursor();
	m_pViewer->update();
}

void MainWindow::on_actionMesh_torus_triggered()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);
	m_pScene->mesh_torus(m_angle, m_sizing, m_approximation);
	QApplication::restoreOverrideCursor();
	m_pViewer->update();
}

void MainWindow::on_actionMesh_ellipsoid_triggered()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);
	m_pScene->mesh_ellipsoid(m_angle, m_sizing, m_approximation);
	QApplication::restoreOverrideCursor();
	m_pViewer->update();
}

void MainWindow::on_actionMesh_with_implicit_triggered(){
	m_pScene->implicit_function();
	m_pViewer->update();
}
