#include "MainWindow.h"
#include "ui_MainWindow.h"

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
	setAcceptDrops(true);

	setScene();
	setEventFilterManagerGroup();
	setViewControl();
	setPolygonInputDevice();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setScene()
{
	QGraphicsScene* scene = new QGraphicsScene();
	scene->setItemIndexMethod(QGraphicsScene::NoIndex);
	scene->setSceneRect(-100, -100, 100, 100);
	ui->SceneGraphicsView->setScene(scene);
	ui->SceneGraphicsView->setMouseTracking(true);
}

void MainWindow::setEventFilterManagerGroup()
{
	manager = new EventFilterManagerGroup();
	manager->addObjectToWatch("viewport", ui->SceneGraphicsView->viewport());
	manager->addObjectToWatch("scene", ui->SceneGraphicsView->scene());
}

void MainWindow::setViewControl()
{
	manager->addFilterWidget("viewport", "navigate", new CGAL::Qt::GraphicsViewNavigation());
}

void MainWindow::setPolygonInputDevice()
{
	//manager->addFilterWidget("scene", "draw",
	//	new PolylineInput(ui->SceneGraphicsView, ui->SceneGraphicsView->scene, 0, false));
}

void MainWindow::setToolbarActionGroup()
{

}

void MainWindow::addToolbarAction()
{

}

void MainWindow::addToolbarActionWithActionGroup()
{

}
