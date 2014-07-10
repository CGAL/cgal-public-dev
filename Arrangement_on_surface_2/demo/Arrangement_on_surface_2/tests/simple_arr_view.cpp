#include <iostream>

#include <QObject>
#include <QApplication>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QMainWindow>
#include <QWidget>
#include <QGridLayout>
#include <QGraphicsRectItem>

#include "PolynomialParser.h"
#include "AlgebraicDemoTraits.h"

typedef AlgebraicDemoTraits::ArrTraitsType ArrTraitsType;
typedef AlgebraicDemoTraits::ArrangementType ArrangementType;
typedef CGAL::Qt::ArrangementGraphicsItem< ArrangementType > ArrGraphicsItemType;

QMainWindow* window;
QGraphicsView* view;
QGraphicsScene* scene;
ArrGraphicsItemType* arrItem;

void SetupUI( )
{
    scene = new QGraphicsScene( -10, -10, 20, 20 );
    view = new QGraphicsView;
    view->setScene( scene );
    QMatrix mat( 1, 0, 0, -1, 0, 0 );
    view->setMatrix( mat );

    // does nothing
    //view->centerOn( -8, -5 );
    //view->centerOn( box );

    view->scale( 10, 10 );

    QMainWindow* window = new QMainWindow;
    window->resize( 800, 600 );
    QWidget* centralWidget = new QWidget(window);
    QGridLayout* layout = new QGridLayout(centralWidget);
    layout->addWidget(view, 0, 0, 1, 1);
    window->setCentralWidget(centralWidget);
    window->show( );
}

int main( int argc, char *argv[] )
{
    QApplication app( argc, argv );

    if ( argc < 2 )
    {
      std::cout << "Usage: " << argv[0] << " alg-dat-file\n";
      return 1;
    }

    ArrangementType arr;
    LoadArrFromFile< AlgebraicDemoTraits > loadArr;
    loadArr( argv[1], &arr );

    SetupUI( );

    arrItem = new ArrGraphicsItemType( &arr );
    arrItem->setScene( scene );
    scene->addItem( arrItem );
    //QGraphicsRectItem* box = new QGraphicsRectItem( 0, 0, 10, 10 ) ;
    //scene->addItem( box );

    return app.exec( );
}
