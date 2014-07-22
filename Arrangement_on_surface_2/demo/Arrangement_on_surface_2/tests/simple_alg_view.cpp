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

void SetupUI( const QRectF& windowRect, double zoom )
{
    scene = new QGraphicsScene( windowRect.left(),
        windowRect.bottom(),
        windowRect.width(),
        windowRect.height() );
    view = new QGraphicsView;
    view->setScene( scene );
    QMatrix mat( 1, 0, 0, -1, 0, 0 );
    view->setMatrix( mat );

    view->scale( zoom, zoom );

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
    std::string curveFileName;
    QRectF windowRect;

    if ( argc < 7 )
    {
      std::cout << "Usage: " << argv[0] << " [curve_file] [x y w h] [zoom]\n";
      return 1;
    }

    curveFileName = argv[1];
    windowRect = QRectF( atof( argv[2] ),
        atof( argv[3] ),
        atof( argv[4] ),
        atof( argv[5] ) );

    double zoom = atof( argv[6] );

    ArrangementType arr;

    // load the curve
    LoadArrFromFile< AlgebraicDemoTraits > loadArr;
    loadArr( curveFileName, &arr );

    SetupUI( windowRect, zoom );

    arrItem = new ArrGraphicsItemType( &arr );
    scene->addItem( arrItem );
    //arrItem->modelChanged( );

    QGraphicsRectItem* box = new QGraphicsRectItem( 0, 10, 10, 10 ) ;
    scene->addItem( box );

    QGraphicsLineItem* x_axis = new QGraphicsLineItem( 0, -10, 0, 10 );
    QGraphicsLineItem* y_axis = new QGraphicsLineItem( -10, 0, 10, 0 );
    scene->addItem( x_axis );
    scene->addItem( y_axis );
    arrItem->modelChanged( );

    return app.exec( );
}
