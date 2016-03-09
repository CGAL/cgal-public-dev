/*******************************************************************************
* Written by Ross Hemsley for INRIA.fr.
* A simple application visualise different walks on Delaunay Triangulations.
*******************************************************************************/

#include <QApplication>
#include "visualiser.h"

/******************************************************************************/

int main(int argc, char **argv)
{
    QApplication app(argc, argv);

    app.setOrganizationDomain( "inria.fr" );
    app.setOrganizationName  ( "INRIA" );
    app.setApplicationName   ( "Delaunay Visualiser" );

    Visualiser visualiser;
    visualiser.show();

    return app.exec();
}

/******************************************************************************/

