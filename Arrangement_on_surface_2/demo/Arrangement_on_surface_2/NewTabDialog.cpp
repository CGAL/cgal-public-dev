// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "NewTabDialog.h"
#include "ArrangementDemoWindow.h"
#include "ui_NewTabDialog.h"

//!	A constructor for the dialog box referencing the radio buttons for each traits to be added onto the dialog box
/*!
  	\param parent A QWidget pointer to the class
  	\param f A QT::WindowFlags config variable
*/
NewTabDialog::NewTabDialog( QWidget* parent, Qt::WindowFlags f ) :
  QDialog( parent, f ),
  ui( new Ui::NewTabDialog ),
  buttonGroup( new QButtonGroup )
{
  this->ui->setupUi( this );

  this->buttonGroup->addButton( this->ui->segmentRadioButton,
								ArrangementDemoWindow::SEGMENT_TRAITS );
  this->buttonGroup->addButton( this->ui->polylineRadioButton,
								ArrangementDemoWindow::POLYLINE_TRAITS );
  this->buttonGroup->addButton( this->ui->conicRadioButton,
								ArrangementDemoWindow::CONIC_TRAITS );
  this->buttonGroup->addButton( this->ui->linearRadioButton,
								ArrangementDemoWindow::LINEAR_TRAITS );
  this->buttonGroup->addButton( this->ui->circularArcRadioButton,
								ArrangementDemoWindow::CIRCULAR_ARC_TRAITS );
  this->buttonGroup->addButton( this->ui->algebraicRadioButton,
								ArrangementDemoWindow::ALGEBRAIC_TRAITS );
  this->buttonGroup->addButton( this->ui->bezierRadioButton,
                                ArrangementDemoWindow::BEZIER_TRAITS );
}

//!	returns the ID of the checked radio button
/*!
  	\return interger value of the checked radio button
*/
int NewTabDialog::checkedId( ) const
{
  return this->buttonGroup->checkedId( );
}
