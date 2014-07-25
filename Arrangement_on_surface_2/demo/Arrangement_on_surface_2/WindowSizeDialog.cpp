// Copyright (c) 2014  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "WindowSizeDialog.h"
#include "ui_WindowSizeDialog.h"

WindowSizeDialog::WindowSizeDialog( QWidget* parent, Qt::WindowFlags f ):
  ui( new Ui::WindowSizeDialog )
{
  this->ui->setupUi( this );
}

QRectF WindowSizeDialog::size( ) const
{
  bool ok1 = false;
  bool ok2 = false;
  bool ok3 = false;
  bool ok4 = false;
  QRectF res;

  double xmin = this->ui->xMin->text( ).toDouble( &ok1 );
  double ymin = this->ui->yMin->text( ).toDouble( &ok2 );
  double xmax = this->ui->xMax->text( ).toDouble( &ok3 );
  double ymax = this->ui->yMax->text( ).toDouble( &ok4 );

  if ( ok1 && ok2 && ok3 && ok4 )
  {
    res = QRectF( xmin, ymin, xmax - xmin, ymax - ymin );
  }

  return res;
}
