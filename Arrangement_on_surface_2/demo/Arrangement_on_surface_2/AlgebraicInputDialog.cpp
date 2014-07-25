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

#include "AlgebraicInputDialog.h"
#include "ui_AlgebraicInputDialog.h"

AlgebraicInputDialog::AlgebraicInputDialog( QWidget* parent, Qt::WindowFlags f ):
  ui( new Ui::AlgebraicInputDialog )
{
  this->ui->setupUi( this );
}

QString AlgebraicInputDialog::line( int index ) const
{
  if ( index < 0 || index >= 10 )
  {
    return "";
  }

  switch ( index )
  {
    case 0: return this->ui->lineEdit->text();
    case 1: return this->ui->lineEdit_2->text();
    case 2: return this->ui->lineEdit_3->text();
    case 3: return this->ui->lineEdit_4->text();
    case 4: return this->ui->lineEdit_5->text();
    case 5: return this->ui->lineEdit_6->text();
    case 6: return this->ui->lineEdit_7->text();
    case 7: return this->ui->lineEdit_8->text();
    case 8: return this->ui->lineEdit_9->text();
    case 9: return this->ui->lineEdit_10->text();
    default: return ""; // never get here
  }
}
