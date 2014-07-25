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

#ifndef WINDOW_SIZE_DIALOG_H
#define WINDOW_SIZE_DIALOG_H

#include <QDialog>

class QButtonGroup;
namespace Ui
{
  class WindowSizeDialog;
}

class WindowSizeDialog : public QDialog
{
public:
  WindowSizeDialog( QWidget* parent = 0, Qt::WindowFlags f = 0 );

  /**
  \return the user-specified window size as a QRectF. It is null if invalid.
  */
  QRectF size( ) const;

protected:
  Ui::WindowSizeDialog* ui;
};
#endif // WINDOW_SIZE_DIALOG_H
