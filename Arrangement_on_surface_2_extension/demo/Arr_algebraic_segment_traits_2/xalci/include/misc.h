// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : AlciX
// File          : demos/xalci/include/misc.h
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.2 $
// Revision_date : $Date: 2009-06-30 13:14:59 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef MISC_H
#define MISC_H

#include <qlistbox.h>
#include <qpushbutton.h>
#include <qlayout.h>
#include <qvbox.h>

#include <sstream>

#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>

/*!\brief
 * declares miscellaneous routines in a separate file (to speed-up
 * compilation)
 */

class Graphic_layer : public CGAL::Qt_widget_layer
{
Q_OBJECT
 
public:
    Graphic_layer(CGAL::Qt_widget *attach_to, int index_,
        int color_index_, QObject* parent = NULL, const char* name = NULL) :
         Qt_widget_layer(parent, name),
          erase(false), index(index_), color_index(color_index_) { 

        attach_to->attach(this);
    }

    void draw();

    bool erase;

protected:
    int index, color_index;
       
};

typedef std::vector<Graphic_layer *> Layers;

#endif // MISC_H
