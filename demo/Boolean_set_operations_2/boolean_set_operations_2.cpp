// Copyright (c) 2012  Tel-Aviv University (Israel).
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
// Author(s) : Saar Katz <kats.saar@gmail.com>


// When in windows, this is required to prevent the exe from being console application (/SUBSYSTEM:windows)
// while keeping the entry point "main" and not "WinMain" as it is the default of windows applications (/ENTRY:mainCRTStartup)
// These lines can be replaced with the line
//		set_target_properties(Target PROPERTIES LINK_FLAGS "/SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
// in a CMakeLists.txt file
// A reason for having this code here is to keep the CMakeLists.txt free of platform specific code.
// A reason for having this code in a CMakeLists.txt file if for the "fix" to be applied for all the demos.
//
// ===============================================++++++++++++++++++
// Sometihng isn't  working with the Release build (for me that is)
// ===============================================++++++++++++++++++
// For now the console is useful for debugging.
//#if defined(_WIN32) && !defined(_DEBUG)
//#include <windows.h>
//#pragma comment(linker, "/SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
//#endif

#include <QApplication>
#include <qmessagebox.h>

#include "MainWindow.h"

void show_warning(std::string aS)
{
  QMessageBox::warning(NULL, "Warning", QString(aS.c_str()));
}

void show_error(std::string aS)
{
  QMessageBox::critical(NULL, "Critical Error", QString(aS.c_str()));
}

void error(std::string aS)
{
  show_error(aS);

  throw std::runtime_error(aS);
}


int main(int argc, char* argv[])
{
  QApplication a(argc, argv);
  try
  {
    MainWindow w;
    w.show();

    return a.exec();
  }
  catch (const std::exception e)
  {
    std::string s = e.what();
    show_error("Exception throne during run of the program:\n" + s);
  }
}
