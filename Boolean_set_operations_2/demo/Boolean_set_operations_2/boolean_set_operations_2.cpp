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

#include "MainWindow.h"
#include <QApplication>
#include <qmessagebox.h>

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
