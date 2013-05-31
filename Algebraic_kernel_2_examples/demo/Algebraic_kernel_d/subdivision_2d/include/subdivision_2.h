
#ifndef SUBDIVISION_2_H
#define SUBDIVISION_2_H

#include <qapplication.h>
#include <qmainwindow.h>

#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>

#include <qlayout.h>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qtextedit.h>
#include <qpixmap.h>
#include <qpainter.h>
#include <qvbox.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>

#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>

typedef CGAL::CORE_arithmetic_kernel AK;
typedef AK::Rational Rational;
typedef AK::Integer Integer;

typedef CGAL::Polynomial_type_generator< Integer, 2 >::Type
            Poly_int_2;

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

    int index, color_index; 
};

class Subdiv_main_window : public QMainWindow
{
    Q_OBJECT

public:

    Subdiv_main_window(int w, int h) {
        setup(w, h);
    }

    void visualize();

public slots:

    void rasterize_click();

    void new_instance()
    {
        widget->lock();
        widget->clear_history();
        widget->set_window(-1.1, 1.1, -1.1, 1.1);
        // set the Visible Area to the Interval
        widget->unlock();
    }

protected slots:

    void about()
    {
        QMessageBox::about(this, "About",
            "This program demonstrates the use of Segment renderern"
            "and space subdivision to visualize algebraic curvesn");
    }
    void aboutQt()
    {
        QMessageBox::aboutQt(this, "About Qt");
    }
    void howto()
    {
        CGAL::Qt_help_window *help =
          new CGAL::Qt_help_window("help/index.html", ".", 0, "help viewer");
        help->resize(400, 400);
        help->setCaption("Demo HowTo");
        help->show();
    }
    void new_window()
    {  }

    void setup(int, int);
        
protected:

    bool input_poly(Poly_int_2& p, const char *ascii);
        
    QPushButton *rasterize_btn;
    QTextEdit *tex_input;
        
    QWidget *central_widget;
    CGAL::Qt_widget *widget;
    CGAL::Qt_widget_standard_toolbar *stoolbar;
};

#endif
