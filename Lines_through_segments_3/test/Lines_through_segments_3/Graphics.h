#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <QtGui/QGraphicsScene>
#include <QtGui/QApplication>
#include <QtGui/QGraphicsLineItem>
#include <QtGui/QGraphicsView>
#include <CGAL/basic.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>

const double CROSS_SIZE = 2;
#define SCREEN_SIZE 50

//forward declerations
template <typename Number_type_,typename Segment_2, typename Point_2>
class Graphics
{
public:
   typedef Number_type_ Number_type;
   
   Graphics(int argc, char* argv[],Number_type min_x = -10, Number_type max_x = 10, Number_type min_y = -10, Number_type max_y = 10)
      :app(argc, argv), scene(-1,-1, 3, 3)
   {
   }

   Graphics()
      :app(0, NULL)
   {
   }
   
   //Graphics();
   
   ~Graphics(void)
   {
   }
   
   void Display()
   {
      static int count = 0;
      std::stringstream str;
      str << "img";
      str << count;
      str << ".png";
      
//      std::cout << str.str() << std::endl;
      QGraphicsView* view = new QGraphicsView(&this->scene);
      CGAL::Qt::GraphicsViewNavigation navigation;
      view->installEventFilter(&navigation);
      view->viewport()->installEventFilter(&navigation);
      view->setRenderHint(QPainter::Antialiasing);

      // a white background
      scene.setBackgroundBrush(Qt::white);
//    view->show();
      

      QImage img(1024,768,QImage::Format_ARGB32_Premultiplied);
      QPainter p(&img);
      scene.render(&p);
      p.end();
      count++;
      img.save(str.str().c_str());

//      app.exec();

      
      return;
   }


   QGraphicsLineItem* draw_edge(const Point_2& p1, const Point_2& p2,char* color)
   {
      QPen pen;
      //      std::cout << "draw_edge (" << p1 << "," << p2 << ")" << std::endl;
      if (strcmp(color,"blue") == 0)
         pen = QColor("blue");//CMYK
      else
         pen = QColor("black");//CMYK
      pen.setWidthF(0.01);
      QGraphicsLineItem * lineItem = scene.addLine(QLineF( CGAL::to_double(p1.x()),
                                                           CGAL::to_double(p1.y()),
                                                           CGAL::to_double(p2.x()),
                                                           CGAL::to_double(p2.y())),pen);
      return lineItem;
   }

   void add_text(const char *text, double x, double y)
   {
      // QGraphicsTextItem* text_item = scene.addText(QString(text));//, const QFont & font = QFont() )
      // text_item->setPos(x,y);
      // QFont font("Times",12,true);
      // font.setItalic(true);
      // text_item->setFont(font);
      // text_item->setDefaultTextColor(QColor("darkRed"));//CMYK
   }
   
   //circle drawing functions
   void draw_cross(double x, double y,double size = CROSS_SIZE)
   {
      draw_edge(Point_2(x - size, y), 
                Point_2(x + size, y), "black");
      draw_edge(Point_2(x, y + size), 
                Point_2(x, y - size), "black");
   }


private:
   QApplication app;
   //QApplication* app_p;
   QGraphicsScene scene;
   
   
};

#endif
