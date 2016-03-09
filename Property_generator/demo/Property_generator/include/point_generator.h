/******************************************************************************/
#ifndef POINT_GENERATOR_H
#define POINT_GENERATOR_H
/******************************************************************************/

#include <QDialog>

/******************************************************************************/

namespace Ui {
class Point_Generator;
}

/******************************************************************************/

class Point_Generator : public QDialog
{
    Q_OBJECT

public:
    explicit Point_Generator(QWidget *parent = 0);
    ~Point_Generator();

signals:
    void build_triangulation(int size);
    void build_triangulation_fractal(int start_size, int iterations);
    void build_triangulation(QString filename);
    void build_triangulation_poisson( double rate, double side );

public slots:
    void on_buttonBox_accepted();

private:
    Ui::Point_Generator *ui;
};

/******************************************************************************/
#endif // POINT_GENERATOR_H
/******************************************************************************/
