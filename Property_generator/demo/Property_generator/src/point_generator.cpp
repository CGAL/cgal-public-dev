/******************************************************************************/
#include "point_generator.h"
#include "ui_point_generator.h"
#include <qdebug>
/******************************************************************************/

Point_Generator::Point_Generator(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Point_Generator)
{
    ui->setupUi(this);
}

/******************************************************************************/

void Point_Generator::on_buttonBox_accepted()
{
    if (ui->radio_fractal->isChecked())
    {
        emit build_triangulation_fractal( ui->spin_fractal_start->value(),
                                          ui->spin_fractal_iterations->value());
    }

    if (ui->radio_uniform->isChecked())
    {
        emit build_triangulation(ui->spin_uniform_points->value());
    }

    if (ui->radio_poisson->isChecked())
    {
        emit build_triangulation_poisson( ui->spin_poisson_rate->value(),
                                          ui->spin_poisson_side->value() );
    }
}

/******************************************************************************/

Point_Generator::~Point_Generator()
{
    delete ui;
}

/******************************************************************************/
