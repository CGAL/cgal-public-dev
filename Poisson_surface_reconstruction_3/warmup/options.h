#ifndef DIALOG_OPTIONS
#define DIALOG_OPTIONS 1

#include "ui_options.h"

class Dialog_options : public QDialog, private Ui::Dialog_options
{
    Q_OBJECT
public:
    Dialog_options(QWidget *parent = 0)
    {
        setupUi(this);
    }
    
    double get_angle() const { return angle_spinbox->value(); }
    void set_angle(const double angle) { angle_spinbox->setValue(angle); }

		double get_sizing() const { return sizing_spinbox->value(); }
    void set_sizing(const double sizing) { sizing_spinbox->setValue(sizing); }

		double get_approximation() const { return approximation_spinbox->value(); }
    void set_approximation(const double approximation) { approximation_spinbox->setValue(approximation); }
};

#endif
