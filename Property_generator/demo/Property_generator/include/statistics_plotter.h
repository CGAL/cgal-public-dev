/******************************************************************************/
#ifndef STATISTICS_PLOTTER_H
#define STATISTICS_PLOTTER_H
/******************************************************************************/

#include <stdlib.h>
#include <vector>

// Qt.
#include <qpen.h>

// Qwt.
#include <qwt_plot_layout.h>
#include <qwt_legend.h>
#include <qwt_legend_item.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_histogram.h>
#include <qwt_column_symbol.h>
#include <qwt_series_data.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_spectrocurve.h>
#include <qwt_color_map.h>

/******************************************************************************/

class ColorMap: public QwtLinearColorMap
{
public:
    ColorMap():
        QwtLinearColorMap(Qt::darkCyan, Qt::red)
    {
        addColorStop(0.2, Qt::cyan);
        addColorStop(0.6, Qt::green);
        addColorStop(0.7, Qt::yellow);
    }
};

/******************************************************************************/

class Histogram: public QwtPlotHistogram
{
public:
    Histogram(const QString &, const QColor &);
    void setValues(std::vector<double> *values);
    void setColor(const QColor &);
};

/******************************************************************************/

class Scatter_Plot: public QwtPlotCurve
{
public:
    Scatter_Plot (const QString &, const QColor &);
    void setColor(const QColor &);
    void setValues(std::vector<double> *values1, std::vector<double> *values2);
};

/******************************************************************************/

class Scatter_Plot_Color: public QwtPlotSpectroCurve
{
private:
    QwtInterval interval;

public:
    Scatter_Plot_Color(const QString &);
    QwtInterval getInterval();
    void setValues( std::vector<double>* values1,
                    std::vector<double>* values2,
                    std::vector<double>* values3  );
};

/******************************************************************************/

class Statistics_plotter : public QwtPlot
{
    Q_OBJECT
public:

    Statistics_plotter(QWidget *parent = NULL);
    void set_titles(QString x_title, QString y_title);
    void add_legend();
    void add_color_bar( QString title, QwtInterval zInterval);

private slots:
    void   showItem(QwtPlotItem *item, bool on);

private:
    QwtLegend *legend;
};

/******************************************************************************/
#endif // STATISTICS_PLOTTER_H
/******************************************************************************/
