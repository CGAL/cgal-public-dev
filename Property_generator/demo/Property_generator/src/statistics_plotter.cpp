/******************************************************************************/

// Local.
#include "statistics_plotter.h"

// Qwt.
#include <qwt_plot_layout.h>
#include <qwt_legend.h>
#include <qwt_legend_item.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_histogram.h>
#include <qwt_column_symbol.h>
#include <qwt_series_data.h>
#include <qwt_color_map.h>
#include <qwt_scale_widget.h>
#include <qwt_symbol.h>

// Qt.
#include <qDebug>
#include <qpen.h>
#include <algorithm>

#include <iostream>

/******************************************************************************/

Scatter_Plot_Color::Scatter_Plot_Color(const QString &title):
    QwtPlotSpectroCurve(title)
{
    // setStyle(QetPlotCurve::NoCurve);
}

/******************************************************************************/

void Scatter_Plot_Color::setValues( std::vector<double>* values1,
                                    std::vector<double>* values2,
                                    std::vector<double>* values3  )
{
    QVector<QwtPoint3D> samples;

    double min = values3->at(0);
    double max;

    for (int i=0; i<values1->size(); i++)
    {
        QwtPoint3D point;

        point.setX( values1->at(i) );
        point.setY( values2->at(i) );
        point.setZ( values3->at(i) );

        if (values3->at(i) > max)
            max = values3->at(i);

        if (values3->at(i) < min)
            min = values3->at(i);

        samples.push_back(point);
    }

    interval.setMinValue(min);
    interval.setMaxValue(max);

    setColorMap(new ColorMap);
    setColorRange(interval);
    setPenWidth(3);
    setSamples(samples);
}

/******************************************************************************/

QwtInterval Scatter_Plot_Color::getInterval()
{
    return interval;
}

/******************************************************************************/

Scatter_Plot::Scatter_Plot(const QString &title, const QColor &symbolColor):
    QwtPlotCurve(title)
{
    setColor(symbolColor);
    setStyle(QwtPlotCurve::NoCurve);
}

/******************************************************************************/

void Scatter_Plot::setColor(const QColor &symbolColor)
{
    QColor color = symbolColor;
    color.setAlpha(180);
    setSymbol(new QwtSymbol( QwtSymbol::Ellipse,
                             color,
                             QPen(color),
                             QSize(4, 4 ) ) );
}

/******************************************************************************/

void Scatter_Plot::setValues( std::vector<double>* values1,
                              std::vector<double>* values2  )
{
    if (values1->size() != values2->size())
    {
        std::cerr << "Value lengths do not match" << endl;
        exit(1);
    }

    setSamples(&(*values1)[0], &(*values2)[0], values1->size());
}

/******************************************************************************/

Histogram::Histogram(const QString &title, const QColor &symbolColor):
    QwtPlotHistogram(title)
{
    setStyle(QwtPlotHistogram::Columns);
    setColor(symbolColor);
}

/******************************************************************************/

void Histogram::setColor(const QColor &symbolColor)
{
    QColor color = symbolColor;
    color.setAlpha(180);

    setPen(QPen(Qt::black));
    setBrush(QBrush(color));

    QwtColumnSymbol *symbol = new QwtColumnSymbol(QwtColumnSymbol::Box);
    symbol->setFrameStyle(QwtColumnSymbol::Raised);
    symbol->setLineWidth(2);
    symbol->setPalette(QPalette(color));
    setSymbol(symbol);
}

/******************************************************************************/

void Histogram::setValues(std::vector<double> *values)
{


    if (values == NULL)
    {
        std::cerr << "Values are null" << endl;
        return;
    }


    // Rule of thumb for choosing the number of bins.
    int num_bins = ceil(sqrt(values->size()));

    if (num_bins == 0) return;

    std::vector<double>::iterator i;

    // Get max and min values.
    i = std::max_element(values->begin(), values->end());
    double max = *i;
    i = std::min_element(values->begin(), values->end());
    double min = *i;

    // Put the data into bins.
    double bin_width = (max-min) / (float)num_bins;
    double bin[num_bins];

    for (int i=0; i<num_bins; i++)
        bin[i] = 0;

    // Ignore 'stupid' plots.
    if (bin_width==0) return;

    for (i=values->begin(); i!= values->end(); ++i)
    {
        int index = (int)floor((*i-min)/bin_width);

        if (! (index < 0 || index >= num_bins) )
            bin[index] ++;
    }

    // Create a set of samples for the histogram.
    QVector<QwtIntervalSample> samples(num_bins);
    for (int i=0; i<num_bins; ++i)
        samples[i] = QwtIntervalSample( bin[i],
                                        min + bin_width* i,
                                        min + bin_width*(i+1) );

    setData( new QwtIntervalSeriesData(samples) );
}

/******************************************************************************/

Statistics_plotter::Statistics_plotter(QWidget *parent) : QwtPlot(parent)
{

    // Destroy this object on close.
    this->setAttribute(Qt::WA_DeleteOnClose);

    // Display this as a window, and keep it on top.
    this->setWindowFlags(Qt::WindowStaysOnTopHint | Qt::Window);

    setTitle("Triangulation Statistics");

    setCanvasBackground(QColor(Qt::white));
    plotLayout()->setAlignCanvasToScales(true);

    QwtPlotGrid *grid = new QwtPlotGrid;
    grid->enableX(false);
    grid->enableY(true);
    grid->enableXMin(false);
    grid->enableYMin(false);
    grid->setMajPen(QPen(Qt::black, 0, Qt::DotLine));
    grid->attach(this);

    setAutoReplot(true);
}

/******************************************************************************/

void Statistics_plotter::add_color_bar( QString title, QwtInterval zInterval )
{
    QwtScaleWidget *rightAxis = axisWidget(QwtPlot::yRight);
    rightAxis->setTitle(title);
    rightAxis->setColorBarEnabled(true);
    rightAxis->setColorMap( zInterval, new ColorMap );

    setAxisScale(QwtPlot::yRight, zInterval.minValue(), zInterval.maxValue() );
    enableAxis(QwtPlot::yRight);
}

/******************************************************************************/

void Statistics_plotter::add_legend()
{
    legend = new QwtLegend;
    legend->setItemMode(QwtLegend::CheckableItem);
    insertLegend(legend, QwtPlot::RightLegend);

    connect(this, SIGNAL(legendChecked(QwtPlotItem *, bool)),
            this, SLOT(showItem(QwtPlotItem *, bool)));

    replot();

    QwtPlotItemList items = itemList(QwtPlotItem::Rtti_PlotItem);

    bool first=true;
    for ( int i = 0; i < items.size(); i++ )
    {
        // Select the first curve or histogram and show it.
        int rtti = items[i]->rtti();
        if (first && (    rtti == (int)QwtPlotItem::Rtti_PlotCurve
                       || rtti == (int)QwtPlotItem::Rtti_PlotSpectroCurve
                       || rtti == (int)QwtPlotItem::Rtti_PlotHistogram ) )
        {
            first=false;

            QwtLegendItem *legendItem;
            legendItem = (QwtLegendItem *)legend->find(items[i]);
            if ( legendItem )
                legendItem->setChecked(true);
            items[i]->setVisible(true);
        }
        else
            items[i]->setVisible(false);
    }
}

/******************************************************************************/

void Statistics_plotter::showItem(QwtPlotItem *item, bool on)
{
    item->setVisible(on);
}

/******************************************************************************/

void Statistics_plotter::set_titles(QString x_title, QString y_title)
{
    setAxisTitle(QwtPlot::yLeft, y_title);
    setAxisTitle(QwtPlot::xBottom, x_title);
}

/******************************************************************************/
