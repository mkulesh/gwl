/*******************************************************************************
 * GWL - Geophysical Wavelet Library
 * *****************************************************************************
 * Copyright (C) 2002-2017 Mikhail Kulesh, Matthias Holschneider
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *  
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#ifndef UTMATRIXOFQWTPLOT_H_
#define UTMATRIXOFQWTPLOT_H_

#include <qwidget.h>
#include <qlayout.h>
#include <qwt_plot.h>
#include <qwt_plot_layout.h>
#include <qwt_plot_canvas.h>

#include <vector.h>

class UTMatrixOfQwtPlot : public QWidget
{
    
private:

    std::vector<QwtPlot *> _plot;

    unsigned _rows, _cols;

    QGridLayout * _layout;

public:

    UTMatrixOfQwtPlot (QWidget * parent, unsigned rows, unsigned cols);

    QwtPlot *
    plot (unsigned row, unsigned col);

    QGridLayout *
    layout ()
    {
        return _layout;
    }
    
    void setRowAxisRange (unsigned row, double ymin, double ymax)
    {
        for (unsigned i = 0; i < _cols; ++i)
        {
            plot(row, i)->setAxisScale(QwtPlot::yLeft, ymin, ymax);
            plot(row, i)->setAxisScale(QwtPlot::yRight, ymin, ymax);
        }
    }
    
    void setColumnAxisRange (unsigned col, double xmin, double xmax)
    {
        for (unsigned i = 0; i < _rows; ++i)
        {
            plot(i, col)->setAxisScale(QwtPlot::xBottom, xmin, xmax);
            plot(i, col)->setAxisScale(QwtPlot::xTop, xmin, xmax);
        }
    }
    
    void
    allignAllPlots (void);

    virtual ~UTMatrixOfQwtPlot ();
};

/**
 *  implementation
 */
#include <qfont.h>
#include <qpen.h>

#include <qwt_scale_widget.h>

#include <cmath>

UTMatrixOfQwtPlot::UTMatrixOfQwtPlot (QWidget * parent, unsigned rows, unsigned cols) :
        QWidget(parent),
        _plot(cols * rows),
        _rows(rows),
        _cols(cols)
{
    _layout = new QGridLayout(this); ///, (int)rows, (int)cols);
    for (unsigned r = 0; r < rows; ++r)
    {
        for (unsigned c = 0; c < cols; ++c)
        {
            _plot[r * cols + c] = new QwtPlot(this);
            _layout->addWidget(_plot[r * cols + c], r, c);
        }
    }
}

QwtPlot *
UTMatrixOfQwtPlot::plot (unsigned row, unsigned col)
{
    return _plot[row * _cols + col];
}

void UTMatrixOfQwtPlot::allignAllPlots (void)
{
    
    for (unsigned i = 0; i < _plot.size(); ++i)
    {
        _plot[i]->replot();
    }
    
    for (unsigned i = 0; i < _plot.size(); ++i)
    {
        _plot[i]->plotLayout()->setAlignCanvasToScales(false);
    }
    
    for (unsigned r = 0; r < _rows - 1; ++r)
    {
        for (unsigned c = 0; c < _cols; ++c)
        {
            QwtScaleDraw *sd = plot(r, c)->axisScaleDraw(QwtPlot::xBottom);
            sd->enableComponent(QwtAbstractScaleDraw::Backbone, false);
            sd->enableComponent(QwtAbstractScaleDraw::Ticks, false);
            sd->enableComponent(QwtAbstractScaleDraw::Labels, false);
            //plot(r,c)->enableAxis(QwtPlot::xBottom,false);
        }
    }
    for (unsigned c = 0; c < _cols; ++c)
    {
        plot(_rows - 1, c)->enableAxis(QwtPlot::xBottom, true);
        plot(_rows - 1, c)->axisScaleDraw(QwtPlot::xBottom)->enableComponent(
                QwtAbstractScaleDraw::Labels, false);
    }
    
    for (unsigned r = 0; r < _rows; ++r)
    {
        for (unsigned c = 1; c < _cols; ++c)
        {
            QwtScaleDraw *sd = plot(r, c)->axisScaleDraw(QwtPlot::yLeft);
            sd->enableComponent(QwtAbstractScaleDraw::Backbone, false);
            sd->enableComponent(QwtAbstractScaleDraw::Ticks, false);
            sd->enableComponent(QwtAbstractScaleDraw::Labels, false);
            //plot(r,c)->enableAxis(QwtPlot::yLeft,false);
        }
    }
    for (unsigned r = 0; r < _rows; ++r)
    {
        plot(r, 0)->enableAxis(QwtPlot::yLeft, true);
    }
    
    int width = 0;
    for (unsigned r = 0; r < _rows; ++r)
    {
        QPen pen;
        pen.setWidth(plot(r, 0)->axisWidget(QwtPlot::yLeft)->penWidth());
        QFont font = plot(r, 0)->axisWidget(QwtPlot::yLeft)->font();
        width = std::max(width, plot(r, 0)->axisScaleDraw(QwtPlot::yLeft)->extent(pen, font));
    }
    width += (int) (0.1 * width);
    for (unsigned r = 0; r < _rows; ++r)
    {
        plot(r, 0)->axisScaleDraw(QwtPlot::yLeft)->setMinimumExtent(width);
    }
    
    for (unsigned i = 0; i < _plot.size(); ++i)
    {
        int margin = 30;
        _plot[i]->plotLayout()->setCanvasMargin(margin, QwtPlot::yLeft);
        _plot[i]->plotLayout()->setCanvasMargin(margin, QwtPlot::yRight);
    }
}

UTMatrixOfQwtPlot::~UTMatrixOfQwtPlot ()
{
}

#endif /*UTMATRIXOFQWTPLOT_H_*/
