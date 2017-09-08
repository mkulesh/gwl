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
#define PPPCONF_USEPARSING true      // using command line parsing
#define PPPCONF_CHECKINDEX true      // index control in all container objects
#define PPPCONF_USEFFTW3 true        // using FFTW3 external library for FFT calculation for linux platform
#define PPPCONF_USEQWT true          // include QWT and related graphical objects
#include "PPPtypes.h"

#define _MIN_SIGNAL_HIGHT 80
#define _MAX_SIGNAL_HIGHT 160
#define _MIN_SPECTRUM_HIGHT 120
#define _N_SAMPLE_WAVELET 1024

/************************************************************************
 * gwlMain
 ***********************************************************************/
class gwlMain : public QMainWindow
{
    Q_OBJECT // is needed since we have slots
    
private:
    
    // options
    UTParsing o_parser;
    UTOption_lit o_nomess;
    UTOption_file o_signal;
    UTOption_int o_chanelsignal;
    UTOption_file o_spectrum;
    UTOption_int o_chanelspectrum;

    // the windows menu
    QPopupMenu *_file_menu;
    QPopupMenu *_layout_menu;
    QPopupMenu *_help_menu;

    // calculated binary data
    PPPSignalPlot _signal;			     // the signal
    PPPSpectrContainer<double> _spec;		             // the wavelet spectrum
    PPPSpectrContainer<double> _specplus, _specminus;        // progressive and regressive spectrum
    PPPSpectrParams _spectrumPar;                        // spectrum parameters
    
    //  plot data
    UTMatrixOfQwtPlot * _centralPlot;
    vector<QwtPlotPicker *> _click;
    unsigned _nPlotRow;
    unsigned _nPlotCol;

private slots:

    void slotAbout()
    {   
        QMessageBox::information(this,
                "About gwlPlot",
                "gwlPlot is an interactive plotting tool for GWL library\n"
                "Authors: Matthias Holschneider and Mikhail Kulesh\n"
                "email: matthias.holschneider@gmail.com mikhail.kulesh@gmail.com");
    };

    void slotAboutQt ()
    {
        QMessageBox::aboutQt(this, "About Qt");
    }
    
    void slotShowWavelet (const QwtDoublePoint &p)
    {
        cout << p.x() << " " << p.y() << endl;
    }
    
public:
    
    gwlMain (const char *aAppName, const char *aModName) :
            QMainWindow(),
            o_parser(aModName, aAppName),
            o_nomess("m", "nomess", "if set, no messave will be printed", false),
            o_signal("s", "signal", "<file>", "signal file (by default signal.dat)", "signal.dat"),
            o_chanelsignal("u", "chanelsignal", "<unsigned>",
                           "number of chanel of signal (by default 0)", 0),
            o_spectrum("t", "spectrum", "<file>",
                       "input file with spectrum data (by default spectrum.dat)", "spectrum.dat"),
            o_chanelspectrum("v", "chanelspectrum", "<unsigned>",
                             "number of chanel of spetrum (by default 0)", 0)
    {
        o_parser.add(o_signal);
        o_parser.add(o_chanelsignal);
        o_parser.add(o_spectrum);
        o_parser.add(o_chanelspectrum);
        ConApplication.setAppName(
                string("\n") + string(aAppName) + string(", ") + string(PPPCONF_VERSION));
        
        // set menu bar
        _file_menu = new QPopupMenu();
        _file_menu->insertItem("&Quit", qApp, SLOT(quit()));
        _layout_menu = new QPopupMenu();
        _help_menu = new QPopupMenu();
        _help_menu->insertItem("&About gwlPlot", this, SLOT(slotAbout()));
        _help_menu->insertItem("About Qt", this, SLOT(slotAboutQt()));
        this->menuBar()->insertItem("&File", _file_menu);
        this->menuBar()->insertSeparator();
        this->menuBar()->insertItem("&Help", _help_menu);
        this->statusBar();
    }
    
    void parse (int argc, char **argv)
    {
        o_parser.parse(argc, argv);
        ConApplication.setMessageMode(!o_nomess.getValue());
    }
    
    void evaluate (void)
    {
        // reading wavelet spectrum
        if (o_spectrum.isOptionGiven()) _readSpectrum(_spec, o_spectrum.getValue(), true);
        // reading signal
        if (o_signal.isOptionGiven()) _signal.read(o_signal.getValue());
        this->_updateGeometry();
        this->_updatePlot();
    }
    
private:
    
    template<class AType>
    void _readSpectrum (AType &aData, const char *aName, bool aCheckType = false)
    {
        FILE *infile = fopen(aName, "rb");
        if (infile == NULL) aData.onError(FILE_ERROPEN + string(aName));
        aData.fread(infile);
        if (aCheckType && strcmp(aData.getObjectVer(), PPPSPECTRCONTAINER_OBJVER) == 0) _spectrumPar
                .fread(infile);
        fclose(infile);
    }
    
    void _addSpectrumToPlot (unsigned actualRow, PPPSpectrContainer<double> &aSpec, unsigned aChann)
    {
        UTPlotRasterData::INTERPTYPE finterp;
        finterp =
                (aSpec.getFreq().getType() == PPPAxis::ATlin) ? UTPlotRasterData::LIN :
                                                                UTPlotRasterData::LOG;
        double aFreqMin, aFreqMax;
        if (aSpec.getFreq().getSign() == PPPAxis::ASplus)
        {
            aFreqMin = aSpec.getFreq().getMin();
            aFreqMax = aSpec.getFreq().getMax();
        }
        else
        {
            aFreqMin = fabs(aSpec.getFreq().getMax());
            aFreqMax = fabs(aSpec.getFreq().getMin());
        }
        UTPlotRasterData * tmp = new UTPlotRasterData(aSpec.getChannel(aChann).begin(),  // the data
                aSpec.points(), aSpec.voices(),                             // the dimensions
                aSpec.getTime().getMin(), aSpec.getTime().getMax(),         // time range
                aFreqMin, aFreqMax,                                      // physical frequency range
                UTPlotRasterData::LIN,	                             // linear intepolation for time
                finterp);		                                    // interpolation for freq
        // generate spectrogram data and attach it to the plot
        QwtPlotSpectrogram * sp = new QwtPlotSpectrogram();
        sp->setData(*tmp);
        sp->attach(_centralPlot->plot(actualRow, 0));
        // fix the axis ranges/types
        _centralPlot->setRowAxisRange(actualRow, aFreqMin, aFreqMax);
        // adapt axis to type of spectrum data
        if (finterp == UTPlotRasterData::LOG) _centralPlot->plot(actualRow, 0)->setAxisScaleEngine(
                QwtPlot::yLeft, new QwtLog10ScaleEngine());
        // set row-stretch factor
        _centralPlot->layout()->setRowStretch(actualRow, 10);
        _centralPlot->plot(actualRow, 0)->setMinimumHeight(_MIN_SPECTRUM_HIGHT);
    }
    
    void _addSignalToPlot (unsigned actualRow, PPPSignalPlot &aSig, unsigned aChann)
    {
        vector<QwtPlotCurve *> si;
        double minval, maxval;
        aSig.getPlotCurves(aChann, si, minval, maxval);
        for (unsigned i = 0; i < si.size(); ++i)
            si[i]->attach(_centralPlot->plot(actualRow, 0));
        _centralPlot->setRowAxisRange(actualRow, minval, maxval);
        _centralPlot->setColumnAxisRange(0, _signal.getTime().getMin(), _signal.getTime().getMax());
        _centralPlot->layout()->setRowStretch(actualRow, 1);
        _centralPlot->plot(actualRow, 0)->setMinimumHeight(_MIN_SIGNAL_HIGHT);
        _centralPlot->plot(actualRow, 0)->setMaximumHeight(_MAX_SIGNAL_HIGHT);
    }
    
    void _updateGeometry ()
    {
        _nPlotRow = 0;
        if (o_signal.isOptionGiven()) ++_nPlotRow;
        if (o_spectrum.isOptionGiven())
        {
            ++_nPlotRow;
            if (_spec.getFreq().getSign() == PPPAxis::ASfull) ++_nPlotRow;
        }
        _nPlotCol = 1;  // for the moment only one colum, later more
        _centralPlot = new UTMatrixOfQwtPlot(this, _nPlotRow, _nPlotCol);
        for (unsigned i = 0; i < _nPlotRow; ++i)
        {
            _click.push_back(
                    new QwtPlotPicker(
                            QwtPlot::xBottom,								// working for time axis
                            QwtPlot::yLeft,									// working for freq axis
                            QwtPicker::PointSelection,				// reacting to click and drag
                            QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                            _centralPlot->plot(i, 0)->canvas()));
            _click[i]->setRubberBandPen(QColor(Qt::green));
            _click[i]->setRubberBand(QwtPicker::CrossRubberBand);
            _click[i]->setTrackerPen(QColor(Qt::white));
        }
        this->setCentralWidget(_centralPlot);
    }
    
    void _updatePlot ()
    {
        unsigned actualRow = 0;
        double _tmin = 0, _tmax = 1;
        //  plot the wavelet spectrum
        if (o_spectrum.isOptionGiven())
        {
            _tmin = _spec.getTime().getMin();
            _tmax = _spec.getTime().getMax();
            if (_spec.getFreq().getSign() != PPPAxis::ASfull)
            {
                _addSpectrumToPlot(actualRow++, _spec, o_chanelspectrum.getValue());
            }
            else
            {
                PPPTransWavelet<double> _trans;
                _trans.WaveletSeparate(_specplus, _specminus, _spec);
                _addSpectrumToPlot(actualRow++, _specplus, o_chanelspectrum.getValue());
                _addSpectrumToPlot(actualRow++, _specminus, o_chanelspectrum.getValue());
            }
        }
        //  plot the signal
        if (o_signal.isOptionGiven())
        {
            _tmin = _signal.getTime().getMin();
            _tmax = _signal.getTime().getMax();
            _addSignalToPlot(actualRow, _signal, o_chanelsignal.getValue());
        }
        // finish the matrix of plots
        _centralPlot->setColumnAxisRange(0, _tmin, _tmax);
        _centralPlot->allignAllPlots();
    }
    
};
// end of object

/************************************************************************
 * main()
 ***********************************************************************/
#include "gwlPlot_moc.cpp"

int main (int argc, char * argv[])
{
    QApplication a(argc, argv);
    gwlMain * plotW = new gwlMain("Plotting of the wavelet spectrum", "gwlPlot");
    plotW->parse(argc, argv);
    ConApplication.onMessage(ConApplication.getAppName());
    plotW->evaluate();
    plotW->setGeometry(40, 40, 800, 400);
    a.setMainWidget(plotW);
    plotW->show();
    return a.exec();
}
