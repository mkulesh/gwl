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
#include "gwlMainObject.h"

typedef gwlMainObject<PPPSpectrContainer<PPPcomplex>, PPPSignalContainer<double> > gwlMainDouble;
typedef gwlMainObject<PPPSpectrContainer<PPPcomplex>, PPPSignalContainer<PPPcomplex> > gwlMainCmpl;

template<class AType, class ATypeMain> class gwlMain : public ATypeMain
{
private:

    UTOption_str o_name;
    UTOption_lit o_ampl;
    UTOption_dbl o_cutoff;

public:
    
    gwlMain (const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut) :
            ATypeMain(aAppName, aModName, aDefIn, aDefOut),
            o_name("n", "name", "<str>", "name of the inverse signal (by default 'inverse signal')",
                   "inverse signal"),
            o_ampl("a", "ampl", "if set, calculate inverse amplitude constant", false),
            o_cutoff("u", "cutoff", "<real>",
                     "wavelet cutoff for slow inverse wavelet transform (by default 0.01)", 0.01)
    {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_progr);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(ATypeMain::o_iscmpl);
        ATypeMain::o_parser.add(ATypeMain::o_wavelet);
        ATypeMain::o_parser.add(ATypeMain::o_wavpar);
        ATypeMain::o_parser.add(o_ampl);
        ATypeMain::o_parser.add(o_cutoff);
    }
    
    void calc (void)
    {
        // prepare transform parameters
        ATypeMain::aSpectrumPar.setInverseParams(ATypeMain::o_wavelet.getValue(),
                                                 ATypeMain::o_wavpar.getValue(),
                                                 o_cutoff.getValue());
        // calculation of inverse signal
        PPPTransWavelet<AType> trans;
        trans.setShowProgress(ATypeMain::o_progr.getValue());
        trans.IWT(ATypeMain::aDest, ATypeMain::aSource, ATypeMain::aSpectrumPar,
                  o_ampl.isOptionGiven());
        ATypeMain::aDest.setObjectName(o_name.getValue());
        strstream str2;
        str2 << ATypeMain::aDest.getInfo() << ends;
        PPPBaseObject::onMessage(str2.str());
    }
    
};
// end of object

main(int argc, char **argv)
{   
    gwlMain<double,gwlMainDouble> WT1("Inverse wavelet transform","gwlIwt","spectr.dat","signal.dat");
    WT1.parse(argc, argv);
    ConApplication.onMessage(ConApplication.getAppName());
    if(WT1.isComplex())
    {   
        gwlMain<PPPcomplex,gwlMainCmpl> WT2("Inverse wavelet transform","gwlIwt","spectr.dat","signal.dat");
        WT2.parse(argc, argv);
        WT2.evaluate();
    }
    else WT1.evaluate();
    return 0;
}

