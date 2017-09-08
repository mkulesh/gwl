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

typedef gwlMainObject<PPPSignalContainer<double>, PPPSignalContainer<PPPcomplex> > gwlMainDouble;
typedef gwlMainObject<PPPSignalContainer<PPPcomplex>, PPPSignalContainer<PPPcomplex> > gwlMainCmpl;

template<class AType, class ATypeMain> class gwlMain : public ATypeMain
{
private:

    UTOption_str o_name;
    UTOption_lit o_shift;

public:

    gwlMain (const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut) :
            ATypeMain(aAppName, aModName, aDefIn, aDefOut),
            o_name("n", "name", "<str>",
                   "name of the output signal (by default 'Fourier spectrum')", "Fourier spectrum"),
            o_shift("s", "shift",
                    "if set, shift Fourier spectrum to the regressive/progressive form", false)
    {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_shift);
    }
    
    void calc (void)
    {
        PPPTransFour<AType> trans;
        trans.setShowProgress(false);
        trans.FT(ATypeMain::aDest, ATypeMain::aSource);
        PPPTransFour<PPPcomplex> transcmpl;
        if (o_shift.isOptionGiven()) transcmpl.FTShift(ATypeMain::aDest);
        ATypeMain::aDest.setObjectName(o_name.getValue());
    }
    
};
// end of object

main(int argc, char **argv)
{   
    gwlMain<double,gwlMainDouble> WT1("Fourier transform","gwlCft","signal.dat","fourier.dat");
    WT1.parse(argc, argv);
    ConApplication.onMessage(ConApplication.getAppName());
    PPPObjectIO::ObjectType head = WT1.read_binheader();
    if(head == PPPObjectIO::SIGD)
    {   
        WT1.evaluate();
    }
    else
    {   
        gwlMain<PPPcomplex,gwlMainCmpl> WT2("Fourier transform","gwlCft","signal.dat","fourier.dat");
        WT2.parse(argc, argv);
        WT2.evaluate();
    }
    return 0;
}

