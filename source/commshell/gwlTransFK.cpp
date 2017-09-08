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

typedef gwlMainObject<PPPSpectrContainer<PPPcomplex>, PPPSpectrContainer<double> > gwlMainDouble;

template<class ATypeMain> class gwlMain : public ATypeMain
{
private:

    UTOption_str o_name;
    UTOption_str o_vel;
    UTOption_str o_inter;
    UTOption_str o_corr;
    UTOption_dbl o_filter;
    UTOption_lit o_norm;
    UTOption_dbl o_dist;

public:
    
    gwlMain (const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut) :
            ATypeMain(aAppName, aModName, aDefIn, aDefOut),
            o_name("n", "name", "<str>", "name of the inverse signal (by default 'f-k image')",
                   "f-k image"),
            o_vel("v", "vel", "<str>", "name of the file with velocity axis (by default 'vel.dat')",
                  "vel.dat"),
            o_inter("e", "inter", "<str>",
                    "type of the velocity interpolation: index, spline (by default 'index')",
                    "index"),
            o_corr("l", "corr", "<str>",
                   "type of component for correlation: abs, arg, cphase (by default 'arg')", "arg"),
            o_filter("f", "filter", "<real>", "parameter of energy filter (by default 0)", 0.0),
            o_norm("g", "norm", "if set, each voice in output f-k image will be normed", false),
            o_dist("d", "dist", "<real>", "propagation distance (by default 0)", 0.0)
    {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_progr);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_vel);
        ATypeMain::o_parser.add(o_inter);
        ATypeMain::o_parser.add(o_corr);
        ATypeMain::o_parser.add(o_filter);
        ATypeMain::o_parser.add(o_norm);
        ATypeMain::o_parser.add(o_dist);
    }
    
    void calc (void)
    {
        map<string, PPPTransFK::InterpType> InterpType;
        InterpType["index"] = PPPTransFK::MMIindex;
        InterpType["spline"] = PPPTransFK::MMIspline;
        if (InterpType.find(o_inter.getValue()) == InterpType.end()) ATypeMain::onError(
        PPPTRANSFK_ERRINT + string("calc()"));
        map<string, PPPTransFK::CorrType> CorrType;
        CorrType["abs"] = PPPTransFK::MMCabs;
        CorrType["arg"] = PPPTransFK::MMCarg;
        CorrType["cphase"] = PPPTransFK::MMCcphase;
        if (CorrType.find(o_corr.getValue()) == CorrType.end()) ATypeMain::onError(
        PPPTRANSFK_ERRCOR + string("calc()"));
        // prepare transform parameters
        PPPTransFK trans;
        trans.setIntType(InterpType[o_inter.getValue()]);
        trans.setCorrType(CorrType[o_corr.getValue()]);
        trans.setFilter(o_filter.getValue());
        trans.setNormVoice(o_norm.isOptionGiven());
        trans.setDistance(o_dist.getValue());
        trans.setShowProgress(ATypeMain::o_progr.getValue());
        // Frequency-wavenumber transformation
        PPPAxis aVelAxis;
        ATypeMain::read_bindata(aVelAxis, o_vel.getValue());
        trans.evalSpectrCorrelation(ATypeMain::aDest, ATypeMain::aSource, aVelAxis);
        ATypeMain::aDest.setObjectName(o_name.getValue());
        // information
        strstream str;
        str << ATypeMain::aDest.getInfo() << ends;
        PPPBaseObject::onMessage(str.str());
    }
    
};
// end of object

main(int argc, char **argv)
{   
    gwlMain<gwlMainDouble> WT1("Frequency-wavenumber transformation","gwlTransFK","spectr.dat","fk.dat");
    WT1.parse(argc, argv);
    ConApplication.onMessage(ConApplication.getAppName());
    WT1.evaluate();
    return 0;
}

