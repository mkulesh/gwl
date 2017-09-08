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

#define gwlET2D_ERR01  "type of the elliptic mashine does not correspont to type of input object"
#define gwlET2D_ERR02  "frequency interval for averaging's method is invalid"
#define gwlET2D_ERR03  "type of polarization filter is undefined: "
#define gwlET2D_ERR04  "number of parameter of polarization filter is invalid"

typedef gwlMainObject<PPPSignalContainer<double>, PPPSignalContainer<PPPellipse3D> > gwlMainSigToSig;
typedef gwlMainObject<PPPSpectrContainer<PPPcomplex>, PPPSpectrContainer<PPPellipse3D> > gwlMainSpecToSpec;

template<class ATypeMain> class gwlMain : public ATypeMain
{
private:

    UTOption_str o_name;
    UTOption_str o_type;
    UTOption_int o_tw;
    UTOption_dbl o_filter;
    PPPTransElli aTrans;

public:
    
    gwlMain (const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut) :
            ATypeMain(aAppName, aModName, aDefIn, aDefOut),
            o_name("n", "name", "<str>",
                   "name of the elliptic parameters (by default 'elliptic parameters')",
                   "elliptic parameters"),
            o_type("y", "type", "<str>",
                   "type of the elliptic mashine: morozov, scovar, acovar (by default 'acovar')",
                   "acovar"),
            o_tw("w", "tw", "<unsigned>",
                 "time window length for scovar and acovar methods (by default 2)", 2),
            o_filter("f", "filter", "<real>", "energy filter by spectral transform (by default 0)",
                     0.0)
    {
        ATypeMain::o_parser.add(ATypeMain::o_nomess);
        ATypeMain::o_parser.add(ATypeMain::o_infile);
        ATypeMain::o_parser.add(ATypeMain::o_outfile);
        ATypeMain::o_parser.add(ATypeMain::o_outtype);
        ATypeMain::o_parser.add(o_name);
        ATypeMain::o_parser.add(o_type);
        ATypeMain::o_parser.add(o_tw);
        ATypeMain::o_parser.add(o_filter);
    }
    
    PPPConstETmachine parseMachine ()
    {
        map<string, PPPConstETmachine> typesMap;
        typesMap["morozov"] = WTEMorozov;
        typesMap["scovar"] = WTESumCovar;
        typesMap["acovar"] = WTECovar;
        if (typesMap.find(o_type.getValue()) == typesMap.end()) PPPBaseObject::onError(
        PPPTRANSELLI_ERRMASH + string("calc()"));
        return typesMap[o_type.getValue()];
    }
    
    void calc (void)
    {
        // prepare transform parameters
        aTrans.setMachine(parseMachine());
        aTrans.setShowProgress(false);
        // calculation of elliptic properties
        calc(ATypeMain::aDest, ATypeMain::aSource);
        ATypeMain::aDest.setObjectName(o_name.getValue());
        // information
        strstream str;
        str << ATypeMain::aDest.getInfo() << ends;
        ATypeMain::onMessage(str.str());
    }
    
    void calc (PPPSignalContainer<PPPellipse3D> &aLocDest, PPPSignalContainer<double> &aLocSource)
    {
        if (aTrans.getMachine() == WTEMorozov)
            aTrans.eval3DSignalPar(aLocDest, aLocSource);
        else
            aTrans.eval3DSignalPar(aLocDest, aLocSource, o_tw.getValue());
    }
    
    void calc (PPPSpectrContainer<PPPellipse3D> &aLocDest,
               PPPSpectrContainer<PPPcomplex> &aLocSource)
    {
        if (aTrans.getMachine() == WTEMorozov)
            aTrans.eval3DSpectrPar(aLocDest, aLocSource, o_filter.getValue());
        else
            aTrans.eval3DSpectrPar(aLocDest, aLocSource, o_filter.getValue(), o_tw.getValue());
    }
    
};
// end of object

main(int argc, char **argv)
{   
    string appname = "3D elliptic transform";
    string appcode = "gwlET3D";
    gwlMain<gwlMainSigToSig> WT1(appname.c_str(),appcode.c_str(),"source.dat","ellipar3D.dat");
    WT1.parse(argc, argv);
    ConApplication.onMessage(ConApplication.getAppName());
    PPPObjectIO::ObjectType head = WT1.read_binheader();
    if(head == PPPObjectIO::SIGD)
    {   
        WT1.evaluate();
    }
    else if(head == PPPObjectIO::SPECC)
    {   
        gwlMain<gwlMainSpecToSpec> WT2(appname.c_str(),appcode.c_str(),"source.dat","ellipar3D.dat");
        WT2.parse(argc, argv);
        WT2.evaluate();
    }
    else
    ConApplication.onError(gwlMain_ERR02);

    return 0;
}

