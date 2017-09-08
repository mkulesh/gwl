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

class gwlMain : public gwlMainObject<PPPAxis, PPPDispersionModel>
{
private:

    UTOption_str o_name;
    UTOption_str o_wn;
    UTOption_str o_wnpar;
    UTOption_str o_atn;
    UTOption_str o_atnpar;
    UTOption_lit o_analyt;

public:
    
    gwlMain (const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut) :
            gwlMainObject<PPPAxis, PPPDispersionModel>(aAppName, aModName, aDefIn, aDefOut),
            o_name("n", "name", "<str>", "name of the model (by default 'dispersion model')",
                   "dispersion model"),
            o_wn("w",
                 "wn",
                 "<str>",
                 "type of the wave number approximation: vel, gauss, polin, bspline, colecole, twogauss (by default 'gauss')",
                 "gauss"),
            o_wnpar("1", "wnpar", "<str>",
                    "parameters of the wave number approximation (by default '')", ""),
            o_atn("a",
                  "atn",
                  "<str>",
                  "type of the attenuation approximation: gauss, polin, bspline (by default 'polin')",
                  "polin"),
            o_atnpar("2", "atnpar", "<str>",
                     "parameters of the attenuation approximation (by default '0')", "0"),
            o_analyt("y", "analyt", "if set, the dispersion model will be analytical", false)
    {
        o_parser.add(o_nomess);
        o_parser.add(o_infile);
        o_parser.add(o_outfile);
        o_parser.add(o_outtype);
        o_parser.add(o_name);
        o_parser.add(o_wn);
        o_parser.add(o_wnpar);
        o_parser.add(o_atn);
        o_parser.add(o_atnpar);
        o_parser.add(o_analyt);
    }
    
    PPPApproximate::ApprType parseApprType (const char *aApprType)
    {
        map<string, PPPApproximate::ApprType> typesMap;
        typesMap["vel"] = PPPApproximate::APTvel;
        typesMap["gauss"] = PPPApproximate::APTgauss;
        typesMap["polin"] = PPPApproximate::APTpolin;
        typesMap["bspline"] = PPPApproximate::APTbspline;
        typesMap["colecole"] = PPPApproximate::APTcolecole;
        typesMap["twogauss"] = PPPApproximate::APTtwogauss;
        if (typesMap.find(aApprType) == typesMap.end()) PPPBaseObject::onError(
        PPPDISPERSIONMODEL_ERRTYPE + string("calc()"));
        return typesMap[aApprType];
    }
    
    void calc (void)
    {
        PPPVectorContainer<double> wnpar, atnpar;
        wnpar.strToVector("{" + string(o_wnpar.getValue()) + "}");
        atnpar.strToVector("{" + string(o_atnpar.getValue()) + "}");
        aDest.setObjectName(o_name.getValue());
        aDest.setAnalytical(true);
        aDest.prepare(parseApprType(o_wn.getValue()), wnpar, parseApprType(o_atn.getValue()),
                      atnpar);
        aDest.evalModel(aSource);
        aDest.setAnalytical(o_analyt.isOptionGiven());
        strstream str;
        str << aDest.getInfo() << ends;
        onMessage(str.str());
    }
    
};
// end of object

main(int argc, char **argv)
{   
    gwlMain AX("Create of a dispersion model","gwlDispModel","frequency.dat","dispmodel.dat");
    AX.parse(argc,argv);
    ConApplication.onMessage(ConApplication.getAppName());
    AX.evaluate();
    return 0;
}

