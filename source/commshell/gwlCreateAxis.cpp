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

class gwlMain : public gwlMainObject< PPPAxis, PPPAxis >
  {
  private:
    UTOption_str        o_name;
    UTOption_int        o_count;
    UTOption_str        o_scale;
    UTOption_str        o_sign;
    UTOption_dbl        o_amin;
    UTOption_dbl        o_amax;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      gwlMainObject< PPPAxis, PPPAxis >(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n",  "name",  "<str>", "name of the axis (by default 'axis')", "axis"),
      o_count  ("u",  "count", "<unsigned>", "count of the axis points (by default 1024)", 1024),
      o_scale  ("s",  "scale", "<str>", "scale of the axis: 'lin'-linear, 'log'-logarithmic (by default 'lin')", "lin"),
      o_sign   ("g",  "sign",  "<str>",  "symmetry of the axis: 'prog'-progressive, 'reg'-regressive, 'full'-full (by default 'prog')", "prog"),
      o_amin   ("1",  "min",   "<real>", "minimum value of axis (by default 0.0)", 0.0),
      o_amax   ("2",  "max",   "<real>", "maximum value of axis (by default 1.0)", 1.0)
        {
        o_parser.add(o_nomess);
        o_parser.add(o_outfile);
        o_parser.add(o_outtype);
        o_parser.add(o_name);
        o_parser.add(o_count);
        o_parser.add(o_scale);
        o_parser.add(o_sign);
        o_parser.add(o_amin);
        o_parser.add(o_amax);
        };

    void calc(void) {
      map < string, PPPAxis::AxisType > scaleMap;
        scaleMap["lin"] = PPPAxis::ATlin;
        scaleMap["log"] = PPPAxis::ATlog;
        scaleMap["loge"] = PPPAxis::ATloge;
      if(scaleMap.find(o_scale.getValue()) == scaleMap.end())
        PPPBaseObject :: onError("scale of the axis is incorrect");

      map < string, PPPAxis::AxisSign > signMap;
        signMap["prog"] = PPPAxis::ASplus;
        signMap["reg"] = PPPAxis::ASminus;
        signMap["full"] = PPPAxis::ASfull;
      if(signMap.find(o_sign.getValue()) == signMap.end())
        PPPBaseObject :: onError("symmetry of the axis is incorrect");

      aDest.setObjectName(o_name.getValue());
      aDest.realloc(o_count.getValue());
      aDest.fill(o_amin.getValue(), o_amax.getValue(), scaleMap[o_scale.getValue()]);
      if(signMap[o_sign.getValue()] != PPPAxis::ASplus) aDest.fillsign(signMap[o_sign.getValue()]);

      strstream str;
      str << aDest.getInfo() << ends;
      onMessage(str.str());
      };

  };  // end of object

main(int  argc, char **argv)
  {
  gwlMain AX("Create of an axes","gwlCreateAxis","","axis.dat");
  AX.parse(argc,argv);
  ConApplication.onMessage(ConApplication.getAppName());
  AX.evaluate();
  return 0;
  }




