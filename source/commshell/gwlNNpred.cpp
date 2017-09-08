/* 
 * This file is a part of GWL - Geophysical Wavelet Library
 * Copyright (C) 2007 Mikhail Kulesh and Matthias Holschneider
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * For more information please visit: http://users.math.uni-potsdam.de/~gwl
 * Email: mkulesh@math.uni-potsdam.de
 * ICQ: 103-405-403
 */

#include "gwlMainObject.h"
#include "PPPnnpred.h"

typedef gwlMainObject< PPPVectorContainer<double>, PPPVectorContainer<double> > gwlMainVector;

class gwlMain : public gwlMainVector
  {
  private:
    UTOption_str        o_name;
    UTOption_int        o_slength;
    UTOption_int        o_plength;
    UTOption_int        o_dim;
    UTOption_int        o_step;
    UTOption_int        o_nncount;
    UTOption_dbl        o_pcoeff;
    UTOption_lit        o_l1;
    UTOption_lit        o_l2;
    UTOption_lit        o_l3;
    UTOption_int        o_trend;
    UTOption_dbl        o_sub;

  public:

    gwlMain(const char *aAppName, const char *aModName, const char *aDefIn, const char *aDefOut):
      gwlMainVector(aAppName, aModName, aDefIn, aDefOut),
      o_name   ("n",  "name",    "<str>",  "name of the predicted signal (by default 'prediction')", "prediction"),
      o_slength("s",  "slength", "<unsigned>", "resize source signal to 'slength' points", 0),
      o_plength("p",  "plength", "<unsigned>", "length of prediction (by default 1)", 1),
      o_dim    ("d",  "dim", "<unsigned>", "embedding o_dimension (by default 10)", 10),
      o_step   ("e",  "step", "<unsigned>", "time delay (by default 1)", 1),
      o_nncount("u",  "nncount", "<unsigned>", "count of nearest neighbors (by default 10)", 10),
      o_pcoeff ("f",  "pcoeff", "<real>", "amplitude factor for calculation of farecasting (by default 1) ", 1.0),
      o_l1     ("1",  "l1", "if set, use L1 norm for nearest neighbors searching", false ),
      o_l2     ("2",  "l2", "if set, use L2 norm for nearest neighbors searching", false ),
      o_l3     ("3",  "l3", "if set, use L1-adaptive norm for nearest neighbors searching (by default 'true')", true),
      o_trend  ("0",  "trend", "<int>", "power of polynomial o_trend; if = 0, then without trend (by default 0)", 0),
      o_sub    ("b",  "sub", "<real>", "substituation in the case of zero minimal element (by default 0) ", 0.0)
        {
        o_parser.add(o_nomess);
        o_parser.add(o_infile);
        o_parser.add(o_outfile);
        o_parser.add(o_outtype);
        o_parser.add(o_name);
        o_parser.add(o_slength);
        o_parser.add(o_plength);
        o_parser.add(o_dim);
        o_parser.add(o_step);
        o_parser.add(o_nncount);
        o_parser.add(o_pcoeff);
        o_parser.add(o_l1);
        o_parser.add(o_l2);
        o_parser.add(o_l3);
        o_parser.add(o_trend);
        o_parser.add(o_sub);
        };

    void read(void) {
      aSource.read(o_infile.getValue());
      if(o_slength.isOptionGiven()) aSource.resize(o_slength.getValue());
      };

    void calc(void) {
      // Detrending of source signal
      PPPNearesNeighbors pred;
      if(o_trend.isOptionGiven() && o_trend.getValue() > 0)
        pred.extractTrend(aSource,o_trend.getValue(),true);

      // Nearest neighbors prediction
      pred.prepare(aSource, o_dim.getValue(), o_step.getValue());
      int norm = (o_l1.isOptionGiven())? 1: (o_l2.isOptionGiven())? 2: (o_l3.isOptionGiven())? 3: 0;
      pred.evalPrediction(aDest,o_plength.getValue(),o_nncount.getValue(),norm,o_pcoeff.getValue(),o_sub.getValue());
      aDest.setObjectName(o_name.getValue());

      // Restoring of o_trend
      if(o_trend.isOptionGiven())
        pred.restoreTend(aDest,aSource.size());
      };

  };  // end of object

main(int  argc, char **argv)
  {
  gwlMain AX("Nearest neighbors forecasting of the time series","gwlNNpred","","serie.asc");
  AX.parse(argc,argv);
  ConApplication.onMessage(ConApplication.getAppName());
  AX.evaluate();
  return 0;
  }




