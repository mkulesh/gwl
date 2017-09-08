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
#ifndef _UTPARSING
#define _UTPARSING

#include <functional>
#include <algorithm>
#include "UToptions.h"

/************************************************************************
 * container class for options
 * UTParsing
 ***********************************************************************/
class UTParsing {

  private:
    vector<UTOption *> _data;
    vector<UTMultiOption * > _dataMult;
    void *              _argtable[100];
    char const *        _progname;
    char const *        _shortdescription;

  public:

    UTParsing (char const * progname, char const * shortdescription):
      _progname(progname),
      _shortdescription(shortdescription) {
    };

    void add (UTOption & newOption) {
      _data.push_back(& newOption);
      return;
    }
    
    void add ( 
    	UTMultiOption & newOption
    ) {
    	_dataMult.push_back ( &newOption );
    	for ( unsigned i=0; i<=newOption.getN();++i) {
    		add ( * newOption.getOption ( i ) );	
    	}
    }

	inline void parse (
    	int argc, char * argv[], bool checkErrors = false){
      /* help menu */
      arg_lit * help = arg_lit0 ( "h", "help", "prints this menu and exits");
      struct arg_end * end = arg_end(40); // TODO please check ????
      fill( &_argtable[0], &_argtable[99], (void*)NULL);
      for (unsigned i=0; i<_data.size(); ++i) {
          _argtable[i]=(void*) _data[i]->getData();
       }
      _argtable[_data.size()]= (void*) help;
      _argtable[_data.size()+1]= end;
      if (arg_nullcheck(_argtable) != 0) {
          /* NULL entries were detected, some allocations must have failed */
          arg_free(_argtable);
          cout << MEM_ERRALLOC + string("parse") << endl;
          exit(-1);
          }
      /* Parse the command line as defined by _argtable[] */
      int nerrors = arg_parse(argc,argv,_argtable);
      if ( help -> count > 0 ) {
          printHelp();
          arg_free(_argtable);
          exit (0);
          }
      /* If the parser returned any errors then display them and exit */
	if(argc == 1) nerrors++;
    
	for ( unsigned i=0; i < _dataMult.size(); ++i ) {
		unsigned n = _dataMult [ i ] -> getN();
      	if ( _dataMult [ i ] -> getOption ( 0 ) -> isOptionGiven() ) {
      		for ( unsigned k=0; k < n; ++k ) {
      			if ( ! _dataMult[i] -> getOption(k) -> isOptionGiven() ) {
      				_dataMult [ i ] -> getOption ( k ) -> copyValueFrom ( 
      					_dataMult [ i ] -> getOption ( n )
      				);
      			}
      		}
      	}	
	}
	
    if (checkErrors && (nerrors > 0)) {
          /* Display the error details contained in the arg_end struct.*/
          arg_print_errors(stdout,end,_progname);
          cout << "Try '" << _progname << " --help' for more information" << endl;
          arg_free(_argtable);
          exit (1);
          }
      for_each ( _data.begin(), _data.end(), mem_fun(&UTOption::setParsed));
      return;
    };

    void printHelp(void) {
      cout << "Usage: " << _progname << " ";
      arg_print_syntax(stdout,_argtable,"\n\n");
      cout << _shortdescription << endl;
      arg_print_glossary(stdout,_argtable,"  %-40s %s\n");
    };

    void info(void) {
      for_each(_data.begin(), _data.end(), mem_fun(&UTOption::print));
      return;
    };

}; // end of object

#endif

