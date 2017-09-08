[![License](https://img.shields.io/badge/license-GNU_GPLv3-orange.svg)](https://github.com/mkulesh/gwl/blob/master/LICENSE)

![Geophysical Wavelet Library](https://github.com/mkulesh/gwl/blob/master/images/title.png)

Geophysical Wavelet Library (GWL) is a software package based on the continuous wavelet transform that allows to perform the direct and inverse continuous wavelet transform, 2C and 3C polarization analysis and filtering, modeling the dispersed and attenuated wave propagation in the time-frequency domain and optimization in signal and wavelet domains with the aim to extract velocities and attenuation parameters from a seismogram. The novelty of this package is that we incorporate the continuous wavelet transform into the library, where the kernel is the time–frequency polarization and dispersion analysis. This library has a wide range of potential applications in the field of signal analysis and may be particularly suitable in geophysical problems that we illustrate by analyzing synthetic, geomagnetic and real seismic data.

We kindly ask you to acknowledge GWL and its authors in any program or publication in which you use GWL. You are not required to do so;  it is up to your common sense to decide whether you want to comply with this request or not. 

If you want to acknowledge GWL, please cite [following paper](http://www.sciencedirect.com/science/article/pii/S0098300408001568): 
M.Kulesh, M.Holschneider, M.S.Diallo (2008) *Geophysical wavelet library: Applications of the continuous wavelet transform to the polarization and dispersion analysis of signals*. Computers & Geosciences. Volume 34, Issue 12.

## GWL structure and implementation technology

GWL includes three logical levels: the library level, the level of command line tools and the interface level:
![GWL Structure](https://github.com/mkulesh/gwl/blob/master/images/library_structure.png)

* The main part of the library level is a C++ hierarchical object library called PPP (Parametric Processing of Pulsations). This library contains an object-based implementation of main data types and mathematical objects and implements all used algorithms. We designed PPP with the assumption that it could be potentially used outside of GWL in any other project related to the time-frequency analysis of the signals. PPP uses the C++ Standard Template Library, an ANSI C command line parser, a C subroutine library for computing the discrete Fourier transform, and Linux GUI components and utility classes (Qwt). 

* GWL command line level is a set of independent C++ modules. These modules are based on the PPP library and provide a command line interface for all methods implemented in this library. After the compilation, we obtain a set of executable modules placed in the GWL/bin directory.

* To perform a calculation using GWL, we run certain modules from this directory in the appointed order. The calculation parameters must be given by command line. The data exchange between different modules is implemented by the data files with binary stream formats. After the calculation process is finished, we collect ASCII or binary files with the calculation results. 

* In general, we can plot these calculation results, saved in ASCII format, using any plotting software. To reduce the required hard drive space and improve plotting performance, we can also store the results as binary files. Toward this end, we developed a special MATLAB package placed on the interface level of GWL in the directory GWL/mshell. This package allows us to read all binary formats supported in GWL and plot GWL objects, such as a multi-channel signal or a multi-channel wavelet spectrum, using high-level subroutines based on the standard MATLAB plotting commands. This package is developed as an M-file library and does not need any installation before use. 

* The second possible use of the GWL/mshell library is an integration of the calculation process with plotting subroutines using an M-file program. In this way, we added some procedures into GWL/mshell to directly execute modules from the GWL command line level. 

* To demonstrate the calculation technique using GWL, we stored a lot of examples into GWL/solutions directory. This folder together with GWL MATLAB tools constitutes a part of GWL interface level. We preciously used all these examples to show the application possibilities of the CWT for dispersion and polarization analysis of synthetic, geomagnetic and real seismic data.

## Supported platforms

Currently, GWL has been tested on following platforms:
* gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-4)
* gcc (GCC) 6.4.1 20170727 (Red Hat 6.4.1-1)

## How to build

See Wiki page [How to build and test GWL](https://github.com/mkulesh/gwl/wiki#how-to-build-and-test-gwl)

## License

This software is published under the *GNU General Public License, Version 3*

Copyright (C) 2014-2017 Mikhail Kulesh

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program.

If not, see [www.gnu.org/licenses](https://www.gnu.org/licenses).

## Dependencies

This software depends or includes the following third-party libraries or code fragments:
* [ANSI C command line parser Argtable2](http://argtable.sourceforge.net/)
* [C subroutine library for computing the discrete Fourier transform FFTW](http://www.fftw.org/)
* [A graphic extension to the Qt GUI: application framework QWT](https://sourceforge.net/projects/qwt/)

## History

This project was previuosly hosted on http://users.math.uni-potsdam.de/∼gwl
