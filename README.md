# Geophysical Wavelet Library

In this project, we consider and summarize applications of the continuous wavelet transform to 2C and 3C polarization analysis and filtering, modeling the dispersed and attenuated wave propagation in the time–frequency domain, and estimation of the phase and group velocity and the attenuation from a seismogram. Along with a mathematical overview of each of the presented methods, we show that all these algorithms are logically combined into one software package “Geophysical Wavelet Library” developed by the authors. The novelty of this package is that we incorporate the continuous wavelet transform into the library, where the kernel is the time–frequency polarization and dispersion analysis. This library has a wide range of potential applications in the field of signal analysis and may be particularly suitable in geophysical problems that we illustrate by analyzing synthetic, geomagnetic and real seismic data.

We kindly ask you to acknowledge GWL and its authors in any program or publication in which you use GWL. You are not required to do so; 
it is up to your common sense to decide whether you want to comply with this request or not. 

If you want to acknowledge GWL, please cite [following paper](http://www.sciencedirect.com/science/article/pii/S0098300408001568): 
M.Kulesh, M.Holschneider, M.S.Diallo (2008) Geophysical wavelet library: Applications of the continuous wavelet transform to the polarization and dispersion analysis of signals. Computers & Geosciences. Volume 34, Issue 12.

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
* [A graphic extension to the Qt GUI application framework QWT](https://sourceforge.net/projects/qwt/)


