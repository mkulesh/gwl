..\..\bin\gwlCreateAxis --count=1024 --min=-10 --max=10 --name="Time" --outfile=time.dat --nomess

..\..\bin\gwlWavelets --infile=time.dat --iscmpl --wavelet=morlet --wavpar=0.7 --time=0 --freq=2 --outtype=1 --outfile=waveletMorlet.dat --four

..\..\bin\gwlWavelets --infile=time.dat --wavelet=morletre --wavpar=2 --time=0 --freq=2 --outtype=1 --outfile=waveletReMorlet.dat --four

..\..\bin\gwlWavelets --infile=time.dat --iscmpl --wavelet=cauchy --wavpar=20 --time=0 --freq=2 --outtype=1 --outfile=waveletCauchy.dat --four

..\..\bin\gwlWavelets --infile=time.dat --wavelet=haar --time=0 --freq=0.75971 --outtype=1 --outfile=waveletHaar.dat --four

..\..\bin\gwlWavelets --infile=time.dat --iscmpl --wavelet=shanon --wavpar=0.4 --time=0 --freq=2.5 --outtype=1 --outfile=waveletShanon.dat --four

del time.dat
