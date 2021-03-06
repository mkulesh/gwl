echo ---------------------------------------------------------
..\..\bin\gwlCreateAxis --scale=lin --count=128 --min=0 --max=127 --outfile=exmDeltaTime.dat --outtype=2 --name="Time" 

echo ---------------------------------------------------------
..\..\bin\gwlCreateAxis --scale=lin --count=128 --min=0.01 --max=1 --outfile=exmDeltaFreq.dat --outtype=2 --name="Frequency"
..\..\bin\gwlConvert --infile=exmDeltaFreq.dat --outfile=exmDeltaFreqASC.dat --outtype=1 --nomess

echo ---------------------------------------------------------
..\..\bin\gwlSignalGen --infile=exmDeltaTime.dat --type=delta --par=64,1 --outfile=exmDeltaSig.dat --outtype=2 --name="Delta function"
..\..\bin\gwlConvert --infile=exmDeltaSig.dat --outfile=exmDeltaSigASC.dat --outtype=1 --nomess

echo 
echo ---------------------------------------------------------
..\..\bin\gwlCwt --infile=exmDeltaSig.dat --wttype=0 --freq=exmDeltaFreq.dat --wavelet=morlet --wavpar=4 --cutoff=0.01 --outfile=exmDeltaCwt.dat --name="Wavelet spectrum"

echo 
echo ---------------------------------------------------------
..\..\bin\gwlConvert --infile=exmDeltaCwt.dat --comp=3 --outfile=exmDeltaCwtMod.dat --outtype=1 --nomess
..\..\bin\gwlConvert --infile=exmDeltaCwt.dat --comp=4 --outfile=exmDeltaCwtArg.dat --outtype=1 --nomess
..\..\bin\gwlConvert --infile=exmDeltaCwt.dat --comp=3,4 --outfile=exmDeltaCwtGnu.dat --outtype=3
..\..\bin\gwlConvert --infile=exmDeltaCwt.dat --comp=3,4 --outfile=exmDeltaCwt.dat --outtype=2

del exmDeltaTime.dat
del exmDeltaFreq.dat
