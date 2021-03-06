..\..\..\bin\gwlCreateAxis --outfile=TwoNoisedFreq.dat --count=1024 --min=9.7704E-4 --max=2047 --name="Frequency" --nomess

..\..\..\bin\gwlDispModel --infile=TwoNoisedFreq.dat --outfile=TwoNoisedMod1.dat --analyt --wn=vel --wnpar=1256,-62,50 --name="1th dispersion model" --nomess

..\..\..\bin\gwlDispModel --infile=TwoNoisedFreq.dat --outfile=TwoNoisedMod2.dat --analyt --wn=vel --wnpar=1275,125,150 --name="2nd dispersion model" --nomess

..\..\..\bin\gwlCreateAxis --outfile=TwoNoisedTime.dat --count=2048 --min=0 --max=1 --name="Time" --nomess

..\..\..\bin\gwlSignalGen --infile=TwoNoisedTime.dat --outfile=TwoNoisedNoise.dat --type=zero --noise=7 --name="Random noise"

..\..\..\bin\gwlWavelets --infile=TwoNoisedTime.dat --outfile=TwoNoisedSour.dat --outtype=2 --wavelet=cauchy --wavpar=5 --time=0.2 --freq=70

..\..\..\bin\gwlSignalSum --infile=TwoNoisedNoise.dat,TwoNoisedSour.dat --outfile=TwoNoisedSour.dat --name="Noised cauchy wavelet" --nomess

..\..\..\bin\gwlDiffeoDisp --infile=TwoNoisedSour.dat --outfile=TwoNoisedSig1.dat --model=TwoNoisedMod1.dat --step=6 --decrs=2 --dist=113 --name="1th propagation mode" --nomess

..\..\..\bin\gwlDiffeoDisp --infile=TwoNoisedSour.dat --outfile=TwoNoisedSig2.dat --model=TwoNoisedMod2.dat --step=6 --decrs=2 --dist=113 --name="2nd propagation mode" --nomess

..\..\..\bin\gwlSignalSum --infile=TwoNoisedSig1.dat,TwoNoisedSig2.dat --outfile=TwoNoisedSig.dat --name="Synthetic seismogram with two modes" --nomess

..\..\..\bin\gwlCreateAxis --outfile=TwoNoisedFreq.dat --count=256 --min=1 --max=250 --name="Frequency" --nomess

..\..\..\bin\gwlCwt --infile=TwoNoisedSig.dat --outfile=TwoNoisedCwt.dat --wttype=2 --freq=TwoNoisedFreq.dat --wavelet=morlet --wavpar=2 --name="Wavelet spectrum"

..\..\..\bin\gwlCreateAxis --outfile=TwoNoisedVel.dat --count=256 --min=1190 --max=1450 --name="Velocity" --nomess

..\..\..\bin\gwlTransFK --infile=TwoNoisedCwt.dat --outfile=TwoNoisedArg.dat --vel=TwoNoisedVel.dat --inter=spline --corr=arg --filter=0.1 --norm --dist=113 --name="Freauency-velocity image" --prog

..\..\..\bin\gwlTransFK --infile=TwoNoisedCwt.dat --outfile=TwoNoisedCphase.dat --vel=TwoNoisedVel.dat --inter=spline --corr=cphase --filter=0.1 --norm --dist=113 --name="Frequency-velocity image" --prog

del TwoNoisedFreq.dat
del TwoNoisedTime.dat
del TwoNoisedSig1.dat
del TwoNoisedSig2.dat
del TwoNoisedSour.dat
del TwoNoisedVel.dat
del TwoNoisedCwt.dat
del TwoNoisedNoise.dat
