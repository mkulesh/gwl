echo
echo --------------------------------------------------
echo -- Generation of two modes synthetic seismogram --
echo --------------------------------------------------
..\..\..\bin\gwlCreateAxis --outfile=SynthSigBTime.dat --count=1024 --min=0 --max=1.023 --name="Time" --nomess

..\..\..\bin\gwlCreateAxis --outfile=SynthSigBFreq.dat --count=1024 --min=0.1 --max=100 --name="Frequency" --nomess

..\..\..\bin\gwlDispModel --infile=SynthSigBFreq.dat --outfile=SynthSigBMod1.dat --analyt --wn=vel --wnpar=1100,400,20 --atn=polin --atnpar=6.204e-04,3.698e-05,-3.200e-06,6.000e-08,-3.336e-10 --name="1th dispersion model" --nomess

..\..\..\bin\gwlSignalGen --infile=SynthSigBTime.dat --outfile=SynthSigBSour1.dat --type=rickdiss --par=120,0,150,1100,400,20 --name="Propagated Ricker wavelet" --nomess

..\..\..\bin\gwlDiffeoDisp --infile=SynthSigBSour1.dat --outfile=SynthSigBSig1.dat --model=SynthSigBMod1.dat --step=4 --dist=170 --name="1th propagation mode"

..\..\..\bin\gwlDispModel --infile=SynthSigBFreq.dat --outfile=SynthSigBMod2.dat --analyt --wn=vel --wnpar=1200,800,50 --atn=polin --atnpar=0.00120,3.898e-005,-3.16e-06,6.03e-08,-3.336e-10 --name="2nd dispersion model" --nomess

..\..\..\bin\gwlSignalGen --infile=SynthSigBTime.dat --outfile=SynthSigBSour2.dat --type=rickdiss --par=70,60,150,1200,800,50 --name="Propagated Ricker wavelet" --nomess

..\..\..\bin\gwlDiffeoDisp --infile=SynthSigBSour2.dat --outfile=SynthSigBSig2.dat --model=SynthSigBMod2.dat --step=4 --dist=170 --name="2nd propagation mode"

..\..\..\bin\gwlSignalSum --infile=SynthSigBSig1.dat,SynthSigBSig2.dat --outfile=SynthSigBSig.dat --name="Synthetic seismogram with two modes" --nomess

echo
echo --------------------------------------------------
echo ------- Two component modulus optimization -------
echo --------------------------------------------------
..\..\..\bin\gwlCreateAxis --outfile=SynthSigBFreq.dat --count=128 --min=5 --max=80 --name="Frequency" --nomess

..\..\..\bin\gwlCwt --infile=SynthSigBSig.dat --outfile=SynthSigBCwt.dat --wttype=1 --freq=SynthSigBFreq.dat --wavelet=morlet --wavpar=3 --name="Wavelet spectrum"

..\..\..\bin\gwlDispModel --infile=SynthSigBFreq.dat --outfile=SynthSigBinit.dat --analyt --wn=polin --wnpar=0.0,+4.6633e-04,+5.5334e-05,-1.9765e-06,+2.7787e-08,-1.3076e-10 --atn=polin --atnpar=1.0133e-03,-8.8029e-06,3.305e-07,3E-10 --name="initial dispersion model"

..\..\..\bin\gwlOptiSP --infile=SynthSigBinit.dat --outfile=SynthSigBOpt1.dat --spec=SynthSigBCwt.dat --ospec=SynthSigBOpt1sp.dat --dist=170 --cmpl=3 --eps=1E-3 --chan=1,2 --prog --name="1st optimized parameters"

echo
echo --------------------------------------------------
echo ------- Multicomponent signal optimization -------
echo --------------------------------------------------
..\..\..\bin\gwlOptiSI --infile=SynthSigBOpt1.dat --outfile=SynthSigBOpt2.dat --sig=SynthSigBSig.dat --osig=SynthSigBOpt2sig.dat --dist=170 --eps=1E-8 --prog --name="2nd optimized parameters"

echo
echo --------------------------------------------------
echo - Multicomponent spectrum optimization: 1st mode -
echo --------------------------------------------------
..\..\..\bin\gwlAutoCorr --infile=SynthSigBSig.dat --outfile=SynthSigBSigCorr.dat --chan1=0,0,0,0,0 --chan2=0,1,2,3,4 --name="Cross correlated seismogram with two modes"

..\..\..\bin\gwlCreateAxis --outfile=SynthSigBFreq.dat --count=128 --min=5 --max=40 --name="First mode frequency" --nomess

..\..\..\bin\gwlCwt --infile=SynthSigBSigCorr.dat --outfile=SynthSigBOpt3corr.dat --wttype=1 --freq=SynthSigBFreq.dat --wavelet=morlet --wavpar=5 --name="Cross correlation wavelet spectrum"

..\..\..\bin\gwlDispModel --infile=SynthSigBFreq.dat --outfile=SynthSigBinit.dat --analyt --wn=polin --wnpar=0.0,+7.6016e-04,-2.5272e-05,+4.1237e-06,-1.4697e-07,+1.5652e-09 --atn=polin --atnpar=+6.567549e-04,+3.079481e-05,-2.470940e-06,+3.199561e-08 --name="initial dispersion model"

..\..\..\bin\gwlOptiSP --infile=SynthSigBinit.dat --outfile=SynthSigBOpt3.dat --spec=SynthSigBOpt3corr.dat --ospec=SynthSigBOpt3sp.dat --dist=170 --cmpl=4 --eps=1E-3 --acorr --prog --name="3rd optimized parameters"

echo
echo --------------------------------------------------
echo - Multicomponent spectrum optimization: 2st mode -
echo --------------------------------------------------
..\..\..\bin\gwlCreateAxis --outfile=SynthSigBFreq.dat --count=128 --min=45 --max=80 --name="First mode frequency" --nomess

..\..\..\bin\gwlCwt --infile=SynthSigBSigCorr.dat --outfile=SynthSigBOpt4corr.dat --wttype=1 --freq=SynthSigBFreq.dat --wavelet=morlet --wavpar=5 --name="Cross correlation wavelet spectrum"

..\..\..\bin\gwlDispModel --infile=SynthSigBFreq.dat --outfile=SynthSigBinit.dat --analyt --wn=polin --wnpar=0.026,-1.4864e-03,+1.1361e-04,-1.9499e-06,+1.5430e-08,-4.9358e-11 --atn=polin --atnpar=-5.0946e-03,+3.1156e-04,-1.02684e-05,+8.2959e-08 --name="initial dispersion model"

..\..\..\bin\gwlOptiSP --infile=SynthSigBinit.dat --outfile=SynthSigBOpt4.dat --spec=SynthSigBOpt4corr.dat --ospec=SynthSigBOpt4sp.dat --dist=170 --cmpl=4 --eps=1E-3 --acorr --prog --name="4th optimized parameters"

del SynthSigBFreq.dat
del SynthSigBinit.dat
del SynthSigBSig1.dat
del SynthSigBSig2.dat
del SynthSigBSour1.dat
del SynthSigBSour2.dat
del SynthSigBTime.dat
