echo
echo --------------------------------------------------
echo ------------------ Subsection A ------------------
echo --------------------------------------------------
..\..\..\bin\gwlSignalRead --infile=KerpenData.asc --outfile=KerpenDataAsig.dat --format=ASCII --type=seis --smplfreq=4000 --chan=17,18,19,20 --to2p --name="Subsection A" 
 
..\..\..\bin\gwlAutoCorr --infile=KerpenDataAsig.dat --outfile=KerpenDataAcorr.dat --chan1=0,0,0,0 --chan2=0,1,2,3 --name="Cross correlation of subsection A"

..\..\..\bin\gwlCreateAxis --outfile=KerpenDataFreq.dat --count=128 --min=25 --max=45 --name="Frequency" --nomess

..\..\..\bin\gwlCwt --infile=KerpenDataAcorr.dat --outfile=KerpenDataAcwt.dat --wttype=1 --freq=KerpenDataFreq.dat --wavelet=morlet --wavpar=3 --name="Wavelet spectrum of subsection A"

..\..\..\bin\gwlDispModel --infile=KerpenDataFreq.dat --outfile=KerpenDataInit.dat --analyt --wn=polin --wnpar=0.0,-1.6123e-02,+3.7980e-03,-1.6939e-04,+2.2788e-06 --atn=polin --atnpar=-1.1466e+00,+1.2960e-01,-9.4122e-03,+1.6986e-04 --name="initial dispersion model of subsection A" --nomess

..\..\..\bin\gwlOptiSP --infile=KerpenDataInit.dat --outfile=KerpenDataAoptpar.dat --spec=KerpenDataAcwt.dat --ospec=KerpenDataAcwt.dat --dist=2 --cmpl=4 --eps=1E-4 --acorr --prog --name="Optimized parameters of subsection A"

..\..\..\bin\gwlIwt --infile=KerpenDataAcwt.dat --outfile=KerpenDataAoptcorr.dat --wavelet=delta --name="Optimazed cross correlation of subsection A"

del KerpenDataAcwt.dat
del KerpenDataAsig.dat

echo
echo --------------------------------------------------
echo ---------- Subsection B: multi signals -----------
echo --------------------------------------------------
..\..\..\bin\gwlSignalRead --infile=KerpenData.asc --outfile=KerpenDataBsig.dat --format=ASCII --type=seis --smplfreq=4000 --chan=23,24,25,26,27 --to2p --name="Subsection B" 
 
..\..\..\bin\gwlAutoCorr --infile=KerpenDataBsig.dat --outfile=KerpenDataBcorr.dat --chan1=0,0,0,0,0 --chan2=0,1,2,3,4 --name="Cross correlation of subsection B"

..\..\..\bin\gwlCreateAxis --outfile=KerpenDataFreq.dat --count=128 --min=5 --max=100 --name="Frequency" --nomess

..\..\..\bin\gwlCwt --infile=KerpenDataBsig.dat --outfile=KerpenDataBcwt.dat --wttype=1 --freq=KerpenDataFreq.dat --wavelet=morlet --wavpar=3 --name="Wavelet spectrum of subsection B"

..\..\..\bin\gwlDispModel --infile=KerpenDataFreq.dat --outfile=KerpenDataInit.dat --analyt --wn=bspline --wnpar=+1.103e-01,-8.284e-02,+1.299e-01,+1.787e-01,+1.731e-01,+3.063e-01,+5.452e-01 --atn=bspline --atnpar=+1.058e-01,+1.237e-01,-5.394e-02,+3.857e-01,-1.319e-02,+1.078e-01,-3.512e-02 --name="initial dispersion model of subsection B" --nomess

..\..\..\bin\gwlOptiSI --infile=KerpenDataInit.dat --outfile=KerpenDataB1optpar.dat --sig=KerpenDataBsig.dat --osig=KerpenDataB1optsig.dat --dist=2 --eps=1E-8 --prog --name="Optimized parameters of subsection B"

echo
echo --------------------------------------------------
echo ----------- Subsection B: first mode -------------
echo --------------------------------------------------
..\..\..\bin\gwlCreateAxis --outfile=KerpenDataFreq.dat --count=128 --min=20 --max=30 --name="Frequency" --nomess

..\..\..\bin\gwlCwt --infile=KerpenDataBcorr.dat --outfile=KerpenDataB2cwt.dat --wttype=1 --freq=KerpenDataFreq.dat --wavelet=morlet --wavpar=2 --name="Wavelet spectrum of subsection B"

..\..\..\bin\gwlDispModel --infile=KerpenDataFreq.dat --outfile=KerpenDataInit.dat --analyt --wn=bspline --wnpar=+4.749e-02,+5.088e-02,+6.264e-02,+7.722e-02,+9.396e-02,+1.121e-01,+1.147e-01 --atn=bspline --atnpar=+8.406e-02,+4.393e-02,+2.207e-02,+4.280e-02,+4.592e-02,+3.121e-02,+1.987e-01 --name="initial dispersion model of subsection B" --nomess

..\..\..\bin\gwlOptiSP --infile=KerpenDataInit.dat --outfile=KerpenDataB2optpar.dat --spec=KerpenDataB2cwt.dat --ospec=KerpenDataB2cwt.dat --dist=2 --cmpl=4 --eps=1E-4 --acorr --prog --name="Optimized parameters of subsection B"

..\..\..\bin\gwlIwt --infile=KerpenDataB2cwt.dat --outfile=KerpenDataB2optcorr.dat --wavelet=delta --name="Optimazed cross correlation of subsection B"

del KerpenDataB2cwt.dat

echo
echo --------------------------------------------------
echo ----------- Subsection B: second mode ------------
echo --------------------------------------------------
..\..\..\bin\gwlCreateAxis --outfile=KerpenDataFreq.dat --count=128 --min=30 --max=45 --name="Frequency" --nomess

..\..\..\bin\gwlCwt --infile=KerpenDataBcorr.dat --outfile=KerpenDataB3cwt.dat --wttype=1 --freq=KerpenDataFreq.dat --wavelet=morlet --wavpar=2 --name="Wavelet spectrum of subsection B"

..\..\..\bin\gwlDispModel --infile=KerpenDataFreq.dat --outfile=KerpenDataInit.dat --analyt --wn=bspline --wnpar=+8.502e-02,+1.047e-01,+1.163e-01,+1.429e-01,+1.691e-01,+2.025e-01,+2.276e-01 --atn=bspline --atnpar=-1.755e-01,+1.056e-01,+9.912e-02,+8.441e-02,+8.757e-02,+1.062e-01,+9.459e-02 --name="initial dispersion model of subsection B" --nomess

..\..\..\bin\gwlOptiSP --infile=KerpenDataInit.dat --outfile=KerpenDataB3optpar.dat --spec=KerpenDataB3cwt.dat --ospec=KerpenDataB3cwt.dat --dist=2 --cmpl=4 --eps=1E-4 --acorr --prog --name="Optimized parameters of subsection B"

..\..\..\bin\gwlIwt --infile=KerpenDataB3cwt.dat --outfile=KerpenDataB3optcorr.dat --wavelet=delta --name="Optimazed cross correlation of subsection B"

del KerpenDataB3cwt.dat
del KerpenDataFreq.dat
del KerpenDataInit.dat
