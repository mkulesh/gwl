..\..\..\bin\gwlSignalRead --infile="WorthBasinSyntx.asc,WorthBasinSyntz.asc" --outfile=WorthBasinSyntsig.dat --type=seis2D --smplfreq=500 --to2p --iscmpl --name="2C elastic synthetic seismogram"

..\..\..\bin\gwlET2D --infile=WorthBasinSyntsig.dat --outfile=WorthBasinSyntelli.dat --type=complex --name="2D elliptic properties"

..\..\..\bin\gwlConvert --infile=WorthBasinSyntelli.dat --outfile=WorthBasinSyntratio.dat --comp=4 --name="Reciprocal ellipticity"

..\..\..\bin\gwlConvert --infile=WorthBasinSyntelli.dat --outfile=WorthBasinSynttilt.dat --comp=14 --name="Rise angle"

..\..\..\bin\gwlCreateAxis --outfile=freq.dat --count=64 --min=0.001 --max=50 --sign=full --name="Frequency"

..\..\..\bin\gwlET2DFilter --infile=WorthBasinSyntsig.dat --outfile=WorthBasinSyntfilt.dat --filter=linhor,0.7,0.15,linvert,0.7,0.15,ellihor,0.7,0.15,ellivert,0.7,0.15 --freq=freq.dat --wavelet=cauchy --wavpar=6 --type=complex --name="Filtered seismogram" --prog

del freq.dat
del WorthBasinSyntelli.dat

