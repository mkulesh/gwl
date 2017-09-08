..\..\..\bin\gwlSignalRead --infile="WorthBasinRealx.asc,WorthBasinRealz.asc" --outfile=WorthBasinRealsig.dat --type=seis2D --smplfreq=500 --to2p --iscmpl --name="2C elastic real seismogram"

..\..\..\bin\gwlCreateAxis --outfile=freq.dat --sign=full --count=64 --min=0.1 --max=100 --name="Frequency"

..\..\..\bin\gwlET2DFilter --infile=WorthBasinRealsig.dat --outfile=WorthBasinRealfilt.dat --filter=linhor,1.57,0.15,ellihor,1.57,0.15 --freq=freq.dat --wavelet=cauchy --wavpar=6 --type=complex --name="Filtered seismogram" --prog

del freq.dat
