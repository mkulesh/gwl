..\..\..\bin\gwlCreateAxis --outfile=KerpenDataFreq.dat --count=151 --scale=log --min=13.7083 --max=46.1154 --name="Frequency" --nomess

..\..\..\bin\gwlCreateAxis --outfile=KerpenDataVel.dat --count=231 --scale=loge --min=113.9600 --max=330.5790 --name="Velocity" --nomess

..\..\..\bin\gwlSignalRead --infile=KerpenData.asc --outfile=KerpenDataAsig.dat --format=ASCII --type=seis --smplfreq=4000 --chan=16,17,18,19 --name="Subsection A" 

..\..\..\bin\gwlCwt --infile=KerpenDataAsig.dat --outfile=KerpenDataAcwt.dat --wttype=2 --freq=KerpenDataFreq.dat --wavelet=morlet --wavpar=6 --name="Wavelet spectrum"

..\..\..\bin\gwlTransFK --infile=KerpenDataAcwt.dat --outfile=KerpenDataAarg.dat --vel=KerpenDataVel.dat --inter=spline --corr=arg --filter=0.1 --norm --dist=2 --name="Frequency-velocity image" --prog

..\..\..\bin\gwlSignalRead --infile=KerpenData.asc --outfile=KerpenDataBsig.dat --format=ASCII --type=seis --smplfreq=4000 --chan=23,24,25,26,27 --name="Subsection B" 

..\..\..\bin\gwlCwt --infile=KerpenDataBsig.dat --outfile=KerpenDataBcwt.dat --wttype=2 --freq=KerpenDataFreq.dat --wavelet=morlet --wavpar=6 --name="Wavelet spectrum"

..\..\..\bin\gwlTransFK --infile=KerpenDataBcwt.dat --outfile=KerpenDataBarg.dat --vel=KerpenDataVel.dat --inter=spline --corr=arg --filter=0.1 --norm --dist=2 --name="Frequency-velocity image" --prog

del KerpenDataFreq.dat
del KerpenDataVel.dat
del KerpenDataAcwt.dat
del KerpenDataBcwt.dat
