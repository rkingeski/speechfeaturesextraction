function [T formi1 formi2 formi3 formi4 band1 band2 band3 band4] = formantes(x,fs,frames,overlap)

%{
%%FORMANTS FUNCTION

x vector of signal;
fs sample frequency;
frames = window legth in miliseconds of the frames
overlap = window legth in miliseconds of the overlap

================================================================================================
Example:
filename =['AUDIOFILE.wav'];
[x,fs]=audioread(filename);
frames = 60;
overlap = 50;

use the function:

[T F1 F2 F3 F4 B1 B2 B3 B4] = formantes(x,fs,60,50);
=============================================================================================

Author: Rafael Kingeski
e-mail: rkingeski@hotmail.com

==============================================================================================

Algorithm based on:

The example of MATLAB help page: https://www.mathworks.com/help/signal/ug/formant-estimation-with-lpc-coefficients.html

And Snell algorithm for LPC analysis:
Snell, Roy C., and Fausto Milinazzo. "Formant location from LPC analysis data." IEEE® Transactions on Speech and Audio Processing. Vol. 1, Number 2, 1993, pp. 129-134.


%}

x=x(:,1);
N = length(x);

nframe = round(frames  * fs / 1000); % converte os ms para amostras
noverlap = round(overlap * fs / 1000);

nfs=round(nframe);
novs=round(noverlap);
 
w = hamming(nframe); %janelamento dos frames usando hamming window;
L = 2^nextpow2(nframe); %tamanho da FFT para cada frame

nlpc=round(fs/1000)+2;


  %% ENERGIA CURTO-TERMO
 
framese = 20; %tempo janelamento em ms

nframee = round(framese  * fs / 1000); % converte os ms para amostras

we = hamming(nframee); %janelamento dos frames usando hamming window;


jan=we.^2;
ener=x.^2;
En=conv(ener(:,1),jan);

winlen=nframee;

temp = (1:length(En))/fs;



  pos = 1; i = 1;
while (pos+nframe < N)
     frame = x(pos:pos+nframe-1);
        cn1= frame(:) .* (w(:));
        [A, err] = lpc(cn1,nlpc);
        
        
        rts = roots(A);

        rts = rts(imag(rts)>=0);
        angz = atan2(imag(rts),real(rts));

        [frqs,indices] = sort(angz.*(fs/(2*pi)));
        bw = -1/2*(fs/(2*pi))*log(abs(rts(indices)));
        
        
        nn = 1;
        for kk = 1:length(frqs)
            if (frqs(kk) > 90 && bw(kk) <400)
                formants(nn) = frqs(kk);
                bands(nn) = bw(kk);
                    if formants(nn)>5000
                        formants(nn)=NaN;
                    end
                    if En(pos)<0.03
                    formants(nn)=NaN;
                    end

                Formi(nn,i)=formants(nn);
                Band(nn,i)= bands(nn);
                    nn = nn+1;

            end
        end

        
       
     pos = pos + (nframe - noverlap);
     i = i + 1;
end

t = (0:N-1)/fs;
T = (round(nframe/2):(nframe-noverlap):N-1-round(nframe/2))/fs;

formi1=Formi(1,:);
formi2=Formi(2,:);
formi3=Formi(3,:);
formi4=Formi(4,:);

band1= Band(1,:);
band2= Band(2,:);
band3= Band(3,:);
band4= Band(4,:);