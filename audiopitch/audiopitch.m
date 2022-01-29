function [F0 T] = audiopitch(x,fs,frames,overlap)
%% PITCH FUNCTION
%{
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
=============================================================================================

Author: Rafael Kingeski
e-mail: rkingeski@hotmail.com

==============================================================================================

Algorithm based on:
Cepstrum Pitch Determination
The Journal of the Acoustical Society of America 41, 293 (1967); https://doi.org/10.1121/1.1910339
A. Michael Noll


%}




x=x(:,1); % if the sound is stereo then reduce for mono
N = length(x);

nframe = round(frames  * fs / 1000); %convert time to samples
noverlap = round(overlap * fs / 1000);

%sizes for the cepstral window
t1=floor(fs*0.002); % 2ms 
t2=floor(fs*0.03); % 30ms
 
 
w = hamming(nframe); %hamming window applied in the frames;
L = 2^nextpow2(nframe); %size of FFT

%% energy for threshold
framese = 20; %tempo janelamento em ms

nframee = round(framese  * fs / 1000); % convert ms to samples

we = hamming(nframee); %hamming window for the frames of energy;

jan=we.^2;
ener=x.^2;
En=conv(ener(:,1),jan);

winlen=nframee;

temp = (1:length(En))/fs;
threshold=max(En)/(db2mag(40)); %value of the threshold
 
 
 

%% PITCH 
 
i = 1;
pos=1;
while (pos+nframe < N)
      frame = x(pos:pos+nframe-1);
      cn1= frame(:) .* w(:);
      y = fft(cn1, L);
      cn2 = abs(ifft(log((abs(y)).^2),L));
     %C(:,i) = cn2;
     [~,px]=max(abs(cn2(t1:t2)));
     f0 = fs/(t1+px-1);
     F0(i) = f0;
     
     if En(pos)<threshold
         F0(i)=NaN;
     end

     pos = pos + (nframe - noverlap);
     i = i + 1;
 end
 

T = (round(nframe/2):(nframe-noverlap):N-1-round(nframe/2))/fs;

F0d=diff(F0)/((frames-overlap)/1000);
F02d=diff(F0d)/((frames-overlap)/1000);

Fdd=abs(F0d);
M=size(F0,2);

 for ni=2:M-1
     if Fdd(ni) > 20000
       F0(ni)=NaN;
       F0d(ni)=NaN;
        F02d=NaN;
     end
 end

for nj=2:length(F0)
    if F0(nj)>320
        F0(nj)=NaN;
    end
end
 
F0= smooth(F0,'moving',3);
F0= smooth(F0,'moving',2);  

 
