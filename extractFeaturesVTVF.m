function [signal dataLF dataMF features] = extractFeaturesVTVF(signal,fs)
%Author: F. Plešinger, ISI of the CAS, Czechia, 2018
%extracts features from 3 seconds ECG signal. It extracts more features, but only 5 features were used in final model.


signal=detrend(signal);
signal = signal';
signal = hfilta(signal, fs,1,35);


%normalize

s2 = signal-min(signal(:));
s3 = s2./max(s2(:));
signal=s3;


%%% end norm
                

dataLF = hobalka(signal,fs,1,8);
dataMF = hobalka(signal,fs,5,25);

med = median(signal);




maxLFMF = 0;

startBlock=1;

%maxLFratio

for i=1:length(signal)
    
    if dataLF(i)<dataMF(i)
        lng = i-startBlock;
        
        if lng>maxLFMF
            maxLFMF = lng;
        end
        
        startBlock = i;
    end
end

endLNG = length(signal)-startBlock;

if (endLNG>maxLFMF)
    maxLFMF = endLNG;
end

features.maxLFMFsecs = maxLFMF/fs;

%positive crossings of signal and its median
pcross= [];
ncross = [];
for i=2:length(signal)
    
    if signal(i-1)<med && signal(i)>=med
        pcross(length(pcross)+1)=i;
    end
    
    if signal(i-1)>=med && signal(i)<med
        ncross(length(ncross)+1)=i;
    end
    
    
end

features.numPMC = length(pcross);
features.numNMC = length(ncross);


if doPicture

figure;
hold on

plot(signal);
plot(dataLF);
plot(dataMF);
plot(pcross,signal(pcross),'r*');
plot(ncross,signal(ncross),'go');

line([0 750],[med med]);

end

%--pomìr sum nad a pod medianem

signalMMed = signal-med;

sumSMM = sum(find(signalMMed>0));
sumSMM2 = sum(find(signalMMed<0));

features.ratioSums = sumSMM/sumSMM2;
if isinf(features.ratioSums)
    features.ratioSums = -1;
end


diffPMC = diff(pcross)./fs;
diffNMC = diff(ncross)./fs;

features.stdDiffPMC = std(diffPMC);
features.stdDiffNMC = std(diffNMC);

features.ratioSTDsPNMC = features.stdDiffPMC/features.stdDiffNMC;

if isinf(features.ratioSTDsPNMC)
    features.ratioSTDsPNMC=250;
end

features.ratioLFMF = sum(dataLF)/sum(dataMF);






features.sumRegions = 2/(sum(diffPMC)+sum(diffNMC));

if isinf(features.sumRegions)
    features.sumRegions=250;
end





  Fs = fs;                      % samples per second
   dt = 1/Fs;                     % seconds per sample
   StopTime = 3;                  % seconds
   t = (0:dt:StopTime-dt)';
   N = size(t,1);
 
   x = signal;
   %% Fourier Transform:
   orf = fft(x);
   X = fftshift(orf);
   %% Frequency specifications:
   dF = Fs/N;                      % hertz
   f = -Fs/2:dF:Fs/2-dF;           % hertz

   
   amps = abs(X)/N;

   
   ws = floor((3*fs)/2);
   
   lim2Hz = ws+2*3; %375 je 0Hz, 1Hz=3smpl pøi 750 window a Fs 250hz
   limTopHz = ws+100*3;

   
   [ampZero id]=max(amps);
   
   f=f(lim2Hz:limTopHz);
   amps = amps(lim2Hz:limTopHz);
   
   
   [vals idx]=findpeaks(amps,'SortStr','descend');
 
   frequency=0;
   amplitude=0;
   
   if length(idx)>0
   
   frequency= f(idx(1));
   amplitude= vals(1);
   
   end;
   
   features.freqAb2Hz=frequency;
   features.ampFreqTop = amplitude/ampZero;


   f2 = frequency*2;
   f3 = frequency*3;
   
   amp2=0;
    amp3=0;
      
    if length(idx>0)
   idx2=floor(idx(1)+frequency*3);
   
   idx3=floor(idx2+frequency*3);


   
   if length(vals)>(idx2+2)
   amp2 = max(vals(idx2-1:idx2+1));
   end
   
 
   
  
   if length(vals)>(idx3+3)
   amp3 = max(vals(idx3-2:idx3+2));
   end
    end
   
   features.ampfreq2harm=amp2/ampZero;
   features.ampfreq3harm=amp3/ampZero;
%    



ads = abs(diff(signal));

md = prctile(ads,90);
ms = prctile(signal,90)-(prctile(signal,10));

features.ratioDiff = md/ms;

%nejmenší vr v 0.15sec oknì

vrs=[];

step =floor(0.075*fs);
win = floor(0.15*fs);


for s=1: step:length(signal)-win-1
    window = s:s+win;
    vrs(length(vrs)+1) = max(signal(window))-min(signal(window));
end

features.ratioMin015 = min(vrs)/ms;
features.ratioMed015 = median(vrs)/ms;



%%% pøidáno ze zlomù
n=length(signal);
nt=floor(n/3);
winCThird=nt:nt+nt;

sgn=signal';

acPost=fcr(sgn,sgn(winCThird));

prc=85;

medAcPost=prctile(acPost,prc);

features.medAcPost=medAcPost;

minDist = 0.1*fs; %toto ale omezí detekci šumu 0.15*fs

[pksPost locsPost]=findpeaks(acPost,'MinPeakDistance',minDist,'MinPeakHeight',medAcPost);

[freqPost ampPost]=getFrequency(acPost,fs,1,45);

features.acFreq=freqPost;
features.acRfreqRatio = features.acFreq/features.freqAb2Hz;


diffLocsPost=diff(locsPost)./fs;

features.meanRRPost=0.01;
if ~isempty(diffLocsPost)
    features.meanRRPost=mean(diffLocsPost);
end

features.meanPksPost=0.01;
if length(pksPost)>0
    features.meanPksPost=mean(pksPost);
end




end

