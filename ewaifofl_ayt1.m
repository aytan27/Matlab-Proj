%aytan-unwrap works,
%can use for EWAIF OFL value prediction after simulations determine
%the lowest rms attn(slope)..
%enter center, deltaf etc values, then cut/paste ewaifofl line


clear all
close all
randn('seed',sum(100*clock));
rand('seed',sum(100*clock));


center=1000;
deltaf=40; %40
attn=70; %70 db SPL filter skirt


SR=25000; %25K
sr=1/SR;
T=.5; %200ms
N=T*SR;
t=[0:sr:T-sr];
spec=zeros(1,N);
fo=1/T;
amp = 1000;
pertsigma = 1;

avgthresh = [];
freqharm = [481 501 521];
% freqharm = [496 501 506];

nblcks = 10;

nn = 500;
ncomp = 3;




% sigdB = -5;
% % sigamp = 1000*10^(sigdB/20);
% sigamp_attn = 10 ^ ( sigdB / 20 )*10^(-attn * oct /20); %sig level amp in db

for bigloop=1:nblcks;
    
    
    meanstanfreq= zeros(1,nn);
    meansigfreq = zeros(1,nn);
    
    thresh=zeros(1,nn+3);
    rev=zeros(1,nn);
    down=0;
    sigdB = -10;
    
    amp1 =amp; %? reference component amplitude on low frequency tone
    oct=log2(center/(center-deltaf)); %octave, 1000/960
    amp2= 10^(-attn * oct /20)*amp1; %octave distance from reference tone
    oct = log2((center + deltaf)/(center-deltaf)); %skirt?1040/960
    amp3 = 10^(-attn * oct /20)*amp2;
    oct=log2(center/(center-deltaf)); %1000/960
    
    % sigamp_attn = 10 ^ ( sigdB / 20 );
    sigamp_attn = 10 ^ ( sigdB / 20 )*10^(-attn * oct /20); %attenuate sig level amp in db
%     sigamp_attn = 10^(-attn * oct /20); %attenuate sig level amp in db
    sigamp = amp2*sigamp_attn;
    for kk = 1:nn;
        %%Setting up Attenuation skirt%%%
        
        
        standard = zeros(1,N);
        signal = zeros(1,N);
        
        %STANDARD
        phase=2*pi*rand(1,ncomp);
        pertstandard=pertsigma*randn(1,ncomp);    %adding perturbation to standard
        standard(freqharm(1))=amp1*10.^(pertstandard(1)/20).*(sin(phase(1))+1i*cos(phase(1))); %adding the phases from sin and cosine while amp is
        standard(freqharm(2))=amp2*10.^(pertstandard(2)/20).*(sin(phase(2))+1i*cos(phase(2))); %adding the phases from sin and cosine while amp is
        standard(freqharm(3))=amp3*10.^(pertstandard(3)/20).*(sin(phase(3))+1i*cos(phase(3))); %adding the phases from sin and cosine while amp is
        
        %     standard(freqharm)=amp.*(sin(phase)+i*cos(phase));
        
        
        %SIGNAL
        
        phase=2*pi*rand(1,ncomp);
        pertsig=pertsigma*randn(1,ncomp);            %adding perturbation to signal (change sign)
        phasesig=phase((ncomp+1)/2);
        %     signal(freqharm)=amp.*(sin(phase)+i*cos(phase));
        
        signal(freqharm(1)) = amp1.*10.^(pertsig(1)/20).*(sin(phase(1))+1i*cos(phase(1)));
        signal(freqharm(2)) = amp2.*10.^(pertsig(2)/20).*(sin(phase(2))+1i*cos(phase(2)));
        signal(freqharm(3)) = amp3.*10.^(pertsig(3)/20).*(sin(phase(2))+1i*cos(phase(3)));
        signal(freqharm(2)) = signal(freqharm(2)) + sigamp*(sin(phasesig)+1i*cos(phasesig));
        
        
        standard = (ifft(standard)*N); %time dom
        signal   = (ifft(signal)*N);
        
        
        %%%%%%%%EWAIF for Standard complex
        
        x=imag(standard);
        y = real(standard);
        instphase = atan(y./x); %arctangent
        instphaseu = unwrap (instphase, pi/2);
        W = diff(instphaseu); %angular velocity which is the derivative of the instantaneaous phase
        
        ind=find(W>0); %diff between phases should NOT be positive
        W(ind)=W(ind)-pi; %so find pos values unwrap - pi
        
        index =(W < -pi); %to find deviant values that switched 180 degrees...
        W(index) = W(index)+pi;
        instfreq =abs(W / (sr * 2 * pi)); %2pi to convert to hz
        
        stan_env = sqrt(real(hilbert(y)).^2 + imag(hilbert(x)).^2);
        stan_env(1) = [];
        
        stan_ewaif = sum(instfreq.*stan_env)/sum(stan_env);
        meanstanfreq(kk)=stan_ewaif;
        %%%%%%%%%%%%%%%%%%%%%%%%%Signal EWAIF
        
        x = imag(signal);
        y = real(signal);
        instphase = atan(y./x); %arctangent
        instphaseu = unwrap (instphase, pi/2);
        W = diff(instphaseu); %angular velocity which is the derivative of the instantaneaous phase
        
        ind=find(W>0); %diff between phases should NOT be positive
        W(ind)=W(ind)-pi; %so find pos values unwrap - pi
        
        index =(W < -pi); %to find deviant values that switched 180 degrees...
        W(index) = W(index)+pi;
        instfreq =abs(W / (sr * 2 * pi)); %2pi to convert to hz
        
        sig_env = sqrt(real(hilbert(y)).^2 + imag(hilbert(x)).^2);
        sig_env(1) = [];
        
        sig_ewaif = sum(instfreq.*sig_env)/sum(sig_env);
        
        meansigfreq(kk) = sig_ewaif;
        
        if rand < 0.5;
            stim = 2;
            if sig_ewaif > stan_ewaif
                resp = 2;
            else
                resp = 1;
            end
            
        else
            stim = 1;
            if sig_ewaif > stan_ewaif
                resp = 1;
            else
                resp = 2;
            end;
        end;%if
        
        if sig_ewaif < stan_ewaif
            % record miss
            track(kk)=0;
            % raise signal 2 dB
            sigamp=sigamp*10^(0.1);
            % round signal
            thresh(kk+3)=20*log10(sigamp/amp);
            % reversal
            if kk>2
                if track(kk-2:kk)==[1 1 0]
                    rev(kk)=1;
                end
            end
            down=0;
        else
            track(kk)=1;
            
            if down==1
                %lower signal 2 dB
                sigamp=sigamp*10^(-0.1);
                %round signal
                thresh(kk+3)=20*log10(sigamp/1000);
                down=0;
            else
                thresh(kk+3)=thresh(kk+2);
                down=1;
            end
            %reversal
            if kk > 2
                if track(kk-2:kk)==[0 1 1]
                    rev(kk)=1;
                end
            end
        end
        %######################################################################3
        
        
        
        data(kk,:)=[stim resp pertstandard pertsig];
        
    end;
    
    ind1 = (bigloop-1)*500+1;
    ind2 = bigloop*500;
    DATA(ind1:ind2,:) = data; %indexing start at 1
    
    MEANsigfreq(ind1:ind2) = meansigfreq;
    MEANstanfreq(ind1:ind2) = meanstanfreq;
    avgthresh(bigloop)=mean(thresh);
end
meanavgthresh = mean(avgthresh);

mstan=mean(MEANstanfreq);
meanstanf = ['Mean StanFreq is ' ,num2str(mstan)];
disp(meanstanf)

msig=mean(MEANsigfreq);
meansigf = ['Mean SigFreq is ',num2str(msig)];
disp(meansigf);
StanDev_StandardFreq = std(MEANstanfreq)
StanDev_SigFreq = std(MEANsigfreq)

% keyboard
% save Inter_simpitch1 DATA meanavgthresh MEANsigfreq MEANstanfreq

%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA1= DATA;
ncomp = 3;
prtbstandard=DATA1(:,3:5);                   %all the rows for column 3 through column 5 for pertstd
prtbsig=DATA1(:,6:8);                        %all the rows for column 6 through column 8 pertsig
Int1= (prtbstandard - prtbsig);             %% HIT, Int 1 when signal is in interval 1,  difference from pertsig to pertstd
Int2= (Int1)* -1;
datainterval1 = [DATA1(:,1:2) Int1];            %creating matrix with all the rows for column 1 to Column 2 and 4
datainterval2 = [DATA1(:,1:2) Int2];            %3 prtb difference values for when signal was in the first interval
%and another datainterval2
%%matrix to account for
%%when signal is in 2nd interval

x=find(datainterval1(:,1)==1);               %check datainterval1 matrix for Column 1-STIM, where the value ==1 (sig standard)
data1=datainterval1(x,:);                     %data1 contains when stim ==1, for all columns associated
y=find(datainterval2(:,1)==2);               %looking in matrix datainterval2 for column 1 when value ==2
data2=datainterval2(y,:);
%deleting column 1 from data
data1(:,1)=[];                      %all rows for column one deleted
data2(:,1)=[];                      %all rows for column one from data2 deleted



%R P1 P2 P3
%R   1
%P1     1
%P2        1
%P3          1

r1=corrcoef(data1);                  %get the correlated coefficients for vector data 1
w1=r1(1,2:4);                       %matrix row 1 from 2column to the 4column holds the corcoeff needed
maxw1=max(abs(w1));                  %highest abs value from the corrcoeff matrix
w1=w1/maxw1;                          %normalize over the maximum value
% w1=-1*w1 ;                           %accounting for the different intervals

r2=corrcoef(data2);                 %correlation coefficient for set of data when stim was 2
w2=r2(1,2:4);                       %weights are values found in row 1 from column 2-4.
maxw2=max(abs(w2));                  %normalizing by the largest value weight
w2= w2/maxw2 ;
% w2= -1*w2   ;                        %sign change here for when not above

wgts=(w1+w2)/2;
rms=sqrt((w1-w2)*(w1-w2)'/ncomp);

h1=figure;
set(gcf, 'PaperPositionMode','auto')
ha=gca;
set(ha,'XTick',1:1:3);
axis([1,3,-1,1]);
hold on
plot(w1); %blue
plot([1 3],[0 0], 'k-');
plot(w2,'r');
plot(wgts,'gp');
title('Pitch');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ORIG
% 
% A1=1; %? reference component amplitude on low frequency tone
% oct=log2(center/(center-deltaf)); %octave, 1000/960
% A2= 10^(-attn * oct /20)*A1; %octave distance from reference tone
% oct = log2((center + deltaf)/(center-deltaf)); %skirt?1040/960
% A3 = 10^(-attn * oct /20)*A1;
% oct=log2(center/(center-deltaf)); %1000/960 ?
% 
% sig = 10 ^ ( ampdB / 20 )*10^(-attn * oct /20); %sig level amp in db
% 
% A2 = A2 + sig;
% 
% n=100;
% meanfreq=zeros(1,n);
% 
% for bigloop=1:n;
%      
% phase  = rand(1,3) * 2 * pi; %one by 3 vector of random phases
% 
% a = (center/fo) - (deltaf/fo) +1;
% b = (center/fo) + 1;
% c = (center/fo) + (deltaf/fo) + 1;
% 
% spec(a) = A1/2*(sin(phase(1)) - i*cos(phase(1)));
% spec(b) = A2/2*(sin(phase(2)) - i*cos(phase(2)));
% spec(c) = A3/2*(sin(phase(3)) - i*cos(phase(3)));
% 
% wavet=(ifft(spec)*N);  %time fxn
% 
% x=real(wavet);
% y=imag(wavet);
% instphase = atan(x./y); %arctangent
% instphaseu = unwrap (instphase, pi/2);
% W = diff(instphaseu); %angular velocity which is the derivative of the instantaneaous phase
% 
% ind=find(W>0); %diff between phases should NOT be positive
% W(ind)=W(ind)-pi; %so find pos values unwrap - pi
% 
% index =(W < -pi); %to find deviant values that switched 180 degrees...
% W(index) = W(index)+pi;
% instfreq =abs(W / (sr * 2 * pi)); %2pi to convert to hz
% 
% %x is the real part of the wave so here we use x rather than wave
% %%env = sqrt(wave.*wave + hilbert(wave).*hilbert(wave));
% % env = sqrt(x.*x + imag(hilbert(x)).^2);
% env = sqrt(real(hilbert(x)).^2 + imag(hilbert(x)).^2);
% env(1) = [];
% 
% ewaif = sum(instfreq.*env)/sum(env);
% meanfreq(bigloop)=ewaif;
% end
% 
% mean(meanfreq) %nochange from bgbver
% std(meanfreq)
% 
