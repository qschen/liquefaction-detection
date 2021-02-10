function MIFr = FrequencyDropRatio(acc,fs,windowT,plotType,plotName)

% Author   : Weiwei Zhan
% Contact  : wzhan@g.clemson.edu
% Last edit: Feburary 08, 2021
% Citation: Zhan, W., Chen, Q. (2021). "An accelerogram-based method for 
% quick assessment of liquefaction occurrence", Journal of Geotechnical 
% and Geoenvironmental Engineering.

% Compute the mean instantaneous frequency decrease rate (MIFr) based on 
% the time-frequency analysis results from Short-time Fourier Transform 
% (STFT).

% [INPUT]
% acc           acceleration time history, unit is cm/s2
%   
% fs            sampling frequency
% 
% time          time vector
% 
% windowT       window length for computing frequency change ratio before
%               and after time of PGA
% 
% plotType      1- plot and save; 0 - no plot;
% 
% plotName      name of te plot


% [OUTPUT]
% MIFr          Mean Instaneous frequency decrease rate

%% STFT analysis
dt = 1.0/fs;
time=linspace(0,(length(acc)-1)*dt,length(acc)); %time vector
% define parameters for STFT analysis 
wlen = round(2.56/dt);              % window length (recomended to be power of 2); 2.56s window
hop  = 2;                           % hop size (recomended to be power of 2) 
nfft = wlen;                        % number of fft points (recomended to be power of 2);
win  = hamming(wlen,'periodic');    % window function
[s,f,t,ps] = spectrogram(acc,win,wlen-hop,nfft,fs,'yaxis');  % ps is power spectrum density; wlen-hop is the overlap length
 
%% compute mean Instaneous frequency (MIF) following Equation(4)
df          = f(2)-f(1);
numerator   = (f'*ps).*df;
temp        = sum(ps,1);
denominator = temp.*df;
MIF         = numerator./denominator;
MIF(isnan(MIF))=0;  % remove NaN at two ends where ps==0

        
%% compute mean Instaneous frequency decrease rate(MIFr) following Equation(5)
[pga tindex] = max(abs(acc));
tPGA    = time(tindex);
tindex1 = max(tPGA-windowT,t(1));
tindex2 = min(tPGA+windowT,t(end));
temp1   = find( t <= tPGA & t > tindex1 );
if isempty(temp1)  % PGA occured at the very beginning, tPGA < winT
    temp1 = 1;
end
temp2 = find( t >  tPGA & t < tindex2 );

MIFr  = (mean(MIF(temp1))-mean(MIF(temp2)))/mean(MIF(temp1));
    

%% visualize the MIF calculation results
if plotType==1
    figure
    plot([tPGA,tPGA],[0,1.5*max(MIF)],'r--','linewidth',0.8);hold on
    plot(t,MIF,'k-','linewidth',0.5)
    axis tight
    xlim([tindex1,tindex2])
    ylim([0,1.1*max(MIF)])
    xline(tPGA,'r--','linewidth',1.0)
    xlabel('Time (s)')
    ylabel('MIF (Hz)')  % Mean Instaneous frequency
    mark = sprintf('MIFr = %.2f',MIFr);
    text(0.5*max(t),max(MIF),mark,'fontsize',9,'fontname','times');
    hold off
    set(gca,'fontsize',9,'fontname','times');
    set(gcf, 'Position', [1000 500 300 240]);
    plotName = strcat('output/',plotName,'_MIFr.jpeg');
    print(plotName,'-djpeg','-r300');    
end
