function  RL = LowFrequencyPortion(acc,fs,freq_low,freq_high,plotType,plotName) 

% Author   : Weiwei Zhan
% Contact  : wzhan@g.clemson.edu
% Last edit: Feburary 08, 2021
% Citation: Zhan, W., Chen, Q. (2021). "An accelerogram-based method for 
% quick assessment of liquefaction occurrence", Journal of Geotechnical 
% and Geoenvironmental Engineering.

% Compute the ratio of low-frequency portion to the whole area of the 
% Fourier amplitude spectrum (RL) based on Fast Fourier Transform (FFT).

% [INPUT]
% acc           acceleration time history, unit is cm/s2
%   
% fs            sampling frequency
% 
% freq_low      upper frequency limit for low-frequency portion
% 
% freq_high     upper frequency limit for the whole frequency range 
% 
% plotType      1- plot and save; 0 - no plot;
% 
% plotName      name of te plot


% [OUTPUT]
% RL          Ratio of the low-frequency portion to the whole area of the Fourier amplitude spectrum 

%% Use Fast Fourier transform (FFT) to compute Fourier amplitude spectrum
acc = acc.*981;
L   = length(acc);
Y   = fft(acc);
P2  = abs(Y/L);
P1  = P2(1:floor(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

%% Compute RL following Equation(2)
temp1  = find(f<=freq_high);
fcut1  = f(temp1);
P1cut1 = P1(temp1);
temp2  = find(f<=freq_low);
fcut2  = f(temp2);
P1cut2 = P1(temp2);

RL = sum(P1cut2)/sum(P1cut1);


%% Visualize the computation of RL    
if plotType == 1
    figure
    plot(f,P1,'k-','linewidth',0.5) 
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (cm/s^{2})')
    mark = sprintf('RL = %.2f',RL);
    text(6,max(P1),mark,'fontsize',9,'fontname','times');
    xlim([0,freq_high]);
    ylim([0,1.1*max(P1)]);
    hold on
    idx = f>0 & f<freq_low;
    H   = area(f(idx),P1(idx));
    set(H,'FaceColor',[0 0.0 0]+0.5);
    hold off
    box on
    set(gca,'fontsize',9,'fontname','times');
    set(gcf, 'Position', [1000 500 300 240]);
    plotName = strcat('output/',plotName,'_RL.jpeg');
    print(plotName,'-djpeg','-r300');
end
