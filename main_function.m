%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main code for computing LQI and decide liquefaction class based on the
% logistic regression model and optimal LQI threshold value. It requires
% Signal Processing Toolbox. 
%
% Author   : Weiwei Zhan
% Contact  : wzhan@g.clemson.edu
% Last edit: Feburary 08, 2021
% 
% Cite and credit:
% Zhan, W. and Chen, Q. (2021). "An accelerogram-based method for 
% quick assessment of liquefaction occurrence", Journal of Geotechnical 
% and Geoenvironmental Engineering.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up the input and output directory
clear all; close all; clc; 

plotType = 1;    % plot signal processing results for each accelerogram? 1 for YES; 0 for NO
mkdir('output'); % directory for signal processing figures

dataset = 'input';  % select folder name for input data
datadir = fullfile(dataset);
addpath(genpath(datadir));
tempSUM     = dir(fullfile(datadir,'**','*.csv'));
num_tempSUM = length(tempSUM);

%% Compute the RL and MIFr using signal processing techniques

for i=1:num_tempSUM
    ga=[];  acc=[];   time=[];       % clear variables that may have different length among ground motion records
    filename=tempSUM(i).name;
    
    ga   = csvread(filename);        % read acceleration time history data
    acc  = ga(:,2:3);                % input two horizontal components
    time = ga(:,1);                  % time vector                 
    dt   = time(2)-time(1);          % time step
    fs   = round(1/dt);              % sampling frequency
    L    = length(acc);              % record length 
    
    PGA(i,1:2) = max(abs(acc(:,1:2)));                % compute peak ground acceleration for each horizontal component
    PGA(i,3)   = max(sqrt(acc(:,1).^2+acc(:,2).^2));  % compute peak ground acceleration for the vector sum of the two horizontal components
    
    f1      = 1;                     % upper frequency limit of the low frequency component, see Equation(2)
    f2      = 10;                    % upper frequency limit of the whole-frequency component, see Equation(2)
    windowT = 10;                    % length of the time interval for MIFr computation

    for j=1:2 % loop over two horizontal components
         plotName  = strrep(filename,'.csv',sprintf('_Component%d_',j)); % get plot name of signal processing results for each horizontal component 
         MIFr(i,j) = FrequencyDropRatio(acc(:,j),fs,windowT,plotType,plotName);
         RL(i,j)   = LowFrequencyPortion(acc(:,j),fs,f1,f2,plotType,plotName); 
    end   
    close all    
   
    RecordName{i,1} = strrep(filename,'.csv','');  % get record name for each row
end

% compute the mean of the two horizontal components
MIFr(:,3) = mean(MIFr,2);
RL(:,3)   = mean(RL,2);

temp  = [PGA,MIFr,RL];
data = array2table(temp,'VariableNames',{'PGA_1','PGA_2','PGA',...  
    'MIFr_1','MIFr_2','MIFr','RL_1','RL_2','RL'},'RowNames',RecordName);

%% filter out low-intensity ground motions that could not trigger liquefaction
pgafilter = 1; 
if pgafilter == 1
   lowPGA_index = find(data.PGA < 0.08); % find records with pga less than the PGA threshold 
end

%% Compute LQI using the logistic regression model
data.LQI = 1./(1+exp(6.44-47.61.*data.RL.*data.MIFr)); % logistic regression model, namely Equation (12) 
data.LQI(lowPGA_index) = 0;     % reset LQI of low PGA cases as zero 
for k = 1:length(lowPGA_index)  % display index of low PGA cases
    disp(sprintf("The %dth accelerogram is classified as NonLiquefied because PGA less than 0.08g",lowPGA_index(k)));
end

%% classify liquefaction class as positive when computed LQI exceeds 0.15
data.LQclass(:) = 0;
LQI_th   = 0.15; % LQI threshold for binary classification
LQ_index = find(data.LQI >= LQI_th);
data.LQclass(LQ_index) = 1;

writetable(data,'output/LQI_results.xlsx','WriteRowNames',true);


%% visualize the liquefaction classification results  
figure
%scatter(data.RL,data.MIFr,'ko')
h = gscatter(data.RL,data.MIFr,data.LQclass,'kk','oo',6,'on','RL','MIFr'); hold on
if length(h)>1
set(h(2), 'MarkerFaceColor', 'k'); hold on;
end
hold on; box on
RL_th   = [0.01:0.01:1]; 
MIFr_th = (6.44-log(1/LQI_th-1))/47.61./RL_th;
plot(RL_th,MIFr_th,'k-','linewidth',1.0); % plot the decision boundary line
xlim([0,1])
ylim([-0.5,1.0])
xlabel('RL')
ylabel('MIFr')
legend({'NonLQ','LQ',sprintf('LQI = %.2f',LQI_th)},'location','northeast')
set(gca,'fontsize',12,'fontname','times');
set(gcf, 'Position', [600 200 600 480]);
print('output/LQI_results','-djpeg','-r300');



