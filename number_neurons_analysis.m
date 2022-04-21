%% importing data from txt files
clearvars
file = importdata('chick28_2.20_SPK01b.txt'); % insert the name of the txt file
SPK = file.data(:,1); % select spike timestamps

%select timestamps of single trials
num1radSel = file.data(:,file.colheaders=="num1radSel");
num1radSel = num1radSel(~isnan(num1radSel));

num2radSel = file.data(:,file.colheaders=="num2radSel");
num2radSel = num2radSel(~isnan(num2radSel));

num3radSel = file.data(:,file.colheaders=="num3radSel");
num3radSel = num3radSel(~isnan(num3radSel));

num4radSel = file.data(:,file.colheaders=="num4radSel");
num4radSel = num4radSel(~isnan(num4radSel));

num5radSel = file.data(:,file.colheaders=="num5radSel");
num5radSel = num5radSel(~isnan(num5radSel));

num1areaSel = file.data(:,file.colheaders=="num1areaSel");
num1areaSel = num1areaSel(~isnan(num1areaSel));

num2areaSel = file.data(:,file.colheaders=="num2areaSel");
num2areaSel = num2areaSel(~isnan(num2areaSel));

num3areaSel = file.data(:,file.colheaders=="num3areaSel");
num3areaSel = num3areaSel(~isnan(num3areaSel));

num4areaSel = file.data(:,file.colheaders=="num4areaSel");
num4areaSel = num4areaSel(~isnan(num4areaSel));

num5areaSel = file.data(:,file.colheaders=="num5areaSel");
num5areaSel = num5areaSel(~isnan(num5areaSel));

num1perimSel = file.data(:,file.colheaders=="num1perimSel");
num1perimSel = num1perimSel(~isnan(num1perimSel));

num2perimSel = file.data(:,file.colheaders=="num2perimSel");
num2perimSel = num2perimSel(~isnan(num2perimSel));

num3perimSel = file.data(:,file.colheaders=="num3perimSel");
num3perimSel = num3perimSel(~isnan(num3perimSel));

num4perimSel = file.data(:,file.colheaders=="num4perimSel");
num4perimSel = num4perimSel(~isnan(num4perimSel));

num5perimSel = file.data(:,file.colheaders=="num5perimSel");
num5perimSel = num5perimSel(~isnan(num5perimSel));

num1sel = file.data(:,file.colheaders=="num1sel");
num1sel = num1sel(~isnan(num1sel));

num2sel = file.data(:,file.colheaders=="num2sel");
num2sel = num2sel(~isnan(num2sel));

num3sel = file.data(:,file.colheaders=="num3sel");
num3sel = num3sel(~isnan(num3sel));

num4sel = file.data(:,file.colheaders=="num4sel");
num4sel = num4sel(~isnan(num4sel));

num5sel = file.data(:,file.colheaders=="num5sel");
num5sel = num5sel(~isnan(num5sel));

numAllSel = file.data(:,file.colheaders=="numAllSel");
numAllSel = numAllSel(~isnan(numAllSel));

%% sorting variables
radAll = sortrows([num1radSel; num2radSel; num3radSel; num4radSel; num5radSel]);
radAll = radAll(~isnan(radAll));

areaAll = sortrows([num1areaSel; num2areaSel; num3areaSel; num4areaSel; num5areaSel]);
areaAll = areaAll(~isnan(areaAll));

perimAll = sortrows([num1perimSel; num2perimSel; num3perimSel; num4perimSel; num5perimSel]);
perimAll = perimAll(~isnan(perimAll));

numAllSel = sortrows(numAllSel(~isnan(numAllSel),1));
idxCellsRadius = ismember(numAllSel,radAll);
idxCellsArea = ismember(numAllSel,areaAll);
idxCellsPerim = ismember(numAllSel,perimAll);

idxNUM1log = ismember(numAllSel(~isnan(numAllSel)), num1sel(~isnan(num1sel))); 
idxNUM2log = ismember(numAllSel(~isnan(numAllSel)), num2sel(~isnan(num2sel))); 
idxNUM3log = ismember(numAllSel(~isnan(numAllSel)), num3sel(~isnan(num3sel))); 
idxNUM4log = ismember(numAllSel(~isnan(numAllSel)), num4sel(~isnan(num4sel))); 
idxNUM5log = ismember(numAllSel(~isnan(numAllSel)), num5sel(~isnan(num5sel))); 

cv = who('SPK*')';
spikes2 = cellfun(@eval,cv,'UniformOutput',0);
visualCellsNames = cv;



evtAllStack = table(numAllSel);

numer = repmat(1,size(evtAllStack,1),1);
numer(idxNUM2log) = 2;
numer(idxNUM3log) = 3;
numer(idxNUM4log) = 4;
numer(idxNUM5log) = 5;

evtAllStack.evtType = numer;
evtAllStack.Properties.VariableNames{1} = 'evtTime';

evtAllStack.fixed = repmat(1,size(evtAllStack,1),1);
evtAllStack.fixed(idxCellsArea) = 2;
evtAllStack.fixed(idxCellsPerim) = 3;

%}
%% cut the spike train into trials
cellSpike = {};

for iTrialEvt1 = 1 : size(evtAllStack,1)
   Evt = evtAllStack.evtTime(iTrialEvt1);
   for iCell = 1 : size(spikes2,2)
       cell = (spikes2{:,iCell});
   
      cellSpike{iCell, iTrialEvt1} = cell(cell < (Evt+1.9) & cell > (Evt-0.5),:)-Evt; % Evt +2 = 3 seconds trials

   end
end

%% calculate mean firing rate, select firing rates >= 1Hz

HighFR = cellfun(@(x)  length(x)/2.4, cellSpike, 'UniformOutput', 0);
idxHighFR = (cell2mat(HighFR)>=1)';

cellSpikeStimFR= cellfun(@(x)  length(x(x>0 & x<=0.5))/0.5, cellSpike, 'UniformOutput', 0);

cellSpikeStimFRmat = cell2mat(cellSpikeStimFR)';
visualCellsNamesMat = [visualCellsNames{:}];

fixed = cellstr(num2str(evtAllStack.fixed));
fixed(idxCellsRadius) = {'radius'};
fixed(idxCellsArea) = {'area'};
fixed(idxCellsPerim) = {'perim'};

numer = cellstr(num2str(evtAllStack.evtType));
numer(idxNUM1log) = {'num1'};
numer(idxNUM2log) = {'num2'};
numer(idxNUM3log) = {'num3'};
numer(idxNUM4log) = {'num4'};
numer(idxNUM5log) = {'num5'};

idxsession = numAllSel>=0;% | numAllSel<=2000;

%% perform ANOVA on firing rates during stimulus presentation
% check AnovaOUT variable for calculated p-values
% for number neurons we selected only cells with at least 7 trials per each
% condition*numerosity, only main effect of numerosity is significant
% (p<0.01)

clear AnovaOUT
for iCell = 1 : size(cellSpikeStimFRmat,2)
    try
    [p,t,stats] = anovan(cellSpikeStimFRmat(all([idxHighFR(:,iCell) idxsession],2),iCell),{fixed(all([idxHighFR(:,iCell) idxsession],2)),numer(all([idxHighFR(:,iCell) idxsession],2))},'model','interaction','varnames',{'fixed','numer'},'display','off');
    AnovaOUT(iCell,:) = p';
    end
end

AnovaOUT = array2table(AnovaOUT);
AnovaOUT.Properties.VariableNames = {'fixed','numer','fixedXnumer'};
AnovaOUT.Names = visualCellsNames';

clear Trial7Count
for iCell = 1 : size(idxHighFR,2)
    
    evtAllStack.HighFR = idxHighFR(:,iCell);
    groupsum = groupsummary(evtAllStack,{'evtType','fixed'},'sum','HighFR');
    Trial7Count(iCell) = all(groupsum.sum_HighFR>=7);
end

AnovaOUT.Trial7Count = Trial7Count';

idxNumerCells = AnovaOUT.numer < 0.01 & AnovaOUT.fixed > 0.01 & AnovaOUT.fixedXnumer > 0.01 & AnovaOUT.Trial7Count==1;

idxNumerCells = AnovaOUT.Trial7Count>0;
NumerCells = AnovaOUT(idxNumerCells,:);


cellSpikeNumer = cellSpike(idxNumerCells,:);
for k = 1:numel(cellSpikeNumer)
    cellSpikeNumer{k} = cellSpikeNumer{k}.';
end
visualCellsNamesNumer = visualCellsNames(:,idxNumerCells);

cellSpikeStimFRNumer = cellSpikeStimFRmat(:,idxNumerCells);
idxHighFRNumer = idxHighFR(:,idxNumerCells);
idxHighFRNumer = logical(idxHighFRNumer.*idxsession);

AnovaOUT
%% create a table with trial-by-trial mean firing rates 
Fstat = cell2mat(t(2:4,6)');
Pstat = cell2mat(t(2:4,7)');
statAll = [Fstat Pstat];

[c, m, h, nms] = multcompare(stats, 'Dimension',[2],'Display','off');
results = array2table(c);

for i = 1:size(results,1)
    results.Num1(i) = nms(results.c1(i),1);
    results.Num2(i) = nms(results.c2(i),1);    
end
results = results(:,6:8);

AllTrialsFR = table(cellSpikeStimFRmat(all([idxHighFR(:,iCell) idxsession],2),iCell), 'VariableNames',{'FR'});
AllTrialsFixed = cell2table(fixed(all([idxHighFR(:,iCell) idxsession],2)), 'VariableNames',{'Fixed'});
AllTrialsNumer = cell2table(numer(all([idxHighFR(:,iCell) idxsession],2)), 'VariableNames',{'Numer'});
AllTrials = [AllTrialsFR AllTrialsFixed AllTrialsNumer];


%% plot kernel psth
ymin = 0;
ymax = Inf;
cellsel = 1;
idxCellsFix =[]; %idxCellsRadius;% idxCellsArea; %idxCellsPerim

sigma = 0.1; %Standard deviation of the kernel = 15 ms

trTime = 2.3; 
edgesBins = -0.499 : 0.001 : 1.8;
edgesSigma=[-1*sigma:.001:1*sigma]; %Time ranges form -3*st. dev. to 3*st.dev. 
kernel = normpdf(edgesSigma,0,sigma); %Evaluate the Gaussian kernel
%kernel = kernel*.001; %Multiply by bin width so the probabilities sum to 1
center = ceil(length(edgesSigma)/2); %Find the index of the kernel center


idxNUM1plot = all([idxCellsFix,idxNUM1log],2); 
idxNUM2plot = all([idxCellsFix,idxNUM2log],2); 
idxNUM3plot = all([idxCellsFix,idxNUM3log],2); 
idxNUM4plot = all([idxCellsFix,idxNUM4log],2); 
idxNUM5plot = all([idxCellsFix,idxNUM5log],2); 

color1 = [255 127 0]./255;%[102 194 165]./255;%[0 0.447 0.741];%[255 127 0]./255;
color2 = [152 78 163]./255;%[252 141 98]./255;%[0.85 0.325 0.098];%[152 78 163]./255;
color3 = [77 175 74]./255;%[141 160 203]./255;%[0.929 0.694 0.125];%[77 175 74]./255;
color4 = [55 126 184]./255;%[231 138 195]./255;%[0.494 0.184 0.556];%[55 126 184]./255;
color5 = [228 26 28]./255;%[166 216 84]./255;%[0.466 0.674 0.188];%[228 26 28]./255;


for iFig = cellsel
    clearvars cellKernel1 cellKernel2 cellKernel3 KernelAll cellBinned1 cellBinned2 cellBinned3 cellBinned4 cellBinned5
    visualCellsNamesPlot = visualCellsNamesNumer(:,iFig);
    cellBinned1= cellfun(@(x)  histc(x,edgesBins), cellSpikeNumer(cellsel,all([idxNUM1plot idxHighFRNumer(:,cellsel)],2)), 'UniformOutput', 0);
    cellBinned2= cellfun(@(x)  histc(x,edgesBins), cellSpikeNumer(cellsel,all([idxNUM2plot idxHighFRNumer(:,cellsel)],2)), 'UniformOutput', 0);
    cellBinned3= cellfun(@(x)  histc(x,edgesBins), cellSpikeNumer(cellsel,all([idxNUM3plot idxHighFRNumer(:,cellsel)],2)), 'UniformOutput', 0);
    cellBinned4= cellfun(@(x)  histc(x,edgesBins), cellSpikeNumer(cellsel,all([idxNUM4plot idxHighFRNumer(:,cellsel)],2)), 'UniformOutput', 0);
    cellBinned5= cellfun(@(x)  histc(x,edgesBins), cellSpikeNumer(cellsel,all([idxNUM5plot idxHighFRNumer(:,cellsel)],2)), 'UniformOutput', 0);

    figure();

    zoom on
    for iCell = 1 : size(cellBinned1,1)

        binned1 = cellBinned1(iCell,:)';
        binned1 = cell2mat(binned1);
        binned2 = cellBinned2(iCell,:)';
        binned2 = cell2mat(binned2); 
        binned3 = cellBinned3(iCell,:)';
        binned3 = cell2mat(binned3);
        binned4 = cellBinned4(iCell,:)';
        binned4 = cell2mat(binned4);
        binned5 = cellBinned5(iCell,:)';
        binned5 = cell2mat(binned5);

        for iTrial = 1 : size(binned1,1)
            s = conv(binned1(iTrial,:),kernel); %Convolve spike data with the kernel
            s=s(center:trTime*1000+center-1); %Trim out the relevant portion of the spike density 
            cellKernel1(iTrial,:) = s; 
        end
        for iTrial = 1 : size(binned2,1) 
            s2 = conv(binned2(iTrial,:),kernel); %Convolve spike data with the kernel
            s2=s2(center:trTime*1000+center-1); %Trim out the relevant portion of the spike density 
            cellKernel2(iTrial,:) = s2; 
        end
        for iTrial = 1 : size(binned3,1) 
            s3 = conv(binned3(iTrial,:),kernel); %Convolve spike data with the kernel
            s3=s3(center:trTime*1000+center-1); %Trim out the relevant portion of the spike density 
            cellKernel3(iTrial,:) = s3;
        end
        for iTrial = 1 : size(binned4,1) 
            s4 = conv(binned4(iTrial,:),kernel); %Convolve spike data with the kernel
            s4=s4(center:trTime*1000+center-1); %Trim out the relevant portion of the spike density 
            cellKernel4(iTrial,:) = s4; 
        end
        for iTrial = 1 : size(binned5,1) 
            s5 = conv(binned5(iTrial,:),kernel); %Convolve spike data with the kernel
            s5=s5(center:trTime*1000+center-1); %Trim out the relevant portion of the spike density 
            cellKernel5(iTrial,:) = s5;
        end

    KernelAll = [cellKernel1;cellKernel2;cellKernel3; cellKernel4; cellKernel5];

    [lineOut1, fillOut1] = stdshade1(cellKernel1,0.1); hold on % 1==red
    lineOut1.Color = color1;
    fillOut1.FaceColor = color1;
    [lineOut2, fillOut2] = stdshade1(cellKernel2,0.1); hold on
    lineOut2.Color = color2;
    fillOut2.FaceColor = color2;
    [lineOut3, fillOut3] = stdshade1(cellKernel3,0.1); hold on
    lineOut3.Color = color3;
    fillOut3.FaceColor = color3;
    [lineOut4, fillOut4] = stdshade1(cellKernel4,0.1); hold on
    lineOut4.Color = color4;
    fillOut4.FaceColor = color4;
    [lineOut5, fillOut5] = stdshade1(cellKernel5,0.1); hold off
    lineOut5.Color = color5;
    fillOut5.FaceColor = color5;

    xlim([100 1900]); ylim([ymin ymax]);
    title('PSTH kernel size 100ms'); %!!!check 
    
    ax = gca; 
    ax.TitleFontSizeMultiplier = 0.8;
    
    max1 = max(mean(cellKernel1(:,:),1)); 
    max2 = max(mean(cellKernel2(:,:),1)); 
    max3 = max(mean(cellKernel3(:,:),1)); 
    max4 = max(mean(cellKernel4(:,:),1)); 
    max5 = max(mean(cellKernel5(:,:),1)); 
    maxAll = max([max1, max2, max3, max4, max5]);
    
    
    line([500 500],[0 maxAll+0.2*maxAll], 'LineStyle', '--','LineWidth',1); 
    line([1000 1000],[0 maxAll+0.2*maxAll], 'LineStyle', '--','LineWidth',1); 
     
    ax = gca;
    ax.XLimSpec = 'Tight';
    ax.TitleFontSizeMultiplier = 0.8;
    ax.XTick = [500 : 500 : 2000];
    ax.XTickLabel = [0 : 0.5 : 1.5];

    linkaxes(ax,'x');

end
end

%% plot tuning curve
for iCell = 1% 1:size(cellSpikeStimFRNumer,2)
    
[meanFR,stdFR,semFR, countTr, grpFR] = grpstats(cellSpikeStimFRNumer(idxHighFRNumer(:,iCell),iCell),{fixed(idxHighFRNumer(:,iCell)),numer(idxHighFRNumer(:,iCell))}, {'mean','std','sem','numel','gname'});
[grpFR, sortNUM] = sortrows(grpFR,2);
meanFR = meanFR(sortNUM,:);
stdFR = stdFR(sortNUM,:);
semFR = semFR(sortNUM,:);
countTr = countTr(sortNUM,:)';
 
[meanFRNumer,stdFRNumer,semFRNumer, grpFRNumer] = grpstats(cellSpikeStimFRNumer(idxHighFRNumer(:,iCell),iCell),{numer(idxHighFRNumer(:,iCell))}, {'mean','std','sem','gname'});
[grpFRNumer, sortNUMNumer] = sortrows(grpFRNumer,1);
meanFRNumer = meanFRNumer(sortNUMNumer,:)';
meanFRNumerNorm = (meanFRNumer-min(meanFRNumer))/(max(meanFRNumer)- min(meanFRNumer));
stdFRNumer = stdFRNumer(sortNUMNumer,:);
semFRNumer = semFRNumer(sortNUMNumer,:);

idxMeanFRRadius = ~cellfun('isempty', regexp(grpFR(:,1), 'radius')); 
idxMeanFRArea = ~cellfun('isempty', regexp(grpFR(:,1), 'area')); 
idxMeanFRPerim = ~cellfun('isempty', regexp(grpFR(:,1), 'perim')); 

ngrps = length(grpFR(idxMeanFRRadius));

try
figure()
tuncurv = errorbar((1:ngrps)', meanFR(idxMeanFRRadius), semFR(idxMeanFRRadius),'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle',':');
% legend({'radius'},{'area'},{'perimeter'},{'average'});
hold on
errorbar((1:ngrps)', meanFR(idxMeanFRArea), semFR(idxMeanFRArea),'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle','-.')
hold on
errorbar((1:ngrps)', meanFR(idxMeanFRPerim), semFR(idxMeanFRPerim),'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle','--')
hold on
errorbar((1:ngrps)', meanFRNumer, semFRNumer,'LineWidth',2,'Color','black','LineStyle','-')
hold off
set(gca,'xtick',1:ngrps,'xticklabel',grpFR(idxMeanFRRadius,2))
title('Tuning curve');
legend('radius','area','perimeter','average')

catch
    fprintf('loop number %d failed\n',iCell)
end
end

%% plot raster plot

evtStart = 0;
evtEnd = 0.5;
msize = 5;
figure();
for iCell =1%1 : size(iwant,2)

    visualCellsNamesPlot = visualCellsNamesNumer(:,iCell); % visualCellsNames2(range,:);

    idxCellsFix =[];% idxCellsArea;
    idxNUM1plot = all([idxCellsFix,idxNUM1log],2);
    idxNUM2plot = all([idxCellsFix,idxNUM2log],2); 
    idxNUM3plot = all([idxCellsFix,idxNUM3log],2); 
    idxNUM4plot = all([idxCellsFix,idxNUM4log],2); 
    idxNUM5plot = all([idxCellsFix,idxNUM5log],2); 
    cell1t = cellSpikeNumer(iCell,all([idxNUM1plot idxHighFRNumer(:,iCell)],2))';
    cell2t = cellSpikeNumer(iCell,all([idxNUM2plot idxHighFRNumer(:,iCell)],2))';
    cell3t = cellSpikeNumer(iCell,all([idxNUM3plot idxHighFRNumer(:,iCell)],2))';
    cell4t = cellSpikeNumer(iCell,all([idxNUM4plot idxHighFRNumer(:,iCell)],2))';
    cell5t = cellSpikeNumer(iCell,all([idxNUM5plot idxHighFRNumer(:,iCell)],2))';

    cellt = [cell1t; cell2t; cell3t; cell4t; cell5t];
    
    
    xPoints1 = [ cell1t{:} ];
    xPoints2 = [ cell2t{:} ];
    xPoints3 = [ cell3t{:} ];
    xPoints4 = [ cell4t{:} ];
    xPoints5 = [ cell5t{:} ];
    xPoints = [ cellt{:} ];

        % Getting the trials is trickier. 3 steps:
        % 1. First convert all the spike times into ones.
        trials = cellfun( @(x) {ones(size(x))}, cellt);
        % 2. Then multiply by trial number.
        for trialNum = 1:length(cellt)
            trials{trialNum} = trialNum*trials{trialNum};
        end
        % 3. Finally convert into a vector
        yPoints = [ trials{:} ];
        
        num1end = length(xPoints1);
         num2end = length(xPoints1)+length(xPoints2);
         num3end = num2end+length(xPoints3);
         num4end = num3end+length(xPoints4);
         num5end = num4end+length(xPoints5);
         
        yPoints1 = yPoints(1:num1end);
        yPoints2 = yPoints((num1end+1):num2end);
        yPoints3 = yPoints((num2end+1):num3end);
        yPoints4 = yPoints((num3end+1):num4end);
        yPoints5 = yPoints((num4end+1):num5end);
       
        plot(xPoints1,yPoints1,'.','MarkerSize',msize,'color',color1); hold on
        plot(xPoints2,yPoints2+20,'.','MarkerSize',msize,'color',color2); hold on
        plot(xPoints3,yPoints3+40,'.','MarkerSize',msize,'color',color3); hold on
        plot(xPoints4,yPoints4+60,'.','MarkerSize',msize,'color',color4); hold on
        plot(xPoints5,yPoints5+80,'.','MarkerSize',msize,'color',color5); hold off
    xlim([-0.4 1.4]); ylim([-20 max(yPoints5)+100]);%length(cellt)]);

    set(gca,'YTick', []);
    title('Raster plot');
    line([0 0],[0 max(yPoints5)+150], 'LineStyle', '--','LineWidth',2);
    line([0.5 0.5],[0 max(yPoints5)+150], 'LineStyle', '--','LineWidth',2);

end




%% create variables with firing rates averaged by condition*numerosity

iCell = 1;
    
[meanFRNumer,stdFRNumer,semFRNumer, grpFRNumer] = grpstats(cellSpikeStimFRNumer(idxHighFRNumer(:,iCell),iCell),{numer(idxHighFRNumer(:,iCell))}, {'mean','std','sem','gname'});
[grpFRNumer, sortNUMNumer] = sortrows(grpFRNumer,1);
meanFRNumer = meanFRNumer(sortNUMNumer,:)'
%meanFRNumerNorm = (meanFRNumer-min(meanFRNumer))/(max(meanFRNumer)- min(meanFRNumer));
stdFRNumer = stdFRNumer(sortNUMNumer,:);
semFRNumer = semFRNumer(sortNUMNumer,:);

[meanFRNumerFixed,stdFRNumerFixed,semFRNumerFixed, grpFRNumerFixed] = grpstats(cellSpikeStimFRNumer(idxHighFRNumer(:,iCell),iCell),{numer(idxHighFRNumer(:,iCell)),fixed(idxHighFRNumer(:,iCell))}, {'mean','std','sem','gname'});

meanFRNumerFixed = array2table(meanFRNumerFixed,'VariableNames',{'meanFR'});
meanFRNumerFixed.Numer = grpFRNumerFixed(:,1);
meanFRNumerFixed.Fixed = grpFRNumerFixed(:,2);

meanFRNumerFixed = sortrows(meanFRNumerFixed,[3 2]);
meanFRNumerFixedClip = [table2array(meanFRNumerFixed(1:5,1)) table2array(meanFRNumerFixed(6:10,1)) table2array(meanFRNumerFixed(11:15,1))]';

