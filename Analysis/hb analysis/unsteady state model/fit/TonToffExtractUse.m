%%TonToffExtractUse
% 'Inherit1-AlpOFF-Turn1-TimeBin1-Combin0-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.05 0.05 0.1]-KU[0.05 0.05 0.6]-Alp0.6\'
%Inherit1-Alp0-Turn1-TimeBin1-Combin1-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.01 0.6]-KU[0.10 0.10 0.7]-TimeStart60-KDisLost-fix
FoldName=['Inherit1-AlpOn-Turn1-TimeBin0.5-Combin3-PonLock0-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle' ...
    '-Vel40-KL[0.01 0.12 1]-KU[0.2 0.12 1.3]-TimeStart60-KDisLost[10 10 1 1 1]-Artifical-Filter5-timeBinMove1-PonBoole0-re6-0520fix\'];
FoldNameCycle12='Inherit1-AlpOn-Turn1-TimeBin0.5-Combin3-PonLock0-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-Vel40-KL[0.01 0.12 1]-KU[0.1 0.12 1.4]-KDisLost[10 5 0 1.2 0]FirstON-Artifical-Filter6-timeBinMove1-PonBoole0-Star0731-3\';
    % ['Inherit1-AlpOn-Turn1-TimeBin0.5-Combin3-PonLock0-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle' ...
    % '-Vel40-KL[0.01 0.11 1]-KU[0.2 0.13 1.3]-KDisLost[5 2 0 0 0]FirstOFF-Artifical-Filter6-timeBinMove1-PonBoole0-LowArtK-Star0702\'];
FoldNameCycle13=...
    ['Inherit1-AlpOn-Turn1-TimeBin0.5-Combin3-PonLock0-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle' ...
    '-Vel40-KL[0.01 0.11 1]-KU[0.2 0.13 1.3]-KDisLost[5 5 1 1 1]FirstON-Artifical-Filter5-timeBinMove1-PonBoole0-re7-0520-2\'];
ELLim=[0.2 0.4];
%'Inherit1-AlpOn-Turn1-TimeBin0.5-Combin3-PonLock0-EL[0 0.6]bin0.1Move0.05-TimeWindowAdd-Vel45-KL[0.01 0.1 0.60]-KU[0.2 0.1 0.8]-TimeStart60-KDisLost[10 10 1 1 1]-Artifical-Filter5-timeBinMove1-PonBoole0-re4\';
SaveName='-TimeWindowMingle-PonBoole0-0726-bcd-ori5\';
%'Inherit1-AlpOn-Turn1-TimeBin0.5-Combin3-EL[0 0.6]bin0.1Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[40 10 1 1 5]-CombinPartArtifical-artificalFilter2-new\';
% 'Inherit1-AlpOn-Turn1-TimeBin0.5-Combin2-EL[0 0.6]bin0.2Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.60]-KU[0.10 0.08 0.75]-TimeStart60-KDisLost-CombinPartArtifical\';
%'Inherit1-AlpOn-Turn1-TimeBin0.5-Combin2-EL[0 0.6]bin0.2Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.6]-KU[0.10 0.08 0.7]-TimeStart60-KDisLost-CombinPartArtifical-fix\';
%'Inherit1-AlpOn-Turn1-TimeBin1-Combin2-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.6]-KU[0.10 0.08 0.7]-TimeStart60-KDisLost\';
%% Fix Time
Vel=40;
kfitTime=(3+283+2277/2)/Vel;Tstay=sum([3+283 2277 407+509])/Vel;
FixTime=0;
%%
HkCycleAll=figure('Position',[800 600 1000 400]);
pbaspect([3 1 1])
HkCycleAll2=figure('Position',[800 600 1000 400]);
pbaspect([3 1 1])
HkaCycleAll=figure('Position',[800 600 1000 400]);
pbaspect([3 1 1])
Hk1=figure('Position',[800 600 1000 400]);
Hk2=figure('Position',[800 600 1000 400]);
DivisionTime=[1,1,1]*60;
CycleTime=[7,9,13]*60;
CycleTimePlot=[0,7+1,7+1+9+1]*60;
Fuc=struct();
CycletimeExtractG=[];
pbaspect([3 1 1])

for CycleI=11:13
    CycleLabel=num2str(CycleI);

    
LoadRoadFun=@(x) ['Y:\TimeFit\OutputModel\Ton\Cycle',x,'\',...
    FoldName];
TimeAdd=[0 0 DivisionTime(CycleI-10) DivisionTime(CycleI-10) DivisionTime(CycleI-10)];
switch CycleLabel
    case '11'
        startPoint=[0.008 -0.008 240 120 10]+TimeAdd;%CYCLE11
        Ub=[inf 0 320 150 80]+TimeAdd;
        Lb=[0 -inf 210 80 -50]+TimeAdd;
        timeSeries=0:1:500;
    case '12'
        startPoint=[0.008 -0.008 260 100 0]+TimeAdd;%CYCLE12
        timeSeries=-60:1:600;
        Ub=[inf 0 410 200 50]+TimeAdd;
        Lb=[0 -inf 320 100 -20]+TimeAdd;
        if ~isempty(FoldNameCycle12)
            LoadRoadFun=@(x) ['Y:\TimeFit\OutputModel\Ton\Cycle',x,'\',...
                FoldNameCycle12];
        end
    case '13'
        startPoint=[0.0035 -0.0035 500 80 0]+TimeAdd;%CYCLE13
        Ub=[inf 0 680 120 50]+TimeAdd;
        Lb=[0 -inf 600 50 -50]+TimeAdd;
        timeSeries=-90:1:900;
        if ~isempty(FoldNameCycle13)
            LoadRoadFun=@(x) ['Y:\TimeFit\OutputModel\Ton\Cycle',x,'\',...
                FoldNameCycle13];
        end
end


OutputRoadFun=@(x) ['Y:\TimeFit\OutputModel\Ton\Cycle',x,'\',...
    SaveName];
OutputRoad=OutputRoadFun(CycleLabel);mkdir(OutputRoad);
OutputRoadCycleALL=OutputRoadFun('Plot');
load([LoadRoadFun(CycleLabel),'Results0.mat'])

Hka=figure('Position',[800 600 1000 400]);    
Hk=figure('Position',[800 600 1000 400]);
pbaspect([3 1 1])
% funFit=@(ton,ts1,ts2,toff,rs,ka,kb,x) piecewise(x<=ton,0,x>ton&x<ts1,ka*(x-ton),...
%     x>=ts1&x<=ts2,rs,x>ts2&x<toff,rs+kb*(x-ts2),x>=toff,0);
ft = fittype('piecewiseLine2(x,ton,ts1,toff,ka,kb)');
% 'independent', {'x'},'problem', {'ton','ts1','ts2','toff','ka','kb'}, ...
%  'dependent','y');
% c_map = [0.57, 0.69, 0.30
%      0.89, 0.88, 0.57
%      0.76, 0.49, 0.58
%      0.47, 0.76, 0.81
%      0.21, 0.21, 0.35
%      0.28, 0.57, 0.54
%      0.07, 0.35, 0.40
%      0.41, 0.20, 0.42
%      0.60, 0.24, 0.18
%      0.76, 0.84, 0.65];
c_map = [cool(ceil(size(EL,2)/2));spring(floor(size(EL,2)/2))];

timeExtractG=[];ELlegends={};ELlegendsOri={};

ELindex=mean(EL,1)>=ELLim(1)&mean(EL,1)<=ELLim(2);
ELuse=1:size(EL,2);
ELuse=ELuse(ELindex);
ELiReal=0;
ActiveRA=[];
for ELi=ELuse
    ELiReal=ELiReal+1;
    el=EL(:,ELi);
    ELlabels{ELi}=['EL ',num2str(el(1)),'-',num2str(el(2))];
    ELlegends=[ELlegends;['EL ',num2str(el(1)),'-',num2str(el(2))];['EL ',num2str(el(1)),'-',num2str(el(2))]];
    ELlegendsOri=[ELlegendsOri;['EL ',num2str(el(1)),'-',num2str(el(2))]];
    % PiniG=PiniELs{ELi,1};
    % MeanExp=MeanExpELs{ELi,:};
    
    if FixTime
        KannealG=KannealELs{ELi,1};
        TimeLabel=TimeMeanELs{ELi,:};
        KannealTime=TimeLabel'.*60-[Tstay kfitTime];
        Tstart=TimeLabel(1)*60-60;
        KannealTime(KannealTime(:,1)<Tstart,1)=Tstart;
        [timeBoundaries, featureValues] = createTimelines(KannealTime, KannealG);
        TimeLabelMean=mean([timeBoundaries(1:end-1);timeBoundaries(2:end)],1);
        KannealG=featureValues;
    else
        KannealG=KPlotannealELs{ELi,1};
        TimeLabelMean=TimePlotMeanELs{ELi,:};
    end
    TimeLabelMean=TimeLabelMean+DivisionTime(CycleI-10);
    % 
    ActiveR=KannealG(:,4);
    ActiveRA(:,ELiReal)=KannealG(:,4);
    f1 = fit(TimeLabelMean',ActiveR,ft,...
        'startpoint',startPoint,'upper',Ub,'lower',Lb);
    % ,...
        % 'Lower',[0.0025 -0.0055 400 60 -180 600],...
        % 'Upper',[0.0045 -0.0025 600 90 -130 900]);
        TS2real=-(f1.ka*(f1.ton-f1.ts1)-f1.kb*f1.toff)/f1.kb;
    timeExtractG(ELiReal,:)=[f1.ton f1.toff f1.ts1 TS2real f1.ton-f1.ts1 TS2real-f1.toff f1.ka f1.kb f1.ka*(f1.ton-f1.ts1) f1.toff-f1.ton TS2real-f1.ts1 CycleTime(CycleI-10)-f1.toff CycleTime(CycleI-10)-TS2real];
    

    figure(HkCycleAll)
    hold on
    scatter(TimeLabelMean+CycleTimePlot(CycleI-10),ActiveR,50,'filled',...
              'MarkerFaceColor',c_map(ELi,:))
    plot(timeSeries+CycleTimePlot(CycleI-10),f1(timeSeries),'LineWidth',2,'color',c_map(ELi,:))
    ylabel('ActiveUncleusRatio')
    xlabel('Time/sec');
    title(['CycleAll','Fitting ActiveUncleusRatio',' vs. time with different ELs'])

    figure(HkCycleAll2)
    hold on
    scatter(TimeLabelMean+CycleTimePlot(CycleI-10),ActiveR,50,'filled',...
              'MarkerFaceColor',c_map(ELi,:))
    ylabel('ActiveUncleusRatio')
    xlabel('Time/sec');
    title(['CycleAll','Fitting Mean ActiveUncleusRatio',' vs. time with different ELs'])
    
    figure(Hk)
    hold on
    scatter(TimeLabelMean,ActiveR,50,'filled',...
              'MarkerFaceColor',c_map(ELi,:))
    plot(timeSeries,f1(timeSeries),'LineWidth',2,'color',c_map(ELi,:))
    ylabel('ActiveUncleusRatio')
    xlabel('Time/sec');
    title(['Cycle',CycleLabel,'Fitting ActiveUncleusRatio',' vs. time with different ELs']) 
end
legend(ELlegends)
saveas(gcf,[OutputRoad,'Fitting ActiveUncleusRatio vs Time','.fig']);
saveas(gcf,[OutputRoad,'Fitting ActiveUncleusRatio vs Time','.png']);
mkdir(OutputRoadCycleALL)
saveas(gcf,[OutputRoadCycleALL,'Cycle',CycleLabel,' Fitting ActiveUncleusRatio vs Time','.fig']);
saveas(gcf,[OutputRoadCycleALL,'Cycle',CycleLabel,' Fitting ActiveUncleusRatio vs Time','.png']);

% mean pon
fa1 = fit(TimeLabelMean',mean(ActiveRA,2),ft,...
        'startpoint',startPoint,'upper',Ub,'lower',Lb);
fa2 = fit(repmat(TimeLabelMean',ELiReal,1),ActiveRA(:),ft,...
        'startpoint',startPoint,'upper',Ub,'lower',Lb);
figure(Hka)
hold on
scatter(TimeLabelMean,mean(ActiveRA,2),50,'filled',...
          'MarkerFaceColor',c_map(ELi,:))
plot(timeSeries,[fa1(timeSeries)';fa2(timeSeries)'],'LineWidth',2)
ylabel('ActiveUncleusRatio Mean')
xlabel('Time/sec');
title(['Cycle',CycleLabel,'Fitting ActiveUncleusRatio Mean',' vs. time with different ELs'])
eval(['Fuc.cycle',CycleLabel,'=fa1;'])
saveas(gcf,[OutputRoadCycleALL,'Cycle',CycleLabel,' Fitting Mean ActiveUncleusRatio vs Time','.fig']);
saveas(gcf,[OutputRoadCycleALL,'Cycle',CycleLabel,' Fitting Mean ActiveUncleusRatio vs Time','.png']);

figure(HkaCycleAll)
hold on
scatter(TimeLabelMean+CycleTimePlot(CycleI-10),mean(ActiveRA,2),50,'filled',...
          'MarkerFaceColor',c_map(ELi,:))
plot(timeSeries+CycleTimePlot(CycleI-10),[fa1(timeSeries)';fa2(timeSeries)'],'LineWidth',2)
ylabel('ActiveUncleusRatio Mean')
xlabel('Time/sec');
title(['CycleAll','Fitting ActiveUncleusRatio Mean',' vs. time with different ELs'])

figure(HkCycleAll2)
hold on
scatter(TimeLabelMean+CycleTimePlot(CycleI-10),mean(ActiveRA,2),80,'filled',...
          'MarkerFaceColor','g')
plot(timeSeries+CycleTimePlot(CycleI-10),[fa1(timeSeries)';fa2(timeSeries)'],'LineWidth',3)
% ELlegendsOri = [ELlegendsOri;'EL mean','Fitting'];

TS2real=-(fa1.ka*(fa1.ton-fa1.ts1)-fa1.kb*fa1.toff)/fa1.kb;
CycletimeExtractG(CycleI-10,:)=...
    [fa1.ton fa1.toff fa1.ts1 TS2real fa1.ton-fa1.ts1 TS2real-fa1.toff fa1.ka fa1.kb fa1.ka*(fa1.ton-fa1.ts1) fa1.toff-fa1.ton TS2real-fa1.ts1 CycleTime(CycleI-10)-fa1.toff CycleTime(CycleI-10)-TS2real];



YlabelK=["Ton","Toff","TonS","ToffE","ActivingTime","InactivingTime","Ka","kb","ActiveMax","TimeWindow","TimeWindow2",'End-Toff','End-ToffE'];
ELMean=mean(EL,1);ELMeanUse=ELMean(ELindex);
for ii=1:size(timeExtractG,2)
    figure(Hk1)
    subplot(4,4,ii)
    hold on
    plot(ELMeanUse,timeExtractG(:,ii),'-','marker','o','LineWidth',2)
    Kname=char(YlabelK(ii));
    ylabel(Kname)
    xlabel('EL');
    title([Kname,' vs. EL'])
end


for ii=1:size(timeExtractG,2)
    figure(Hk2)
    subplot(4,4,ii)
    hold on
    [B,TFrm] = rmoutliers(timeExtractG(:,ii));
    plot(ELMeanUse(~TFrm),timeExtractG(~TFrm,ii),'-','marker','o','LineWidth',2)
    Kname=char(YlabelK(ii));
    ylabel(Kname)
    xlabel('EL');
    title([Kname,' vs. EL'])
end
end
figure(HkCycleAll)
legend(ELlegends)
xlim([0 2000])
saveas(gcf,[OutputRoadCycleALL,'CycleAll Fitting ActiveUncleusRatio vs Time','.fig']);
saveas(gcf,[OutputRoadCycleALL,'CycleAll Fitting ActiveUncleusRatio vs Time','.png']);

figure(HkCycleAll2)
legend(ELlegendsOri)
xlim([0 2000])
saveas(gcf,[OutputRoadCycleALL,'CycleAll Fitting Mean2 ActiveUncleusRatio vs Time','.fig']);
saveas(gcf,[OutputRoadCycleALL,'CycleAll Fitting Mean2 ActiveUncleusRatio vs Time','.png']);

figure(HkaCycleAll)
xlim([0 2000])
saveas(gcf,[OutputRoadCycleALL,'CycleAll Fitting Mean ActiveUncleusRatio vs Time','.fig']);
saveas(gcf,[OutputRoadCycleALL,'CycleAll Fitting Mean ActiveUncleusRatio vs Time','.png']);

figure('Position',[200 300 2000 800]);
subplot(3,1,1)
Sindex=[1 2 3 4 10 11];
X = categorical(YlabelK(Sindex));
X = reordercats(X,YlabelK(Sindex));
bar(X,CycletimeExtractG(:,Sindex))
ylabel('Time/sec')
set(gca,'FontSize',20);
subplot(3,1,2)
Sindex=[5 6 12 13];
X = categorical(YlabelK(Sindex));
X = reordercats(X,YlabelK(Sindex));
bar(X,CycletimeExtractG(:,Sindex))
ylabel('Time/sec')
set(gca,'FontSize',20);
subplot(3,1,3)
Sindex=[9];
X = categorical(YlabelK(Sindex));
X = reordercats(X,YlabelK(Sindex));
bar(X,CycletimeExtractG(:,Sindex))
ylabel('Ratio')
ylim([0 1])
set(gca,'FontSize',20);
legend('Cycle11','Cycle12','Cycle13')
saveas(gcf,[OutputRoadCycleALL,'FittingResults Bar','.fig']);
saveas(gcf,[OutputRoadCycleALL,'FittingResults Bar','.png']);

figure('Position',[200 300 2000 800]);
subplot(3,1,1)
Sindex=[1 2 10 11];
X = categorical(YlabelK(Sindex));
X = reordercats(X,YlabelK(Sindex));
bar(X,CycletimeExtractG(:,Sindex))
ylabel('Time/sec')
set(gca,'FontSize',20);
subplot(3,1,2)
Sindex=[5 6];
X = categorical(YlabelK(Sindex));
X = reordercats(X,YlabelK(Sindex));
bar(X,CycletimeExtractG(:,Sindex))
ylabel('Time/sec')
set(gca,'FontSize',20);
subplot(3,1,3)
Sindex=[9];
X = categorical(YlabelK(Sindex));
X = reordercats(X,YlabelK(Sindex));
bar(X,CycletimeExtractG(:,Sindex))
ylabel('Ratio')
ylim([0 1])
set(gca,'FontSize',20);
legend('Cycle11','Cycle12','Cycle13')
saveas(gcf,[OutputRoadCycleALL,'FittingResults Bar Choose','.fig']);
saveas(gcf,[OutputRoadCycleALL,'FittingResults Bar Choose','.png']);


figure(Hk1)
legend('Cycle 11','Cycle 12','Cycle 13')
saveas(gcf,[OutputRoadCycleALL,'Fitting K',' vs EL','.fig']);
saveas(gcf,[OutputRoadCycleALL,'Fitting K',' vs EL','.png']);

figure(Hk2)
legend('Cycle 11','Cycle 12','Cycle 13')
saveas(gcf,[OutputRoadCycleALL,'Fitting K',' vs EL',' DelOutliers.fig']);
saveas(gcf,[OutputRoadCycleALL,'Fitting K',' vs EL',' DelOutliers.png']);

save([OutputRoadCycleALL,'FittingResults.mat'],'Fuc','YlabelK','CycletimeExtractG')