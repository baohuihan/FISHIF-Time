%%TonToffExtract2
% 'Inherit1-AlpOFF-Turn1-TimeBin1-Combin0-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.05 0.05 0.1]-KU[0.05 0.05 0.6]-Alp0.6\'
%Inherit1-Alp0-Turn1-TimeBin1-Combin1-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.01 0.6]-KU[0.10 0.10 0.7]-TimeStart60-KDisLost-fix
FoldName='Inherit1-AlpOn-Turn1-TimeBin0.5-Combin3-PonLock0-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-Vel40-KL[0.01 0.1 0.60]-KU[0.2 0.1 0.8]-TimeStart60-KDisLost[10 10 1 1 1]-Artifical-Filter5-timeBinMove1-PonBoole0-re4-0507\';
%'Inherit1-AlpOn-Turn1-TimeBin0.5-Combin3-EL[0 0.6]bin0.1Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[40 10 1 1 5]-CombinPartArtifical-artificalFilter2-new\';
% 'Inherit1-AlpOn-Turn1-TimeBin0.5-Combin2-EL[0 0.6]bin0.2Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.60]-KU[0.10 0.08 0.75]-TimeStart60-KDisLost-CombinPartArtifical\';
%'Inherit1-AlpOn-Turn1-TimeBin0.5-Combin2-EL[0 0.6]bin0.2Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.6]-KU[0.10 0.08 0.7]-TimeStart60-KDisLost-CombinPartArtifical-fix\';
%'Inherit1-AlpOn-Turn1-TimeBin1-Combin2-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.6]-KU[0.10 0.08 0.7]-TimeStart60-KDisLost\';
%% Fix Time
Vel=40;
kfitTime=(3+283+2277/2)/Vel;Tstay=sum([3+283 2277 407+509])/Vel;
FixTime=0;
%%
Hk1=figure('Position',[800 600 1000 400]);
Hk2=figure('Position',[800 600 1000 400]);
DivisionTime=[1,1,1.5]*60;
CycleTime=[7,9,13]*60;

pbaspect([3 1 1])
for CycleI=11:13
    CycleLabel=num2str(CycleI);
OutputRoadFun=@(x) ['Y:\TimeFit\OutputModel\Ton\Cycle',x,'\',...
    FoldName];
OutputRoad=OutputRoadFun(CycleLabel);
OutputRoadCycleALL=OutputRoadFun('Plot');
load([OutputRoad,'Results0.mat'])

TimeAdd=[0 0 DivisionTime(CycleI-10) DivisionTime(CycleI-10) DivisionTime(CycleI-10)];
switch CycleLabel
    case '11'
        startPoint=[0.008 -0.008 240 120 10]+TimeAdd;%CYCLE11
        Ub=[inf 0 320 200 80]+TimeAdd;
        Lb=[0 -inf 210 100 -50]+TimeAdd;
        timeSeries=-100:1:500;
    case '12'
        startPoint=[0.008 -0.008 260 100 0]+TimeAdd;%CYCLE12
        timeSeries=-100:1:600;
        Ub=[inf 0 410 200 50]+TimeAdd;
        Lb=[0 -inf 220 50 -20]+TimeAdd;
    case '13'
        startPoint=[0.0035 -0.0035 500 80 0]+TimeAdd;%CYCLE13
        Ub=[inf 0 680 120 50]+TimeAdd;
        Lb=[0 -inf 600 50 -50]+TimeAdd;
        timeSeries=-100:1:900;
end
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

timeExtractG=[];ELlegends={};

for ELi=1:size(EL,2)
    el=EL(:,ELi);
    ELlabels{ELi}=['EL ',num2str(el(1)),'-',num2str(el(2))];
    ELlegends=[ELlegends;['EL ',num2str(el(1)),'-',num2str(el(2))];['EL ',num2str(el(1)),'-',num2str(el(2))]];

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
    f1 = fit(TimeLabelMean',ActiveR,ft,...
        'startpoint',startPoint,'upper',Ub,'lower',Lb);
    % ,...
        % 'Lower',[0.0025 -0.0055 400 60 -180 600],...
        % 'Upper',[0.0045 -0.0025 600 90 -130 900]);
        TS2real=-(f1.ka*(f1.ton-f1.ts1)-f1.kb*f1.toff)/f1.kb;
    timeExtractG(ELi,:)=[f1.ton f1.toff f1.ts1 TS2real f1.ka f1.kb f1.ka*(f1.ton-f1.ts1) f1.toff-f1.ton TS2real-f1.ts1 CycleTime(CycleI-10)-f1.toff CycleTime(CycleI-10)-TS2real];
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

YlabelK=["Ton","Toff","TonS","ToffE","Ka","kb","ActiveMax","TimeWindow","TimeWindow2",'End-Toff','End-ToffE'];
ELMean=mean(EL,1);
for ii=1:size(timeExtractG,2)
    figure(Hk1)
    subplot(3,4,ii)
    hold on
    plot(ELMean,timeExtractG(:,ii),'-','marker','o','LineWidth',2)
    Kname=char(YlabelK(ii));
    ylabel(Kname)
    xlabel('EL');
    title([Kname,' vs. EL'])
end


for ii=1:size(timeExtractG,2)
    figure(Hk2)
    subplot(3,4,ii)
    hold on
    [B,TFrm] = rmoutliers(timeExtractG(:,ii));
    plot(ELMean(~TFrm),timeExtractG(~TFrm,ii),'-','marker','o','LineWidth',2)
    Kname=char(YlabelK(ii));
    ylabel(Kname)
    xlabel('EL');
    title([Kname,' vs. EL'])
end
end
figure(Hk1)
legend('Cycle 11','Cycle 12','Cycle 13')
saveas(gcf,[OutputRoadCycleALL,'Fitting K',' vs EL','.fig']);
saveas(gcf,[OutputRoadCycleALL,'Fitting K',' vs EL','.png']);

figure(Hk2)
legend('Cycle 11','Cycle 12','Cycle 13')
saveas(gcf,[OutputRoadCycleALL,'Fitting K',' vs EL',' DelOutliers.fig']);
saveas(gcf,[OutputRoadCycleALL,'Fitting K',' vs EL',' DelOutliers.png']);