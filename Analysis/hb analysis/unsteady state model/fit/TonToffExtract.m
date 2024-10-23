%%TonToffExtract
% 'Inherit1-AlpOFF-Turn1-TimeBin1-Combin0-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.05 0.05 0.1]-KU[0.05 0.05 0.6]-Alp0.6\'
%Inherit1-Alp0-Turn1-TimeBin1-Combin1-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.01 0.6]-KU[0.10 0.10 0.7]-TimeStart60-KDisLost-fix
FoldName='Inherit1-AlpOn-Turn1-TimeBin0.5-Combin2-EL[0 0.6]bin0.2Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.60]-KU[0.10 0.08 0.75]-TimeStart60-KDisLost-CombinPartArtifical\';
%'Inherit1-AlpOn-Turn1-TimeBin0.5-Combin2-EL[0 0.6]bin0.2Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.6]-KU[0.10 0.08 0.7]-TimeStart60-KDisLost-CombinPartArtifical-fix\';
%'Inherit1-AlpOn-Turn1-TimeBin1-Combin2-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.6]-KU[0.10 0.08 0.7]-TimeStart60-KDisLost\';
%% Fix Time
kfitTime=(3+283+2277/2)/25;Tstay=sum([3+283 2277 407+509])/25;
FixTime=0;
%%
Hk1=figure('Position',[800 600 1000 400]);
pbaspect([3 1 1])
for CycleI=11:13
    CycleLabel=num2str(CycleI);
OutputRoadFun=@(x) ['Y:\TimeFit\OutputModel\Ton\Cycle',x,'\',...
    FoldName];
OutputRoad=OutputRoadFun(CycleLabel);
OutputRoadCycleALL=OutputRoadFun('Plot');
load([OutputRoad,'Results0.mat'])
switch CycleLabel
    case '11'
        startPoint=[0.008 -0.008 240 120 -20 320];%CYCLE11
        Ub=[inf 0 260 200 20 400];
        Lb=[0 -inf 200 50 -50 300];
        timeSeries=-100:1:500;
    case '12'
        startPoint=[0.008 -0.008 290 120 -20 380];%CYCLE12
        timeSeries=-100:1:500;
        Ub=[inf 0 310 140 50 400];
        Lb=[0 -inf 270 50 -20 360];
    case '13'
        startPoint=[0.0035 -0.0035 530 80 0 650];%CYCLE13
        Ub=[inf 0 580 120 50 700];
        Lb=[0 -inf 480 50 -100 600];
        timeSeries=-100:1:800;
end
    
Hk=figure('Position',[800 600 1000 400]);
pbaspect([3 1 1])
% funFit=@(ton,ts1,ts2,toff,rs,ka,kb,x) piecewise(x<=ton,0,x>ton&x<ts1,ka*(x-ton),...
%     x>=ts1&x<=ts2,rs,x>ts2&x<toff,rs+kb*(x-ts2),x>=toff,0);
ft = fittype('piecewiseLine(x,ton,ts1,ts2,toff,ka,kb)');
% 'independent', {'x'},'problem', {'ton','ts1','ts2','toff','ka','kb'}, ...
%  'dependent','y');
c_map = [0.57, 0.69, 0.30
     0.89, 0.88, 0.57
     0.76, 0.49, 0.58
     0.47, 0.76, 0.81
     0.21, 0.21, 0.35
     0.28, 0.57, 0.54
     0.07, 0.35, 0.40
     0.41, 0.20, 0.42
     0.60, 0.24, 0.18
     0.76, 0.84, 0.65];

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
    % 
    ActiveR=KannealG(:,4);
    f1 = fit(TimeLabelMean',ActiveR,ft,...
        'startpoint',startPoint,'upper',Ub,'lower',Lb);
    % ,...
        % 'Lower',[0.0025 -0.0055 400 60 -180 600],...
        % 'Upper',[0.0045 -0.0025 600 90 -130 900]);
        TS2real=(f1.ka*(f1.ton-f1.ts1)+f1.kb*f1.toff)/f1.kb;
    timeExtractG(ELi,:)=[f1.ton f1.toff f1.ts1 TS2real f1.ka f1.kb f1.ka*(f1.ton-f1.ts1)];
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

YlabelK=["Ton","Toff","TonS","ToffE","Ka","kb","ActiveMax"];
ELMean=mean(EL,1);
for ii=1:size(timeExtractG,2)
    figure(Hk1)
    subplot(2,4,ii)
    hold on
    plot(ELMean,timeExtractG(:,ii),'-o','LineWidth',2)
    Kname=char(YlabelK(ii));
    ylabel(Kname)
    xlabel('EL');
    title([Kname,' vs. EL'])
end
end
legend('Cycle 11','Cycle 12','Cycle 13')

saveas(gcf,[OutputRoadCycleALL,'Fitting K',' vs EL','.fig']);
saveas(gcf,[OutputRoadCycleALL,'Fitting K',' vs EL','.png']);