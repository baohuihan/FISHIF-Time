%%CyclePlotUse
%% Ini
ModelPlot='Ton';
InheritSwitchPlot=1;
RepeatNum=1;
TimeBinPlot=0.5;
CombinNum=3;
LabelTipSave='Use0417-ELcat\';
LabelTip='Use\';
%'-EL[0 0.6]bin0.1Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[40 10 1 1 5]-CombinPartArtifical-artificalFilter2-new\';
%'-EL[0 0.6]bin0.1Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[40 10 1 1 5]-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.1Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[25 10 1 1 5]-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.1Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.001 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[25 10 1 2]-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.1Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.001 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[25 10 1 1]-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.1Move0.1-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.001 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[25 10 1 1]-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.1Move0.1-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.05 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[25 10 1 1]-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.2Move0.1-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.05 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[25 10 1 1]-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.1Move0.1-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.05 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[25 10 1 1]-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.1Move0.1-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.05 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[50 10 1 1]-CombinPartArtifical-artificalFilter2\';
% '-EL[0 0.6]bin0.1Move0.1-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.05 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[50 10 1 1]-CombinPartArtifical-artificalFilter2\';
% '-EL[0 0.6]bin0.1Move0.1-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.0001 0.0001 0.60]-KU[0.0020 0.0020 0.75]-TimeStart60-KDisLost-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.1Move0.1-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.001 0.001 0.60]-KU[0.010 0.010 0.75]-TimeStart60-KDisLost-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.2Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.60]-KU[0.10 0.08 0.75]-TimeStart60-KDisLost-CombinPartArtifical\';
%'-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.08 0.08 0.6]-KU[0.1 0.1 0.7]-TimeStart60-KDisLost-fix\';
OutFolderMain=@(Cycle) ['Y:\TimeFit\OutputModel\',ModelPlot,'\Cycle',Cycle,...
    '\',LabelTip];
OutputRoad=['Y:\TimeFit\OutputModel\',ModelPlot,'\CyclePlot',...
    '\',LabelTipSave];mkdir(OutputRoad);
%% Loading
Color=["#19CAAD","#F4606C","#9999CC"];

% c_map = [0.57, 0.69, 0.30
%      0.89, 0.88, 0.57
%      0.76, 0.49, 0.58
%      0.47, 0.76, 0.81
%      0.21, 0.21, 0.35
%      0.28, 0.57, 0.54
%      0.07, 0.35, 0.40
%      0.41, 0.20, 0.42
%      0.60, 0.24, 0.18
%      0.76, 0.84, 0.65
%      0.46, 0.54, 0.15
%      0.16, 0.24, 0.85];

CycleTime=[0,7,9+7];
DivisionTime=[1,1,1.5];
Ytime=CycleTime+DivisionTime;

h1=figure('Position',[800 600 1000 400]);
h2=figure('Position',[800 600 1000 400]);
h3=figure('Position',[800 600 1000 400]);
h4=figure('Position',[800 600 1000 400]);
h5=figure('Position',[800 600 1000 400]);

for CycleI=11:13
    CycleLabel=num2str(CycleI);
    RoadCycle=OutFolderMain(CycleLabel);
    load([RoadCycle,'Results0.mat'])
ELlegends={};
c_map = [cool(ceil(size(EL,2)/2));spring(floor(size(EL,2)/2))];
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    
    ELlabels{ELi}=['EL ',num2str(el(1)),'-',num2str(el(2))];
    ELlegends=[ELlegends;['EL ',num2str(el(1)),'-',num2str(el(2))]];
    KannealG=KPlotannealELs{ELi,1};
    PiniG=PstatePlotMeanELs{ELi,1};
    TimeLabelMean=TimePlotMeanELs{ELi,:}+Ytime(CycleI-10)*60;

%% DataMean
% try
%     figure('Position',[800 600 1000 400])
%     for ELi=1:size(EL,2)
%         el=EL(:,ELi);
%         ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
%         nexttile
%         hold on
%         plot(Ytime,[Kplot11.MeanData(:,ELi);Kplot12.MeanData(:,ELi);Kplot13.MeanData(:,ELi)],'-o','LineWidth',2,'color',Color(1))
%         plot(Ytime,[Kplot11.MeanFit(:,ELi);Kplot12.MeanFit(:,ELi);Kplot13.MeanFit(:,ELi)],'-o','LineWidth',2,'color',Color(2))
%         plot(Ytime,[Kplot11.MeanFitSteady(:,ELi);Kplot12.MeanFitSteady(:,ELi);Kplot13.MeanFitSteady(:,ELi)],'-o','LineWidth',2,'color',Color(3))
%         plot([7 7],ylim,'-k','LineWidth',2)
%         plot([16 16],ylim,'-k','LineWidth',2)
%         ylabel('MeanExpression')
%         xlabel('Time/min');title(ELlabel)
%         legend('Data','Fit','FitInSteady')
%     end
%     pbaspect([3 1 1])
%     saveas(gcf,[OutputRoad,'Mean vs Time','.fig']);
%     saveas(gcf,[OutputRoad,'Mean vs Time','.png']);
% catch
%     close gcf
%     disp('Not Find: mean information')
% end
%% K
try
    figure(h1)
        hold on
        plot(TimeLabelMean,KannealG(:,3),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Kini');ylim([0 1]);
        yyaxis right
        plot(TimeLabelMean,KannealG(:,5),'-.','LineWidth',1,'color',c_map(ELi,:))
       
        ylabel('Alpha');ylim([0 1]);
        xlabel('Time/sec');title('Kini&Alpha vs. Time')
        % legend('Kini','Alpha')
catch
    close gcf
    disp('Not Find: Kini information')
end

try
    figure(h2)
        hold on
        plot(TimeLabelMean,KannealG(:,1),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Kon')
 
        xlabel('Time/sec');title('Kon vs. Time')
        ylim([0 0.1])
        % legend('Kon','Koff')
    
    % ylim([0 0.05])
    
catch
    close gcf
    disp('Not Find: Kon/off information')
end

try
    figure(h3)
        hold on
        plot(TimeLabelMean,KannealG(:,2),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Koff')

        xlabel('Time/sec');title('Koff vs. Time')
        ylim([0 0.1])
        % legend('Kon','Koff')
    
    % ylim([0 0.05])
    
catch
    close gcf
    disp('Not Find: Kon/off information')
end
%% Ton
try
    figure(h4)
        hold on
        plot(TimeLabelMean,KannealG(:,4),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Active nucleus ratio')
        xlabel('Time/sec');title('Pon vs. Time')
        ylim([0 1])
catch
    close gcf
    disp('Not Find: Ton information')
end
%% State
try
    figure(h5)
        hold on
        plot(TimeLabelMean,PiniG(:,2),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('ON state frequence')
        xlabel('Time/sec');title('Son vs. Time')
        ylim([0 1])
catch
    close gcf
    disp('Not Find: Son information')
end

end
end
%% 
figure(h1)
plot([7 7]*60,ylim,'-k','LineWidth',2)
plot([16 16]*60,ylim,'-k','LineWidth',2)
legend(ELlegends,'Location','northeastoutside')
pbaspect([3 1 1])
saveas(gcf,[OutputRoad,'Kini vs Time','.fig']);
saveas(gcf,[OutputRoad,'Kini vs Time','.png']);
figure(h2)
plot([7 7]*60,ylim,'-k','LineWidth',2)
plot([16 16]*60,ylim,'-k','LineWidth',2)
legend(ELlegends,'Location','northeastoutside')
pbaspect([3 1 1])
saveas(gcf,[OutputRoad,'Kon vs Time','.fig']);
saveas(gcf,[OutputRoad,'Kon vs Time','.png']);
figure(h3)
plot([7 7]*60,ylim,'-k','LineWidth',2)
plot([16 16]*60,ylim,'-k','LineWidth',2)
legend(ELlegends,'Location','northeastoutside')
pbaspect([3 1 1])
saveas(gcf,[OutputRoad,'Koff vs Time','.fig']);
saveas(gcf,[OutputRoad,'Koff vs Time','.png']);
figure(h4)
plot([7 7]*60,ylim,'-k','LineWidth',2)
plot([16 16]*60,ylim,'-k','LineWidth',2)
legend(ELlegends,'Location','northeastoutside')
pbaspect([3 1 1])
saveas(gcf,[OutputRoad,'Ton vs Time','.fig']);
saveas(gcf,[OutputRoad,'Ton vs Time','.png']);
figure(h5)
plot([7 7]*60,ylim,'-k','LineWidth',2)
plot([16 16]*60,ylim,'-k','LineWidth',2)
legend(ELlegends,'Location','northeastoutside')
pbaspect([3 1 1])
saveas(gcf,[OutputRoad,'Son vs Time','.fig']);
saveas(gcf,[OutputRoad,'Son vs Time','.png']);
%% END