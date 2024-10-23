%%CyclePlotOri
%% Ini
ModelPlot='Ton';
InheritSwitchPlot=1;
RepeatNum=1;
TimeBinPlot=0.5;
CombinNum=3;
PonLock=1;
LabelTip='-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-Vel40-KL[0.01 0.11 1]-KU[0.2 0.13 1.3]-KDisLost[10 15 0 0 0]FirstON-Artifical-Filter6-timeBinMove1-PonBoole0-LowArtK-Star0708\';
%'-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-Vel40-KL[0.01 0.12 1]-KU[0.2 0.12 1.3]-KDisLost[5 5 1 1 1]FOFF-Artifical-Filter5-timeBinMove1-PonBoole0-re6-0523-Star3\';
% '-TimeWindowMingle-PonBoole0-0512\';
FoldNameCycle11=...
    '-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-Vel40-KL[0.01 0.12 1]-KU[0.2 0.12 1.3]-TimeStart60-KDisLost[10 10 1 1 1]-Artifical-Filter5-timeBinMove1-PonBoole0-re6-0521-Star2\';
FoldNameCycle12=...
    '-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-Vel40-KL[0.01 0.12 1]-KU[0.2 0.12 1.4]-KDisLost[10 5 0 0 0]FirstON-Artifical-Filter6-timeBinMove1-PonBoole0-Star0713\';
FoldNameCycle13=...
    '-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-Vel40-KL[0.01 0.12 1]-KU[0.2 0.12 1.3]-TimeStart60-KDisLost[10 10 1 1 1]-Artifical-Filter5-timeBinMove1-PonBoole0-re6-0521-Star2\';

LabelTipSave='-TimeWindowMingle-PonBoole0-0726-bcd-2d\';
CycleTime=[0,7,9+7];
DivisionTime=[1,1,1]*60;
CycleTimePlot=[0,7+1,7+1+9+1]*60;
% Ytime=CycleTime+cumsum(DivisionTime);
Ytime=[0 0 0];
%'-EL[0 0.6]bin0.1Move0.05-TimeWindowAdd-Vel45-KL[0.01 0.1 0.60]-KU[0.2 0.1 0.8]-TimeStart60-KDisLost[10 10 1 1 1]-Artifical-Filter5-timeBinMove1-PonBoole0-re4\';
%'-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-Vel50-KL[0.01 0.1 0.60]-KU[0.2 0.1 0.8]-TimeStart60-KDisLost[25 10 1 1 5]-Artifical-Filter5-timeBinMove1-PonBoole0-re4-Star\';
%'-EL[0 0.6]bin0.1Move0.05-DataFilter-KL[0.01 0.12 0.60]-KU[0.2 0.12 0.75]-TimeStart60-KDisLost[10 10 1 1 5]-Artifical-Filter5-timeBinMove1-PonBoole-re4\';
% '-EL[0 0.6]bin0.1Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[40 10 1 1 5]-Artifical-Filter4-timeBinMove1-re\';
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
%'-EL[0 0.6]bin0.1Move0.1-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.05 0.08 0.60]-KU[0.1 0.08 0.75]-TimeStart60-KDisLost[50 10 1 1]-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.1Move0.1-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.0001 0.0001 0.60]-KU[0.0020 0.0020 0.75]-TimeStart60-KDisLost-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.1Move0.1-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.001 0.001 0.60]-KU[0.010 0.010 0.75]-TimeStart60-KDisLost-CombinPartArtifical-artificalFilter2\';
%'-EL[0 0.6]bin0.2Move0.05-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.01 0.08 0.60]-KU[0.10 0.08 0.75]-TimeStart60-KDisLost-CombinPartArtifical\';
%'-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.08 0.08 0.6]-KU[0.1 0.1 0.7]-TimeStart60-KDisLost-fix\';
OutFolderMain=@(Cycle) ['Y:\TimeFit\OutputModel\',ModelPlot,'\Cycle',Cycle,...
    '\Inherit',num2str(InheritSwitchPlot),'-AlpOn-Turn',num2str(RepeatNum),'-TimeBin',num2str(TimeBinPlot),'-Combin',num2str(CombinNum),'-PonLock',num2str(PonLock),LabelTip];
OutputRoad=['Y:\TimeFit\OutputModel\',ModelPlot,'\CyclePlot',...
    '\Inherit',num2str(InheritSwitchPlot),'-AlpOn-Turn',num2str(RepeatNum),'-TimeBin',num2str(TimeBinPlot),'-Combin',num2str(CombinNum),'-PonLock',num2str(PonLock),LabelTipSave];mkdir(OutputRoad);
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



h1=figure('Position',[800 600 1000 400]);
h2=figure('Position',[800 600 1000 400]);h21=figure('Position',[800 600 1000 400]);
h3=figure('Position',[800 600 1000 400]);
h4=figure('Position',[800 600 1000 400]);
h5=figure('Position',[800 600 1000 400]);

XTimeAll=[];
KGAAll=[];KoffGAAll=[];
TonGAAll=[];
KiniGAAll=[];

for CycleI=11:13
    KGA=[];TGA=[];TonGA=[];KiniGA=[];KoffGA=[];
    CycleLabel=num2str(CycleI);
    switch CycleLabel
        case '11'
            if ~isempty(FoldNameCycle11)
                LabelTip = FoldNameCycle11;
            end
        case '12'
            if ~isempty(FoldNameCycle12)
                LabelTip = FoldNameCycle12;
            end
        case '13'
            if ~isempty(FoldNameCycle13)
                LabelTip = FoldNameCycle13;
            end
    end
    OutFolderMain=@(Cycle) ['Y:\TimeFit\OutputModel\',ModelPlot,'\Cycle',Cycle,...
    '\Inherit',num2str(InheritSwitchPlot),'-AlpOn-Turn',num2str(RepeatNum),'-TimeBin',num2str(TimeBinPlot),'-Combin',num2str(CombinNum),'-PonLock',num2str(PonLock),LabelTip];

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
    TimeLabelMean=TimePlotMeanELs{ELi,:}+Ytime(CycleI-10);

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

KonChoose=KannealG(:,4)+1e-3>=max(KannealG(:,4),[],1)|[diff(KannealG(:,4));0]>0;

try
    figure(h2)
        hold on
        plot(TimeLabelMean,KannealG(:,1),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Kon')
 
        xlabel('Time/sec');title('Kon vs. Time')
        ylim([0 0.2])
        % legend('Kon','Koff')
    
    % ylim([0 0.05])
    figure(h21)
        hold on
        plot(TimeLabelMean(KonChoose),KannealG(KonChoose,1),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Kon credible')
 
        xlabel('Time/sec');title('Kon vs. Time')
        ylim([0 0.2])
catch
    close gcf
    disp('Not Find: Kon/off information')
end

if 1%CycleI==12
    % KGA=[KGA,KannealG(KonChoose,1)];
    % TGA=[TGA,TimeLabelMean(KonChoose)'];
    KGA=[KGA,KannealG(:,1)];
    KoffGA=[KoffGA,KannealG(:,2)];
    KiniGA=[KiniGA,KannealG(:,3)];
    TonGA=[TonGA,KannealG(:,4)];
    TGA=[TGA,TimeLabelMean'];
end


try
    figure(h3)
        hold on
        plot(TimeLabelMean,KannealG(:,2),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Koff')

        xlabel('Time/sec');title('Koff vs. Time')
        ylim([0 0.2])
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
figure('Position',[800 600 800 600]);
XTime=mean(TGA,2);
YEL=mean(EL,1);
imagesc(XTime',YEL,KGA')
c=colorbar;
c.Label.String = 'Kon';
axis square
ylabel('EL')
xlabel('Time/sec');title('2D Kon vs. Time&EL')
saveas(gcf,[OutputRoad,'Cycle',CycleLabel,' Kon vs Time&EL','.fig']);
saveas(gcf,[OutputRoad,'Cycle',CycleLabel,' Kon vs Time&EL','.png']);

figure('Position',[800 600 800 600]);
XTime=mean(TGA,2);
YEL=mean(EL,1);
imagesc(XTime',YEL,KoffGA')
c=colorbar;
c.Label.String = 'Koff';
clim([0.1 0.15]);
axis square
ylabel('EL')
xlabel('Time/sec');title('2D Koff vs. Time&EL')
saveas(gcf,[OutputRoad,'Cycle',CycleLabel,' Koff vs Time&EL','.fig']);
saveas(gcf,[OutputRoad,'Cycle',CycleLabel,' Koff vs Time&EL','.png']);

figure('Position',[800 600 800 600]);
XTime=mean(TGA,2);
YEL=mean(EL,1);
imagesc(XTime',YEL,KiniGA')
c=colorbar;
c.Label.String = 'Kini';
clim([0.5 1.5]);
axis square
ylabel('EL')
xlabel('Time/sec');title('2D Kini vs. Time&EL')
saveas(gcf,[OutputRoad,'Cycle',CycleLabel,' Kini vs Time&EL','.fig']);
saveas(gcf,[OutputRoad,'Cycle',CycleLabel,' Kini vs Time&EL','.png']);

figure('Position',[800 600 800 600]);
XTime=mean(TGA,2);
YEL=mean(EL,1);
imagesc(XTime',YEL,TonGA')
c=colorbar;
c.Label.String = 'Pon';
axis square
ylabel('EL')
xlabel('Time/sec');title('2D Pon vs. Time&EL')
saveas(gcf,[OutputRoad,'Cycle',CycleLabel,' Pon vs Time&EL','.fig']);
saveas(gcf,[OutputRoad,'Cycle',CycleLabel,' Pon vs Time&EL','.png']);
save([OutputRoad,'Results',CycleLabel,'.mat'],'KGA','KoffGA','TonGA','KiniGA','XTime','YEL');

XTimeAll=[XTimeAll;XTime+CycleTimePlot(CycleI-10)+DivisionTime(CycleI-10)];
KGAAll=[KGAAll;KGA];
KoffGAAll=[KoffGAAll;KoffGA];
TonGAAll=[TonGAAll;TonGA];
KiniGAAll=[KiniGAAll;KiniGA];
end
%% 
figure('Position',[800 600 1400 600]);
XTime=mean(TGA,2);
YEL=mean(EL,1);
imagesc(XTimeAll',YEL,KGAAll')
c=colorbar;
c.Label.String = 'Kon';
clim([0 0.15]);
% axis square
ylabel('EL')
xlabel('Time/sec');title('2D Kon vs. Time&EL')
saveas(gcf,[OutputRoad,'CycleAll',' Kon vs Time&EL','.fig']);
saveas(gcf,[OutputRoad,'CycleAll',' Kon vs Time&EL','.png']);

figure('Position',[800 600 1400 600]);
XTime=mean(TGA,2);
YEL=mean(EL,1);
imagesc(XTimeAll',YEL,KoffGAAll')
c=colorbar;
c.Label.String = 'Koff';
clim([0.1 0.15]);
% axis square
ylabel('EL')
xlabel('Time/sec');title('2D Koff vs. Time&EL')
saveas(gcf,[OutputRoad,'CycleAll',' Koff vs Time&EL','.fig']);
saveas(gcf,[OutputRoad,'CycleAll',' Koff vs Time&EL','.png']);

figure('Position',[800 600 1400 600]);
XTime=mean(TGA,2);
YEL=mean(EL,1);
imagesc(XTimeAll',YEL,KiniGAAll')
c=colorbar;
c.Label.String = 'Kini';
clim([0.5 1.5]);
% axis square
ylabel('EL')
xlabel('Time/sec');title('2D Kini vs. Time&EL')
saveas(gcf,[OutputRoad,'CycleAll',' Kini vs Time&EL','.fig']);
saveas(gcf,[OutputRoad,'CycleAll',' Kini vs Time&EL','.png']);

figure('Position',[800 600 1400 600]);
XTime=mean(TGA,2);
YEL=mean(EL,1);
imagesc(XTimeAll',YEL,TonGAAll')
c=colorbar;
c.Label.String = 'Pon';
% axis square
ylabel('EL')
xlabel('Time/sec');title('2D Pon vs. Time&EL')
saveas(gcf,[OutputRoad,'CycleAll',' Pon vs Time&EL','.fig']);
saveas(gcf,[OutputRoad,'CycleAll',' Pon vs Time&EL','.png']);
save([OutputRoad,'Results','CycleAll','.mat'],'KGAAll','KoffGAAll','TonGAAll','KiniGAAll','XTimeAll','YEL');



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
figure(h21)
plot([7 7]*60,ylim,'-k','LineWidth',2)
plot([16 16]*60,ylim,'-k','LineWidth',2)
legend(ELlegends,'Location','northeastoutside')
pbaspect([3 1 1])
saveas(gcf,[OutputRoad,'Kon credible vs Time','.fig']);
saveas(gcf,[OutputRoad,'Kon credible vs Time','.png']);
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