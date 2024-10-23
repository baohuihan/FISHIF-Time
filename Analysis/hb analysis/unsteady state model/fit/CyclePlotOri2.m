%%CyclePlotOri2
%% Ini
ModelPlot='Ton';
InheritSwitchPlot=1;
RepeatNum=1;
TimeBinPlot=0.5;
CombinNum=3;
PonLock=0;

LabelTip='-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-Vel25-KL[0.01 0.01 0.60]-KU[0.2 0.2 0.8]-TimeStart60-KDisLost[0 0 0 0 0]-Artifical-Filter5-timeBinMove1-PonBoole0-re4\';
Cycle=[11];

kName=["k01" "k02" "k10" "k12" "k20" "k21" "ki0" "ki1" "ki2" "Ton" "dataintensity1" "alpha0"];
OutFolderMain=@(Cycle) ['Y:\TimeFit\OutputModel\',ModelPlot,'\Cycle',Cycle,...
    '\Inherit',num2str(InheritSwitchPlot),'-AlpOn-Turn',num2str(RepeatNum),'-TimeBin',num2str(TimeBinPlot),'-Combin',num2str(CombinNum),'-PonLock',num2str(PonLock),LabelTip];
OutputRoad=['Y:\TimeFit\OutputModel\',ModelPlot,'\CyclePlot',...
    '\Inherit',num2str(InheritSwitchPlot),'-AlpOn-Turn',num2str(RepeatNum),'-TimeBin',num2str(TimeBinPlot),'-Combin',num2str(CombinNum),'-PonLock',num2str(PonLock),LabelTip];mkdir(OutputRoad);
%% Loading
Color=["#19CAAD","#F4606C","#9999CC"];

CycleTime=[0,7,9+7];
DivisionTime=[1,1,1.5];
Ytime=CycleTime+DivisionTime;

h1=figure('Position',[800 600 1000 400]);
h2=figure('Position',[800 600 1000 400]);h21=figure('Position',[800 600 1000 400]);
h3=figure('Position',[800 600 1000 400]);
h4=figure('Position',[800 600 1000 400]);
h5=figure('Position',[800 600 1000 400]);
h6=figure('Position',[800 600 1000 400]);
h7=figure('Position',[800 600 1000 400]);
h8=figure('Position',[800 600 1000 400]);

for CycleI=Cycle
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
    LldG=LldPlotMeanELs{ELi,1};
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
        plot(TimeLabelMean,KannealG(:,5),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        plot(TimeLabelMean,KannealG(:,6),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Kini');ylim([0 2]);
        yyaxis right
        plot(TimeLabelMean,KannealG(:,end),'-.','LineWidth',1,'color',c_map(ELi,:))
       
        ylabel('Alpha');ylim([0 1]);
        xlabel('Time/sec');title('Kini&Alpha vs. Time')
        % legend('Kini','Alpha')
catch
    close gcf
    disp('Not Find: Kini information')
end

KonChoose=KannealG(:,4)>=max(KannealG(:,4),[],1)|[diff(KannealG(:,4));0]>0;

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


try
    figure(h6)
        hold on
        plot(TimeLabelMean,KannealG(:,3),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Kon')
 
        xlabel('Time/sec');title('Kon 12 vs. Time')
        ylim([0 0.2])
        % legend('Kon','Koff')
catch
    close gcf
    disp('Not Find: Kon/off information')
end



try
    figure(h7)
        hold on
        plot(TimeLabelMean,KannealG(:,4),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Koff')

        xlabel('Time/sec');title('Koff 21 vs. Time')
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
        plot(TimeLabelMean,KannealG(:,7),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
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
    subplot(2,1,1)
        hold on
        plot(TimeLabelMean,PiniG(:,2),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('S1 state frequence')
        xlabel('Time/sec');title('S1 vs. Time')
        ylim([0 0.5])
    subplot(2,1,2)
        hold on
        plot(TimeLabelMean,PiniG(:,3),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('S2 state frequence')
        xlabel('Time/sec');title('S2 vs. Time')
        ylim([0 0.5])
catch
    close gcf
    disp('Not Find: Son information')
end

%% Lost
try
    figure(h8)
    subplot(2,1,1)
        hold on
        plot(TimeLabelMean,LldG(:,1),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Likelihood state frequence')
        xlabel('Time/sec');title('Likelihood vs. Time')
        % ylim([0 1])
catch
    close gcf
    disp('Not Find: Likelihood information')
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
% legend(ELlegends,'Location','northeastoutside')
pbaspect([6 1 1])
saveas(gcf,[OutputRoad,'Son vs Time','.fig']);
saveas(gcf,[OutputRoad,'Son vs Time','.png']);

figure(h6)
plot([7 7]*60,ylim,'-k','LineWidth',2)
plot([16 16]*60,ylim,'-k','LineWidth',2)
legend(ELlegends,'Location','northeastoutside')
pbaspect([3 1 1])
saveas(gcf,[OutputRoad,'k12 vs Time','.fig']);
saveas(gcf,[OutputRoad,'k12 vs Time','.png']);
figure(h7)
plot([7 7]*60,ylim,'-k','LineWidth',2)
plot([16 16]*60,ylim,'-k','LineWidth',2)
legend(ELlegends,'Location','northeastoutside')
pbaspect([3 1 1])
saveas(gcf,[OutputRoad,'k21 vs Time','.fig']);
saveas(gcf,[OutputRoad,'k21 vs Time','.png']);
figure(h8)
plot([7 7]*60,ylim,'-k','LineWidth',2)
plot([16 16]*60,ylim,'-k','LineWidth',2)
legend(ELlegends,'Location','northeastoutside')
pbaspect([3 1 1])
saveas(gcf,[OutputRoad,'lld vs Time','.fig']);
saveas(gcf,[OutputRoad,'lld vs Time','.png']);
%% END