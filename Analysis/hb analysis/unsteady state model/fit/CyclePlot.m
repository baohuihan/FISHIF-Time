%%CyclePlot
%% Ini
ModelPlot='Ton';
InheritSwitchPlot=1;
RepeatNum=1;
TimeBinPlot=1;
CombinNum=1;
% LabelTip='-EL[0.2 0.4]bin0.2-DataFilter-ShowInherit--Label[0 0 1]-kfitTime[0 1 1]\';
% LabelTip='-EL[0.2 0.4]bin0.2-DataFilter-ShowInherit--Label[0 0 1]-kfitTime[0 1 1]-KL[0.001 0.001 0.01]-KU[0.02 0.02 0.55]-TestBug3\';
LabelTip='-EL[0 0.6]bin0.2-DataFilter-ShowInherit--Label[0 0-1 1]-kfitTime[0 0.5 1]-KL[0.08 0.08 0.6]-KU[0.1 0.1 0.7]-TimeStart60-KDisLost-fix\';
OutFolderMain=@(Cycle) ['Y:\TimeFit\OutputModel\',ModelPlot,'\Cycle',Cycle,...
    '\Inherit',num2str(InheritSwitchPlot),'-AlpOn-Turn',num2str(RepeatNum),'-TimeBin',num2str(TimeBinPlot),'-Combin',num2str(CombinNum),LabelTip];

OutputRoad=['Y:\TimeFit\OutputModel\',ModelPlot,'\CyclePlot',...
    '\Inherit',num2str(InheritSwitchPlot),'-AlpOn-Turn',num2str(RepeatNum),'-TimeBin',num2str(TimeBinPlot),'-Combin',num2str(CombinNum),LabelTip];mkdir(OutputRoad);
%% Loading
Color=["#19CAAD","#F4606C","#9999CC"];
RoadCycle12=OutFolderMain('12');
RoadCycle13=OutFolderMain('13');
RoadCycle11=OutFolderMain('11');
load([RoadCycle12,'Results.mat'])
Ytime12=TimeLabelMean';
Kplot12=Kplot;
load([RoadCycle13,'Results.mat'])
Ytime13=TimeLabelMean';
Kplot13=Kplot;
load([RoadCycle11,'Results.mat'])
Ytime11=TimeLabelMean';
Kplot11=Kplot;
Xel=mean(EL,1)';

Ytime=[Ytime11;7+Ytime12;9+7+Ytime13];

%% DataMean
try
    figure('Position',[800 600 1000 400])
    for ELi=1:size(EL,2)
        el=EL(:,ELi);
        ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
        nexttile
        hold on
        plot(Ytime,[Kplot11.MeanData(:,ELi);Kplot12.MeanData(:,ELi);Kplot13.MeanData(:,ELi)],'-o','LineWidth',2,'color',Color(1))
        plot(Ytime,[Kplot11.MeanFit(:,ELi);Kplot12.MeanFit(:,ELi);Kplot13.MeanFit(:,ELi)],'-o','LineWidth',2,'color',Color(2))
        plot(Ytime,[Kplot11.MeanFitSteady(:,ELi);Kplot12.MeanFitSteady(:,ELi);Kplot13.MeanFitSteady(:,ELi)],'-o','LineWidth',2,'color',Color(3))
        plot([7 7],ylim,'-k','LineWidth',2)
        plot([16 16],ylim,'-k','LineWidth',2)
        ylabel('MeanExpression')
        xlabel('Time/min');title(ELlabel)
        legend('Data','Fit','FitInSteady')
    end
    pbaspect([3 1 1])
    saveas(gcf,[OutputRoad,'Mean vs Time','.fig']);
    saveas(gcf,[OutputRoad,'Mean vs Time','.png']);
catch
    close gcf
    disp('Not Find: mean information')
end
%% K
try
    figure('Position',[800 600 1000 400])
    for ELi=1:size(EL,2)
        el=EL(:,ELi);
        ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
        nexttile
        hold on
        plot(Ytime,[Kplot11.Kini(:,ELi);Kplot12.Kini(:,ELi);Kplot13.Kini(:,ELi)],'-o','LineWidth',2,'color',Color(1))
        ylabel('Kini');ylim([0 1]);
        yyaxis right
        plot(Ytime,[Kplot11.Alpha(:,ELi);Kplot12.Alpha(:,ELi);Kplot13.Alpha(:,ELi)],'--d','LineWidth',2,'color',Color(2))
        plot([7 7],ylim,'-k','LineWidth',2)
        plot([16 16],ylim,'-k','LineWidth',2)
        ylabel('Alpha');ylim([0 1]);
        xlabel('Time/min');title(ELlabel)
        legend('Kini','Alpha')
    end
    pbaspect([3 1 1])
    saveas(gcf,[OutputRoad,'Kini vs Time','.fig']);
    saveas(gcf,[OutputRoad,'Kini vs Time','.png']);
catch
    close gcf
    disp('Not Find: Kini information')
end

try
    figure('Position',[800 600 1000 400])
    for ELi=1:size(EL,2)
        el=EL(:,ELi);
        ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
        nexttile
        hold on
        plot(Ytime,[Kplot11.Kon(:,ELi);Kplot12.Kon(:,ELi);Kplot13.Kon(:,ELi)],'-o','LineWidth',2,'color',Color(1))
        ylabel('K')
        plot(Ytime,[Kplot11.Koff(:,ELi);Kplot12.Koff(:,ELi);Kplot13.Koff(:,ELi)],'--d','LineWidth',2,'color',Color(2))
        plot([7 7],ylim,'-k','LineWidth',2)
        plot([16 16],ylim,'-k','LineWidth',2)
        xlabel('Time/min');title(ELlabel)
    %     ylim([0 0.1])
        legend('Kon','Koff')
    end
    ylim([0 0.05])
    pbaspect([3 1 1])
    saveas(gcf,[OutputRoad,'K vs Time','.fig']);
    saveas(gcf,[OutputRoad,'K vs Time','.png']);
catch
    close gcf
    disp('Not Find: Kon/off information')
end

try
    figure('Position',[800 600 1000 400])
    for ELi=1:size(EL,2)
        el=EL(:,ELi);
        ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
        nexttile
        hold on
        plot(Ytime,[Kplot11.ON(:,ELi);Kplot12.ON(:,ELi);Kplot13.ON(:,ELi)],'-o','LineWidth',2,'color',Color(1))
        ylabel('K')
        plot(Ytime,[Kplot11.OFF(:,ELi);Kplot12.OFF(:,ELi);Kplot13.OFF(:,ELi)],'--d','LineWidth',2,'color',Color(2))
        plot([7 7],ylim,'-k','LineWidth',2)
        plot([16 16],ylim,'-k','LineWidth',2)
        xlabel('Time/min');title(ELlabel)
    %     ylim([0 0.1])
        legend('ON','OFF')
    end
    pbaspect([3 1 1])
    saveas(gcf,[OutputRoad,'ONOFF vs Time','.fig']);
    saveas(gcf,[OutputRoad,'ONOFF vs Time','.png']);
catch
    close gcf
    disp('Not Find: ON/OFF information')
end
%% Ton
try
    figure('Position',[800 600 1000 400])
    for ELi=1:size(EL,2)
        el=EL(:,ELi);
        ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
        nexttile
        hold on
        plot(Ytime,[Kplot11.ActiveR(:,ELi);Kplot12.ActiveR(:,ELi);Kplot13.ActiveR(:,ELi)],'-o','LineWidth',2,'color',Color(1))
        plot([7 7],ylim,'-k','LineWidth',2)
        plot([16 16],ylim,'-k','LineWidth',2)
        ylabel('Active nucleus ratio')
        xlabel('Time/min');title(ELlabel)
        ylim([0 1])
    end
    pbaspect([3 1 1])
    saveas(gcf,[OutputRoad,'Ton vs Time','.fig']);
    saveas(gcf,[OutputRoad,'Ton vs Time','.png']);
catch
    close gcf
    disp('Not Find: Ton information')
end
%% State
try
    figure('Position',[800 600 1000 400])
    for ELi=1:size(EL,2)
        el=EL(:,ELi);
        ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
        nexttile
        hold on
        plot(Ytime,[Kplot11.Sf(:,ELi);Kplot12.Sf(:,ELi);Kplot13.Sf(:,ELi)],'-o','LineWidth',2,'color',Color(1),'MarkerSize',3)
        plot(Ytime,[Kplot11.S0(:,ELi);Kplot12.S0(:,ELi);Kplot13.S0(:,ELi)],'-o','LineWidth',2,'color',Color(2),'MarkerSize',3)
        plot(Ytime,[Kplot11.S1(:,ELi);Kplot12.S1(:,ELi);Kplot13.S1(:,ELi)],'-o','LineWidth',2,'color',Color(3),'MarkerSize',3)
        plot(Ytime,[Kplot11.Kon(:,ELi)./(Kplot11.Kon(:,ELi)+Kplot11.Koff(:,ELi));...
            Kplot12.Kon(:,ELi)./(Kplot12.Kon(:,ELi)+Kplot12.Koff(:,ELi));...
            Kplot13.Kon(:,ELi)./(Kplot13.Kon(:,ELi)+Kplot13.Koff(:,ELi))],'--','LineWidth',2,'color','#EDB120','Marker','o','MarkerSize',3)
        plot([7 7],ylim,'-k','LineWidth',2)
        plot([16 16],ylim,'-k','LineWidth',2)
        ylabel('State')
        xlabel('Time/min');title(ELlabel)
        ylim([0 1])
        legend({'S0','S1','S2','Steady On'},'Location','northeastoutside')
    end
    pbaspect([3 1 1])
    % daspect([4,2,1])
    saveas(gcf,[OutputRoad,'State vs Time','.fig']);
    saveas(gcf,[OutputRoad,'State vs Time','.png']);
catch
    close gcf
    disp('Not Find: State information')
end
%% 
%% END