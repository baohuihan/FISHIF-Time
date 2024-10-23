%%CyclePlot
OutputRoad='Y:\TimeFit\OutputNew\CDS\CyclePlot\Inherit0-Three-AlpON-Turn10-TimeBin0.5-NewIniModel1-EL[0.2 0.4]bin0.2-Lock[Koff 0.0001 0.1]-Kini[0.15 0.55]-7.5minFN-k20-TimeCombin\';mkdir(OutputRoad);
Color=["#19CAAD","#F4606C","#9999CC"];
RoadCycle12='Y:\TimeFit\OutputNew\CDS\Cycle12\Inherit0-Three-AlpON-Turn10-TimeBin0.5-NewIniModel1-EL[0.2 0.4]bin0.2-Lock[Koff 0.0001 0.1]-Kini[0.15 0.55]-7.5minFN-k20-TimeCombin\';
RoadCycle13='Y:\TimeFit\OutputNew\CDS\Cycle13\Inherit0-Three-AlpON-Turn10-TimeBin0.5-NewIniModel1-EL[0.2 0.4]bin0.2-Lock[Koff 0.0001 0.1]-7.5minFN-k20-TimeCombin\';
RoadCycle11='Y:\TimeFit\OutputNew\CDS\Cycle11\Inherit0-Three-AlpON-Turn10-TimeBin0.5-NewIniModel1-EL[0.2 0.4]bin0.2-Lock[Koff 0.0001 0.1]-Kini[0.15 0.55]-7.5minFN-k20-TimeCombin\';
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

figure%2D
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    nexttile
    hold on
    plot(Ytime,[Kplot11.MeanData(:,ELi);Kplot12.MeanData(:,ELi);Kplot13.MeanData(:,ELi)],'-o','LineWidth',2,'color',Color(1))
    plot(Ytime,[Kplot11.MeanFit(:,ELi);Kplot12.MeanFit(:,ELi);Kplot13.MeanFit(:,ELi)],'-o','LineWidth',2,'color',Color(2))
    plot(Ytime,[Kplot11.MeanFitSteady(:,ELi);Kplot12.MeanFitSteady(:,ELi);Kplot13.MeanFitSteady(:,ELi)],'-o','LineWidth',2,'color',Color(3))
    ylabel('MeanExpression')
    xlabel('Time/min');title(ELlabel)
    legend('Data','Fit','FitInSteady')
end
saveas(gcf,[OutputRoad,'Mean vs Time','.fig']);
saveas(gcf,[OutputRoad,'Mean vs Time','.png']);

figure%2D
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    nexttile
    hold on
    plot(Ytime,[Kplot11.Kini(:,ELi);Kplot12.Kini(:,ELi);Kplot13.Kini(:,ELi)],'-o','LineWidth',2,'color',Color(1))
    ylabel('Kini');ylim([0 1]);
    yyaxis right
    plot(Ytime,[Kplot11.Alpha(:,ELi);Kplot12.Alpha(:,ELi);Kplot13.Alpha(:,ELi)],'--d','LineWidth',2,'color',Color(2))
    ylabel('Alpha');ylim([0 1]);
    xlabel('Time/min');title(ELlabel)
    legend('Kini','Alpha')
end
saveas(gcf,[OutputRoad,'Kini vs Time','.fig']);
saveas(gcf,[OutputRoad,'Kini vs Time','.png']);


figure%2D
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    nexttile
    hold on
    plot(Ytime,[Kplot11.Kon(:,ELi);Kplot12.Kon(:,ELi);Kplot13.Kon(:,ELi)],'-o','LineWidth',2,'color',Color(1))
    ylabel('K')
    plot(Ytime,[Kplot11.Koff(:,ELi);Kplot12.Koff(:,ELi);Kplot13.Koff(:,ELi)],'--d','LineWidth',2,'color',Color(2))
    xlabel('Time/min');title(ELlabel)
%     ylim([0 0.1])
    legend('Kon','Koff')
end
saveas(gcf,[OutputRoad,'K vs Time','.fig']);
saveas(gcf,[OutputRoad,'K vs Time','.png']);

figure%2D
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    nexttile
    hold on
    plot(Ytime,[Kplot11.ON(:,ELi);Kplot12.ON(:,ELi);Kplot13.ON(:,ELi)],'-o','LineWidth',2,'color',Color(1))
    ylabel('K')
    plot(Ytime,[Kplot11.OFF(:,ELi);Kplot12.OFF(:,ELi);Kplot13.OFF(:,ELi)],'--d','LineWidth',2,'color',Color(2))
    xlabel('Time/min');title(ELlabel)
%     ylim([0 0.1])
    legend('ON','OFF')
end
saveas(gcf,[OutputRoad,'ONOFF vs Time','.fig']);
saveas(gcf,[OutputRoad,'ONOFF vs Time','.png']);

% figure%2D
% for ELi=1:size(EL,2)
%     el=EL(:,ELi);
%     ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
%     nexttile
%     hold on
%     plot(Ytime,[Kplot12.Ton(:,ELi);Kplot13.Ton(:,ELi)],'-o','LineWidth',2,'color',Color(1))
%     ylabel('Ton')
%     xlabel('Time/min');title(ELlabel)
%     ylim([0 1])
% end
% saveas(gcf,[OutputRoad,'Ton vs Time','.fig']);
% saveas(gcf,[OutputRoad,'Ton vs Time','.png']);
