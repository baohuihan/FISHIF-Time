function ModuleFitTimeSaveTwoState(TimeLabelMean,EL,KannealELs,PiniELs,MeanExpELs,OutFolderMain,Color)
Kplot=struct('ON',[],'OFF',[],'Kon',[],'Koff',[],'Sf',[],'S0',[],'S1',[],'Kratio',[],'Kini',[],'Alpha',[],'MeanData',[],'MeanFit',[],'MeanFitSteady',[],'OnRatio',[]);
for ELi=1:size(EL,2)
    KannealG=KannealELs{ELi,1};
    PiniG=PiniELs{ELi,1};
    MeanExp=MeanExpELs{ELi,:};
    Kplot.Kon=[Kplot.Kon,KannealG(:,1)];
    Kplot.Koff=[Kplot.Koff,KannealG(:,2)];
    Kplot.Sf=[Kplot.Sf,PiniG(:,1)];Kplot.S0=[Kplot.S0,PiniG(:,2)];Kplot.S1=[Kplot.S1,PiniG(:,3)];
    Kplot.Kratio=[Kplot.Kratio,KannealG(:,1)./KannealG(:,2)];
    Kplot.OnRatio=[Kplot.OnRatio,KannealG(:,1)./(KannealG(:,1)+KannealG(:,2))];
    Kplot.Kini=[Kplot.Kini,KannealG(:,3)];
    Kplot.Alpha=[Kplot.Alpha,KannealG(:,4)];
    Kplot.MeanData=[Kplot.MeanData,MeanExp(:,1)];
    Kplot.MeanFit=[Kplot.MeanFit,MeanExp(:,2)];
    Kplot.MeanFitSteady=[Kplot.MeanFitSteady,MeanExp(:,3)];
end
save([OutFolderMain,['Results','.mat']],'KannealELs','MeanExpELs','EL','TimeLabelMean','Kplot')

Ytime=TimeLabelMean';
Xel=mean(EL,1)';

figure%2D
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    nexttile
    hold on
    plot(Ytime,Kplot.MeanData(:,ELi),'-o','LineWidth',2,'color',Color(1))
    plot(Ytime,Kplot.MeanFit(:,ELi),'-o','LineWidth',2,'color',Color(2))
    plot(Ytime,Kplot.MeanFitSteady(:,ELi),'-o','LineWidth',2,'color',Color(3))
%     plot(Ytime,Kplot.Kratio(:,ELi),'-s','LineWidth',2,'color','c')
    ylabel('MeanExpression')
%     yyaxis right
%     plot(Ytime,Kplot.OnRatio(:,ELi),'--d','LineWidth',2,'color','m')
%     ylabel('ON');ylim([0 1]);
    xlabel('Time/sec');title(ELlabel)
    legend('Data','Fit','FitInSteady')
end
saveas(gcf,[OutFolderMain,'Mean vs Time','.fig']);
saveas(gcf,[OutFolderMain,'Mean vs Time','.png']);

figure%2D
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    nexttile
    hold on
    plot(Ytime,Kplot.Kini(:,ELi),'-o','LineWidth',2,'color',Color(1))
    ylabel('Kini');ylim([0 1]);
    yyaxis right
    plot(Ytime,Kplot.Alpha(:,ELi),'--d','LineWidth',2,'color',Color(2))
    ylabel('Alpha');ylim([0 1]);
    xlabel('Time/sec');title(ELlabel)
    legend('Kini','Alpha')
end
saveas(gcf,[OutFolderMain,'Kini vs Time','.fig']);
saveas(gcf,[OutFolderMain,'Kini vs Time','.png']);


figure%2D
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    nexttile
    hold on
    plot(Ytime,Kplot.Kon(:,ELi),'-o','LineWidth',2,'color',Color(1))
    ylabel('K')
    plot(Ytime,Kplot.Koff(:,ELi),'--d','LineWidth',2,'color',Color(2))
    xlabel('Time/sec');title(ELlabel)
    ylim([0 0.05])
    legend('Kon','Koff')
end
saveas(gcf,[OutFolderMain,'K vs Time','.fig']);
saveas(gcf,[OutFolderMain,'K vs Time','.png']);


% figure%2D
% for ELi=1:size(EL,2)
%     el=EL(:,ELi);elre=ELre(:,ELi);
%     ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
%     nexttile
%     hold on
%     plot(Ytime,Kplot.ON(:,ELi),'-o','LineWidth',2,'color',Color(1))
%     ylabel('K')
%     plot(Ytime,Kplot.OFF(:,ELi),'--d','LineWidth',2,'color',Color(2))
%     xlabel('Time/sec');title(ELlabel)
%     ylim([0 0.2])
%     legend('ON','OFF')
% end
% saveas(gcf,[OutFolderMain,'NF vs Time','.fig']);
% saveas(gcf,[OutFolderMain,'NF vs Time','.png']);

figure%2D
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    nexttile
    hold on
    plot(Ytime,Kplot.Sf(:,ELi),'-o','LineWidth',2,'color',Color(1))
    plot(Ytime,Kplot.S0(:,ELi),'-o','LineWidth',2,'color',Color(2))
    plot(Ytime,Kplot.S1(:,ELi),'-o','LineWidth',2,'color',Color(3))
    xlabel('Time/sec');title(ELlabel)
    ylim([0 1])
    legend('Soff','S0','S1')
end
saveas(gcf,[OutFolderMain,'Pstate vs Time','.fig']);
saveas(gcf,[OutFolderMain,'Pstate vs Time','.png']);

if length(Xel)>1
    [X,Y] = meshgrid(Xel,Ytime);
    figure
    set(gcf,'Position',[100 100 2200 900])
    subplot(1,4,1)
    mesh(X,Y,Kplot.Kratio)
    xlabel('EL');ylabel('Time');zlabel('Kratio')
    subplot(1,4,2)
    mesh(X,Y,Kplot.Kini)
    xlabel('EL');ylabel('Time');zlabel('Kini')
    subplot(1,4,3)
    mesh(X,Y,Kplot.Alpha)
    xlabel('EL');ylabel('Time');zlabel('Alpha')
    subplot(1,4,4)
    mesh(X,Y,Kplot.MeanData)
    xlabel('EL');ylabel('Time');zlabel('Exp')
    sgtitle('K vs EL & Time')
    saveas(gcf,[OutFolderMain,'K vs EL & Time','.fig']);
    saveas(gcf,[OutFolderMain,'K vs EL & Time','.png']);
end