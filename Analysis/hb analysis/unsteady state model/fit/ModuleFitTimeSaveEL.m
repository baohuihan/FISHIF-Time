function ModuleFitTimeSaveEL(TimeMeanELs,EL,KannealELs,PiniELs,MeanExpELs,Label,OutFolderMain)
Hk1=figure;Hk2=figure;Hk3=figure;Hk4=figure;Hk5=figure;
Hp1=figure;Hp2=figure;Hp3=figure;
YlabelK=["Kon","Koff","Kini","ActiveR","Double"];
YlabelP=["S0","S1","S2"];
ELlabels={};
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    ELlabels{ELi}=['EL ',num2str(el(1)),'-',num2str(el(2))];
    KannealG=KannealELs{ELi,1};
    PiniG=PiniELs{ELi,1};
    MeanExp=MeanExpELs{ELi,:};
    TimeLabelMean=TimeMeanELs{ELi,:};
    
    for ii=1:size(KannealG,2)
        eval(['figure(Hk',num2str(ii),')'])
        hold on
        plot(TimeLabelMean,KannealG(:,ii),'-o','LineWidth',2)
        Kname=char(YlabelK(ii));
        ylabel(Kname)
        xlabel('Time/min');
        title([Kname,' vs. time with different ELs'])
    end

    % for ii=1:size(PiniG,2)
    %     eval(['figure(Hp',num2str(ii),')'])
    %     hold on
    %     plot(TimeLabelMean,PiniG(:,ii),'-o','LineWidth',2)
    %     Pname=char(YlabelP(ii));
    %     ylabel(Pname)
    %     xlabel('Time/min');
    %     title([Pname,' vs. time with different ELs'])
    % end
end
for ii=1:size(YlabelK,2)
    eval(['figure(Hk',num2str(ii),')'])
    % legend(ELlabels)
    saveas(gcf,[OutFolderMain,char(YlabelK(ii)),' vs Time ',Label,' ELs.fig']);
    saveas(gcf,[OutFolderMain,char(YlabelK(ii)),' vs Time ',Label,' ELs.png']);
end
for ii=1:size(YlabelP,2)
    eval(['figure(Hp',num2str(ii),')'])
    legend(ELlabels)
    saveas(gcf,[OutFolderMain,char(YlabelP(ii)),' vs Time ',Label,' ELs.fig']);
    saveas(gcf,[OutFolderMain,char(YlabelP(ii)),' vs Time ',Label,' ELs.png']);
end