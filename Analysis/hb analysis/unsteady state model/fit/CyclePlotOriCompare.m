%%CyclePlotOriCompare
%% Ini
ModelPlot='Ton';
InheritSwitchPlot=1;
RepeatNum=1;
TimeBinPlot=0.5;
CombinNum=3;
PonLock=0;
ELLim=[0.2 0.4];
CompareNum='Vel Compare\0520Mingle[0.2 0.4]\';
% ForNum=["Vel25","Vel30","Vel32.5","Vel35","Vel37.5","Vel40","Vel42.5","Vel45","Vel50","Vel55","Vel60","Vel80"];
% ForNum=["Vel20","Vel25","Vel30","Vel35","Vel40","Vel45","Vel50","Vel55","Vel80"];
ForNum=["Vel30","Vel40","Vel50"];
LabelTipG=@(x) ...
    ['-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-',x,'-KL[0.01 0.12 1]-KU[0.2 0.12 1.3]-TimeStart60-KDisLost[10 10 1 1 1]-Artifical-Filter5-timeBinMove1-PonBoole0-re6-0520fix-Vel\'];
    % ['-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-',x,'-KL[0.01 0.1 0.60]-KU[0.2 0.1 0.8]-TimeStart60-KDisLost[10 10 1 1 1]-Artifical-Filter5-timeBinMove1-PonBoole0-re4\'];
Cycle=[11];
h8=figure('Position',[800 600 1000 400]);
h9=figure('Position',[800 600 1000 400]);
LldAbsG=[];
for iii=1:size(ForNum,2)
    LabelTip=LabelTipG(char(ForNum(iii)));
    titleNum=char(ForNum(iii));

kName=["k01" "k02" "k10" "k12" "k20" "k21" "ki0" "ki1" "ki2" "Ton" "dataintensity1" "alpha0"];
OutFolderMain=@(Cycle) ['Y:\TimeFit\OutputModel\',ModelPlot,'\Cycle',Cycle,...
    '\Inherit',num2str(InheritSwitchPlot),'-AlpOn-Turn',num2str(RepeatNum),'-TimeBin',num2str(TimeBinPlot),'-Combin',num2str(CombinNum),'-PonLock',num2str(PonLock),LabelTip];
OutputRoad=['Y:\TimeFit\OutputModel\',ModelPlot,'\CyclePlot\',...
    CompareNum];mkdir(OutputRoad);
%% Loading
Color=["#19CAAD","#F4606C","#9999CC"];

CycleTime=[0,7,9+7];
DivisionTime=[1,1,1.5];
Ytime=CycleTime+DivisionTime;



for CycleI=Cycle
    CycleLabel=num2str(CycleI);
    RoadCycle=OutFolderMain(CycleLabel);
    load([RoadCycle,'Results0.mat'])
    Kstart=[0.05,0.10,0.65,0,0];
    
    switch CycleI
        case 11
            CombinPart=[3 2 2 5];
            MeanLld=[1 2 3];
        case 12'
            CombinPart=[5 2 3];
        case 13
            CombinPart=[6 2 2 2 2 2 6];
    end
ELlegends={};
c_map = [cool(ceil(size(EL,2)/2));spring(floor(size(EL,2)/2))];
ELuse=mean(EL,1)>=ELLim(1)&mean(EL,1)<=ELLim(2);
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    
    ELlabels{ELi}=['EL ',num2str(el(1)),'-',num2str(el(2))];
    ELlegends=[ELlegends;['EL ',num2str(el(1)),'-',num2str(el(2))]];
    KannealG=KPlotannealELs{ELi,1};
    PiniG=PstatePlotMeanELs{ELi,1};
    LldG=LldPlotMeanELs{ELi,1};
    TimeLabelMean=TimePlotMeanELs{ELi,:}+Ytime(CycleI-10)*60;
% calculate the k lost
[C,ia,ic]=unique(KannealG(:,1),'stable');
KannealGUnique=[Kstart;KannealG(ia,:)];
KannealDiff=KannealGUnique(2:end,:)-KannealGUnique(1:end-1,:);
KdistanceG=[];
for ii=1:(size(KannealDiff,1)-1)
    Kdistance=1e-1*dist([10,10,1,1,1]'.*[KannealDiff(ii+1,1:5)',KannealDiff(ii,1:5)']);
    KdistanceG(ii)=Kdistance(1,2);
end
KdistancePart=[];
Starti=1;
for ii=1:size(CombinPart,2)
    Endi=Starti+CombinPart(ii);
    Endi=min([Endi;size(KdistanceG,2)+1]);
    KdistancePart(ii)=sum(KdistanceG(Starti:(Endi-1)),2);
    Starti=Endi;
end
%% Lost
try
    figure(h8)
    subplot(1,size(ForNum,2),iii)
        hold on
        plot(TimeLabelMean,abs(LldG(:,1)),'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
        ylabel('Lost state frequence')
        xlabel('Time/sec');title('Likelihood vs. Time')
        title(titleNum)
        ylim([0 15])

    % figure(h9)
    % [C,ia,ic]=unique(abs(LldG(:,1)),'stable');
    LldAbs=real(LldG(:,1))+abs(imag(LldG(:,1)));
    % for ii=1:max(ic)
    %     LldAbs(ic==ii)=LldAbs(ic==ii)-KdistancePart(ii);
    % end
    % subplot(1,size(ForNum,2),iii)
    %     hold on
    %     plot(TimeLabelMean,LldAbs,'-o','MarkerSize',3,'LineWidth',1.5,'color',c_map(ELi,:))
    %     ylabel('Lost state frequence')
    %     xlabel('Time/sec');title('Likelihood distrbution vs. Time')
    %     title(titleNum)
    %     ylim([0 15])
    % 
    % LldAbsU=LldAbs(ia)./CombinPart';
    % LldAbsG(ELi,iii)=mean(LldAbsU(MeanLld));
    LldAbsG(ELi,iii)=mean(LldAbs);
catch
    close gcf
    disp('Not Find: Lost information')
end

end
end
end

figure('Position',[200 400 1200 600])
X = categorical(ForNum);
X = reordercats(X,ForNum);
plot(X,mean(LldAbsG(ELuse,:),1)./mean(LldAbsG(ELuse,:),'all'),'Color','#A2142F','LineWidth',2)
title('Lost vs. Vel');
ylabel('Lost')
grid on
set(gca,'FontSize',15);
saveas(gcf,[OutputRoad,'lld vs Vel','.fig']);
saveas(gcf,[OutputRoad,'lld vs Vel','.png']);
%% 
figure(h8)
% legend(ELlegends,'Location','northeastoutside')
% pbaspect([3 1 1])
saveas(gcf,[OutputRoad,'lld vs Time','.fig']);
saveas(gcf,[OutputRoad,'lld vs Time','.png']);
figure(h9)
% legend(ELlegends,'Location','northeastoutside')
% pbaspect([3 1 1])
saveas(gcf,[OutputRoad,'lld distrbution vs Time','.fig']);
saveas(gcf,[OutputRoad,'lld distrbution vs Time','.png']);
%% END