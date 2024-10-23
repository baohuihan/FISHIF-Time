function TimeFitMainCycle12
%%% Main File of time fit process
%%% Main frame: WorkPath DataExtract TimeFit VisualOutput
%% Ini set
EL=[[0:0.1:0.4];[0.2:0.1:0.6]];
ELre=1-EL;ELre=[ELre(2,:);ELre(1,:)];
OutFolderMain='Y:\TimeFit\Output\CDS\Cycle12\NoInherit-Two-AlpON-Ton-Turn1-TimeBin2-IniModel1-EL[0 0.6]bin0.2-Lock[Kini0.25 0.55][Koff 0.05][Kon 0.05]\';
mkdir(OutFolderMain);
Color=["#19CAAD","#F4606C","#9999CC"];
TimeHead='0.83_M';
%% Embryo Path
TeamPath='X:\bhh-fish\oreR_wjy\number_hb_12new2.xlsx';%Excel of work list.
Sheets={'Heng','wjy' 'wjy2' 'bhh'}; 
RnaOrder=[1 2 3 1];
SheetPaths={'F:\','X:\bhh-fish\oreR_wjy\','X:\bhh-fish\oreR_wjy\OreR_hb1hb7CDS_60X\','X:\bhh-fish\oreR_wjy\'};
TimeLabel=[];EmNumber=0;
for Si=[1 2 4]
    Read=readcell(TeamPath,'sheet',Sheets{Si});
    Nomissing=cell2mat(cellfun(@(x) mean(~ismissing(x)),Read(:,1),'UniformOutput',false))&...
        cell2mat(cellfun(@(x) mean(~ismissing(x)),Read(:,strcmpi(TimeHead,Read(1,:))),'UniformOutput',false));
    TimeLabel=cat(1,TimeLabel,cell2mat(Read(Nomissing,strcmpi(TimeHead,Read(1,:)))));%Time(cycle) of embryo /min
    ReadWork=Read(Nomissing,:);
    for Ei=1:size(ReadWork,1)
        EmPath{EmNumber+Ei,1}=[SheetPaths{Si},num2str(ReadWork{Ei,1}),'\Results\',ReadWork{Ei,2},'_new.mat'];%Embryo results path
        EmName{EmNumber+Ei,1}=ReadWork{Ei,2};%embryo name
        RnaSelect(EmNumber+Ei,1)=RnaOrder(Si);%Rna choose
    end
    EmNumber=EmNumber+size(ReadWork,1);
end
disp(['Working embryo number: ',num2str(EmNumber)])
%% For of EL
TimeLabelOri=TimeLabel;
KannealELs={};MeanExpELs={};
for ELi=1:size(EL,2)
    el=EL(:,ELi);elre=ELre(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    OutFolder=[OutFolderMain,ELlabel,'\'];
    mkdir(OutFolder);
%% Data Extract
[TimeLabelSort,SortIndex]=sort(TimeLabelOri);
EmPathSort=EmPath(SortIndex);
EmNameSort=EmName(SortIndex);
RnaSelectSort=RnaSelect(SortIndex);%Sort by time
[RnaSignal,ErrorIndex]=RnaSignalExtract(EmPathSort,RnaSelectSort,el);%Extract
[RnaSignalReverse,~]=RnaSignalExtract(EmPathSort,RnaSelectSort,elre);%Deal reverse el
for ii=1:size(RnaSignal,1)
    [~,I]=max([mean(RnaSignal{ii}),mean(RnaSignalReverse{ii})]);
    if I==1,RnaSignalUse(ii,1)=RnaSignal(ii);
    else,RnaSignalUse(ii,1)=RnaSignalReverse(ii);
    end
end
RnaSignalInput=RnaSignalUse(~ErrorIndex,:);
TimeLabelInput=TimeLabelSort(~ErrorIndex,:);
EmNameInput=EmNameSort(~ErrorIndex,:);%delet missing data
%% Time Fit
MeanExp=MeanPlot(TimeLabelInput*60,RnaSignalInput,'');%per process
xlabel('Time/sec');ylabel('MeanExpr')
hold on
TFrm=zeros(size(RnaSignal,1),1);
% TFrm([22 24 27 29 30 31 33 36])=1;
TFrm=logical(TFrm);%% Delet outliers
DataCell=RnaSignalInput;
TimeLabel=TimeLabelInput';
if sum(TFrm)~=0
    scatter(TimeLabelInput(TFrm)*60,MeanExp(TFrm),500,'x','LineWidth',3)
    DataCell=RnaSignalInput(~TFrm');
    TimeLabel=TimeLabelInput(~TFrm)';
end
xlim([0 max(TimeLabelInput*60)+60]);title(ELlabel)
saveas(gcf,[OutFolder,'Mean-Time','.png']);

TimeLabelMean=[];
DataCellMean={};
TimeBin=[[0:1:7];[1:1:8]];%combine data in 2 min
TimeIndex=1:size(TimeLabel,2);
for TBi=1:size(TimeBin,2)
    TimeInBin=TimeLabel>=TimeBin(1,TBi)&TimeLabel<TimeBin(2,TBi);
    if sum(TimeInBin,2)==0
        continue
    end
    if sum(TimeInBin,2)<5
        IndexU=TimeIndex(TimeInBin);
        for ii=1:size(IndexU,2)
            TimeLabelMean=[TimeLabelMean,TimeLabel(IndexU(ii))];
            DataCellMean=cat(1,DataCellMean,DataCell{IndexU(ii)});
        end
    else
        TimeLabelMean=[TimeLabelMean,mean(TimeLabel(TimeInBin),2)];
        DataCellMean=cat(1,DataCellMean,cat(1,DataCell{TimeInBin'}));
    end
end
MeanExp=MeanPlot(TimeLabelMean*60,DataCellMean,'');%data combin
xlabel('Time/sec');ylabel('MeanExpr');xlim([0 max(TimeLabelMean*60)+60]);title(ELlabel)
saveas(gcf,[OutFolder,'Mean-Time-combine','.png']);
[KannealG,MeanExpFitG]=ModuleFitTime_Ton(DataCellMean,TimeLabelMean,OutFolder);
KannealELs{ELi,1}=KannealG;
MeanExpELs{ELi,1}=[MeanExp,MeanExpFitG];%[MeanData MeanFit MeanFitSteady]
end
%%Plot EL&Time
Kplot=struct('Kon',[],'Koff',[],'Ton',[],'Kratio',[],'Kini',[],'Alpha',[],'MeanData',[],'MeanFit',[],'MeanFitSteady',[],'OnRatio',[]);
for ELi=1:size(EL,2)
    KannealG=KannealELs{ELi,1};
    MeanExp=MeanExpELs{ELi,:};
    Kplot.Kon=[Kplot.Kon,KannealG(:,1)];
    Kplot.Koff=[Kplot.Koff,KannealG(:,2)];
    Kplot.Ton=[Kplot.Ton,KannealG(:,4)];
    Kplot.Kratio=[Kplot.Kratio,KannealG(:,1)./KannealG(:,2)];
    Kplot.OnRatio=[Kplot.OnRatio,KannealG(:,1)./(KannealG(:,1)+KannealG(:,2))];
    Kplot.Kini=[Kplot.Kini,KannealG(:,3)];
    Kplot.Alpha=[Kplot.Alpha,KannealG(:,5)];
    Kplot.MeanData=[Kplot.MeanData,MeanExp(:,1)];
    Kplot.MeanFit=[Kplot.MeanFit,MeanExp(:,2)];
    Kplot.MeanFitSteady=[Kplot.MeanFitSteady,MeanExp(:,3)];
end
save([OutFolderMain,['Results','.mat']],'KannealELs','MeanExpELs','EL','TimeLabelMean','Kplot')

Ytime=TimeLabelMean';
Xel=mean(EL,1)';

figure%2D
for ELi=1:size(EL,2)
    el=EL(:,ELi);elre=ELre(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    nexttile
    hold on
    plot(Ytime,Kplot.MeanData(:,ELi),'-o','LineWidth',2,'color',Color(1))
    plot(Ytime,Kplot.MeanFit(:,ELi),'-o','LineWidth',2,'color',Color(2))
    plot(Ytime,Kplot.MeanFitSteady(:,ELi),'-o','LineWidth',2,'color',Color(3))
    plot(Ytime,Kplot.Kratio(:,ELi),'-s','LineWidth',2,'color','c')
    ylabel('MeanExpression')
    yyaxis right
    plot(Ytime,Kplot.OnRatio(:,ELi),'--d','LineWidth',2,'color','m')
    ylabel('ON');ylim([0 1]);
    xlabel('Time/sec');title(ELlabel)
    legend('Data','Fit','FitInSteady','Kratio: Kon/Koff','ON')
end
saveas(gcf,[OutFolderMain,'Mean vs Time','.fig']);
saveas(gcf,[OutFolderMain,'Mean vs Time','.png']);

figure%2D
for ELi=1:size(EL,2)
    el=EL(:,ELi);elre=ELre(:,ELi);
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
    el=EL(:,ELi);elre=ELre(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    nexttile
    hold on
    plot(Ytime,Kplot.Kon(:,ELi),'-o','LineWidth',2,'color',Color(1))
    ylabel('K')
    plot(Ytime,Kplot.Koff(:,ELi),'--d','LineWidth',2,'color',Color(2))
    xlabel('Time/sec');title(ELlabel)
    ylim([0 0.1])
    legend('Kon','Koff')
end
saveas(gcf,[OutFolderMain,'K vs Time','.fig']);
saveas(gcf,[OutFolderMain,'K vs Time','.png']);

figure%2D
for ELi=1:size(EL,2)
    el=EL(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    nexttile
    hold on
    plot(Ytime,Kplot.Ton(:,ELi),'-o','LineWidth',2,'color',Color(1))
    xlabel('Time/sec');title(ELlabel)
    ylim([0 1])
    legend('Ton')
end
saveas(gcf,[OutFolderMain,'Ton vs Time','.fig']);
saveas(gcf,[OutFolderMain,'Ton vs Time','.png']);

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