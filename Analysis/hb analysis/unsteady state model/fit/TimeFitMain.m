function TimeFitMain(Cycle,Model,TimeBin,InheritSwitch,CombinNum)
%%% Main File of time fit process
%%% Main frame: WorkPath DataExtract TimeFit VisualOutput
%% Ini set
% load('Y:\TimeFit\OutputModel\Ton\CyclePlot\Use0417-ELcat\FittingResults.mat')
PonLock=0;
PonBoole=0;
Vel=40;
% Cycle='12';
% Model='TwoState';
% TimeBin=0.5;
% InheritSwitch=0;
LabelTip='-EL[0 0.6]bin0.1Move0.05-TimeWindowMingle-Vel40-KL[0.01 0.12 1]-KU[0.1 0.12 1.4]-KDisLost[10 5 0 1.2 0]FirstON-Artifical-Filter6-timeBinMove1-PonBoole0-Star0731-3\';
RepeatNum=1;
% EL=[[0:0.1:0.4];[0.2:0.1:0.6]];
EL=[[0:0.05:0.5];[0.1:0.05:0.6]];
% EL=[[0:0.05:0.4];[0.2:0.05:0.6]];
% EL=[[0.2:0.1:0.2];[0.4:0.1:0.4]];
switch Model
    case 'Ton'
        Bound=[...
            [1.00 0 1.00 0 0 0 0 1.0 0 1 0 1];...
            [0.10,0,0.12,0,0,0,0,1.4,0,1,1,1];...
            [0.01,0,0.12,0,0,0,0,1.0,0,0,1,0]];
    case 'TwoState'
        Bound=[...
            [1.000 0 1.000 0 0 0 0 1.00 0 0 0 1];... 
            [1.000,0,1.000,0,0,0,0,0.55,0,0,1,1];...
            [0.001,0,0.001,0,0,0,0,0.01,0,0,1,0]];
    case 'Ton&ThreeState'
        Bound=[...
            [1.00 0 1.00 1.00 0 1.00 0 1.00 1 1 0 1];...
            [0.20,0,0.20,0.20,0,0.20,0,0.90,1.6,1,1,0];...
            [0.01,0,0.01,0.01,0,0.10,0,0.60,1.3,0,1,0]];
end
PonLocation=[0 0 0 0 0 0 0 0 0 1 0 0];

% RnaLenth=[497+2725+3+283 2277 407+509];
RnaLenth=[3+283 2277 407+509];
kfitTime=(3+283+2277/2)/Vel;
%k=[k01 k02 k10 k12 k20 k21 ki0 ki1 ki2 Ton dataintensity1 alpha0]
%% IniDeal
switch Cycle
    case '11'
        DataDeal=[1 1 1];
        TimeBins=[[0:TimeBin:7-TimeBin];[TimeBin:TimeBin:7]];
        % artificalFilter=[6 24 26];
        artificalFilter=[21 23];
        artificalTimeCombin=[];
        % artificalFilter=[6 23 25 27];
        artificalFilterHold=[];
        CombinPart=[5 5];
        StdDel=25;
        % CombinPart=[3 3];
        % FucPon=Fuc.cycle11;
        DivisionTime=1*60;
    case '12'
        DataDeal=[1 1 1];
        TimeBins=[[0:TimeBin:9-TimeBin];[TimeBin:TimeBin:9]];
        % TimeBins=[[0:TimeBin/2:9-TimeBin];[TimeBin:TimeBin/2:9]];
        artificalFilter=[2 6 34 35];
        artificalFilterHold=[1 36];
        artificalTimeCombin=[];
        % CombinPart=[5 2 3];
        CombinPart=[5 2 3];
        StdDel=25;
        % CombinPart=[6 8];
        % CombinPart=[3 4];
        % FucPon=Fuc.cycle12;
        DivisionTime=1*60;
    case '13'
        DataDeal=[1 1 1];
        TimeBins=[[0:TimeBin:14-TimeBin];[TimeBin:TimeBin:14]];
        % artificalFilter=[4 7 56];
        artificalFilter=[55];
        % artificalFilter=[];
        artificalFilterHold=[];
        artificalTimeCombin=[];
        CombinPart=[6 4 4 2 6];
        StdDel=25;
        % CombinPart=[5 3 2 2 2 2 4];
        % CombinPart=[4 3 2 4];
        % FucPon=Fuc.cycle13;
        DivisionTime=1*60;
end
% PonT=FucPon(0:0.01:1000);
% clear FucPon
% addAttachedFiles(gcp, 'FuncPon.m')

ELre=1-EL;ELre=[ELre(2,:);ELre(1,:)];
OutFolderMain=['Y:\TimeFit\OutputModel\',Model,'\Cycle',Cycle,...
    '\Inherit',num2str(InheritSwitch),'-AlpOn-Turn',num2str(RepeatNum),'-TimeBin',num2str(TimeBin),'-Combin',num2str(CombinNum),'-PonLock',num2str(PonLock),LabelTip];
mkdir(OutFolderMain);
Color=["#19CAAD","#F4606C","#9999CC"];
%% Embryo Path
TeamPath=['Y:\TimeFit\number_hb_',Cycle,'.xlsx'];%Excel of work list.
% Sheets={'Heng','wjy1' 'wjy2' 'bhh'}; 
Sheets={'all'};
% RnaOrder=[1 2 3 1];
% SheetPaths={'F:\','X:\bhh-fish\oreR_wjy\','X:\bhh-fish\oreR_wjy\OreR_hb1hb7CDS_60X\','X:\bhh-fish\oreR_wjy\'};
TimeLabel=[];EmNumber=0;RnaChannels=[];
for Si=1:size(Sheets,2)
    Read=readcell(TeamPath,'sheet',Sheets{Si});
    Nomissing=cell2mat(cellfun(@(x) mean(~ismissing(x)),Read(:,1),'UniformOutput',false))&...
        cell2mat(cellfun(@(x) mean(~ismissing(x)),Read(:,strcmpi('noverT',Read(1,:))),'UniformOutput',false));
    TimeLabel=cat(1,TimeLabel,cell2mat(Read(Nomissing,strcmpi('noverT',Read(1,:)))));%Time(cycle) of embryo /min
    RnaChannels=cat(1,RnaChannels,cell2mat(Read(Nomissing,strcmpi('Rna Channel',Read(1,:)))));
    ReadWork=Read(Nomissing,:);
    for Ei=1:size(ReadWork,1)
        EmPath{EmNumber+Ei,1}=[num2str(ReadWork{Ei,1}),num2str(ReadWork{Ei,2}),'\Results\',ReadWork{Ei,3},'_new.mat'];%Embryo results path
        EmName{EmNumber+Ei,1}=ReadWork{Ei,2};%embryo name
        RnaSelect(EmNumber+Ei,1)=RnaChannels(Ei,1);%Rna choose
    end
    EmNumber=EmNumber+size(ReadWork,1);
end
disp(['Working embryo number: ',num2str(EmNumber)])
%% For of EL
TimeLabelOri=TimeLabel;KannealELs={};MeanExpELs={};PiniELs={};TimeMeanELs={};KPlotannealELs={};TimePlotMeanELs={};PstatePlotMeanELs={};LldPlotMeanELs={};
PoolCreate(size(EL,2))
parfor ELi=1:size(EL,2)
    
    %% Data Extract
    el=EL(:,ELi);elre=ELre(:,ELi);
    ELlabel=['EL ',num2str(el(1)),'-',num2str(el(2))];
    OutFolder=[OutFolderMain,ELlabel,'\'];
    mkdir(OutFolder);
    [TimeLabelSort,SortIndex]=sort(TimeLabelOri);
    EmPathSort=EmPath(SortIndex);
    EmNameSort=EmName(SortIndex);
    RnaSelectSort=RnaSelect(SortIndex);%Sort by time
    [RnaSignal,ErrorIndex]=RnaSignalExtract(EmPathSort,RnaSelectSort,el);%Extract
    [RnaSignalReverse,~]=RnaSignalExtract(EmPathSort,RnaSelectSort,elre);%Extract
    [RnaSignalNoReverseT,~]=RnaSignalExtract(EmPathSort,RnaSelectSort,[0.2 0.4]');%Deal reverse el
    [RnaSignalReverseT,~]=RnaSignalExtract(EmPathSort,RnaSelectSort,[0.6 0.8]');%Deal reverse el
    
    RnaSignalUse={};
    for ii=1:size(RnaSignal,1)
        % RnaSignalUse(ii,1)=RnaSignal(ii);
        [~,I]=max([mean(RnaSignalNoReverseT{ii}),mean(RnaSignalReverseT{ii})]);
        if I==1
            RnaSignalUse{ii,1}=RnaSignal{ii}*DataDeal(RnaSelectSort(ii));
        else
            RnaSignalUse{ii,1}=RnaSignalReverse{ii};
        end
    end
    RnaSignalInput=RnaSignalUse(~ErrorIndex,:);
    TimeLabelInput=TimeLabelSort(~ErrorIndex,:);
    RnaSelectSortInput=RnaSelectSort(~ErrorIndex,:);
    EmNameInput=EmNameSort(~ErrorIndex,:);%delet missing data
    %% show data%%%Time Fit
    ZeroSignalIndex=cellfun(@isempty,RnaSignalInput);
    RnaSignalInput(ZeroSignalIndex)={[0;0]};

    MPmodel='';
    MeanExp=MeanPlot(TimeLabelInput*60,RnaSignalInput,MPmodel);%per process
    xlabel('Time/sec');ylabel('MeanExpr')
    hold on
    for ii=unique(RnaSelectSortInput)'
        MeanPlotScatter(TimeLabelInput(RnaSelectSortInput==ii)*60,RnaSignalInput(RnaSelectSortInput==ii),MPmodel);
    end
    legend('all','1','2','3')
    
    %% Data filtering
    
    TFrm=zeros(size(RnaSignalInput,1),1);
    artificalTimeCombinE=TFrm;
    artificalTimeCombinE(artificalTimeCombin)=1;
    TFrm(MeanExp>40)=1;%% Delet outliers
    scatter(TimeLabelInput(TFrm==1)*60,MeanExp(TFrm==1),500,'x','LineWidth',3)
    RnaSignalStd=cell2mat(cellfun(@(x) std(x(x~=0)),RnaSignalInput,'UniformOutput',false));
    % figure
    % for ii=1:length(RnaSignalInput)
    %     nexttile;
    %     histogram(RnaSignalInput{ii},0:1:70,'Normalization','probability')
    %     ylim([0 0.05]);xlim([-0.5 70])
    %     title(num2str(RnaSignalStd(ii)))
    % end
    TFrm(RnaSignalStd>StdDel)=2;
    
    scatter(TimeLabelInput(TFrm==2)*60,MeanExp(TFrm==2),500,'x','LineWidth',3)
    TFrm(artificalFilter)=3;
    scatter(TimeLabelInput(TFrm==3)*60,MeanExp(TFrm==3),500,'x','LineWidth',3)
    TFrm=logical(TFrm);
    xlim([0 max(TimeLabelInput*60)+60]);title(ELlabel)
    TFrm(artificalFilterHold)=0;
    saveas(gcf,[OutFolder,'EL',num2str(ELi),' Mean-Time','.png']);
    saveas(gcf,[OutFolder,'EL',num2str(ELi),' Mean-Time','.fig']);
    DataCell=RnaSignalInput(~TFrm');
    TimeLabel=TimeLabelInput(~TFrm)';
    artificalTimeCombinE(TFrm)=[];artificalTimeCombinE=logical(artificalTimeCombinE');
    %% Data Combin
    
    TimeLabelMean=[];
    DataCellMean={};
    %combine data in 2 min
    TimeIndex=1:size(TimeLabel,2);
    for TBi=1:size(TimeBins,2)
        TimeInBinCheck=TimeLabel>=TimeBins(1,TBi)&TimeLabel<TimeBins(2,TBi);
        if ~isempty(artificalTimeCombinE)
            artificalTimeCombinF=TimeIndex(artificalTimeCombinE');
            if sum(ismember(TimeIndex(TimeInBinCheck),artificalTimeCombinF),2)
                TimeLabelMean=[TimeLabelMean,mean(TimeLabel(artificalTimeCombinE),2)];
                DataCellMean=cat(1,DataCellMean,cat(1,DataCell{artificalTimeCombinE'}));
                artificalTimeCombinE=[];
            end
        end
        TimeInBin=TimeLabel>=TimeBins(1,TBi)&TimeLabel<TimeBins(2,TBi)&~ismember(TimeIndex,artificalTimeCombinE);
        if sum(TimeInBin,2)==0
            continue
        end
        if sum(TimeInBin,2)<0
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
    [TimeLabelMean,I]=sort(TimeLabelMean);
    DataCellMean=DataCellMean(I);

    MeanExp=MeanPlot(TimeLabelMean*60,DataCellMean,'');%data combin
    xlabel('Time/sec');ylabel('MeanExpr');xlim([0 max(TimeLabelMean*60)+60]);title(ELlabel)
    saveas(gcf,[OutFolder,'Mean-Time-combine','.fig']);
    saveas(gcf,[OutFolder,'Mean-Time-combine','.png']);
    % [KannealG,MeanExpFitG,PiniG]=ModuleFitTime_FFN_Unlock20(DataCellMean,TimeLabelMean,OutFolder,InheritSwitch);
    if strcmp(Cycle,'13')
        TimeLabelMean=[1 TimeLabelMean];
        DataCellMean=cat(1,zeros(10,1),DataCellMean);
    end
    %% ModuleFitTime
    if PonLock==1
        Tstay=sum(RnaLenth)/Vel;
        KannealTime=TimeLabelMean'*60-mean([Tstay kfitTime],2)+DivisionTime;
        PonTime=FuncPon(Cycle,KannealTime);
        if PonBoole==1
            PonTime=PonTime./max(PonTime,[],"all");
        end
        % PonTime=PonT(round((TimeLabelMean*60+DivisionTime)/0.01));
    else
        PonTime=[];
    end
    if CombinNum==0
        [KannealG,MeanExpFitG,PiniG,TimeLabelPlot,KannealGPlot,PstatePlot]=...
            ModuleFitTime(DataCellMean,TimeLabelMean,RnaLenth,kfitTime,OutFolder,InheritSwitch,RepeatNum,Bound,Model,PonLock,PonTime);
    else
        [KannealG,MeanExpFitG,PiniG,TimeLabelPlot,KannealGPlot,PstatePlot,LldPlot]=...
            ModuleCombinFitTime(DataCellMean,TimeLabelMean,RnaLenth,kfitTime,OutFolder,InheritSwitch,RepeatNum,Bound,Model,CombinNum,CombinPart,PonLock,PonTime,PonLocation,Vel);
    end
    % load([OutFolder,['Results','.mat']]);MeanExpFitG=[];
    % ModuleFitTimeShow(DataCellMean,TimeLabelMean,DataCell,TimeLabel,RnaLenth,kfitTime,KannealG,PiniG,OutFolder,InheritSwitch,RepeatNum,Bound,Model)
    KPlotannealELs{ELi,1}=KannealGPlot;
    TimePlotMeanELs{ELi,1}=TimeLabelPlot;
    PstatePlotMeanELs{ELi,1}=PstatePlot;
    LldPlotMeanELs{ELi,1}=LldPlot;
    KannealELs{ELi,1}=KannealG;
    PiniELs{ELi,1}=PiniG;
    MeanExpELs{ELi,1}=[MeanExp,MeanExpFitG];%[MeanData MeanFit MeanFitSteady]
    TimeMeanELs{ELi,1}=TimeLabelMean;
    
end
%%Plot EL&Time
save([OutFolderMain,['Results0','.mat']],'KannealELs','MeanExpELs','EL','TimeMeanELs','PiniELs','KPlotannealELs','TimePlotMeanELs','PstatePlotMeanELs','LldPlotMeanELs','Bound')
ModuleFitTimeSaveEL(TimeMeanELs,EL,KannealELs,PiniELs,MeanExpELs,'Kdivi',OutFolderMain)
ModuleFitTimeSaveEL(TimePlotMeanELs,EL,KPlotannealELs,PstatePlotMeanELs,MeanExpELs,'Kmean',OutFolderMain)
ModuleFitTimeSave(TimeMeanELs,EL,KannealELs,PiniELs,MeanExpELs,OutFolderMain,Model,Color)