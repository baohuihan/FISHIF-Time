function [KannealG,MeanExpFitG,PiniG,TimeLabelPlot,KannealGPlot,PstatePlot,LldPlot]=ModuleCombinFitTime(RnaSignal,TimeLabel,RnaLenth,kfitTime,OutFolder,InheritSwitch,RepeatNum,Bound,Model,CombinNum,CombinPart,PonLock,PonTime,PonLocation,Vel)
%%% A program to fit the experiment data by model

%% Data Initialization
TimeLabel=TimeLabel*60; % Convert time label from minutes to seconds
TimeNum=size(TimeLabel,2);
DataCell=RnaSignal'; % Transpose RnaSignal for processing
UnLockIndex=logical(Bound(1,:));
Ub=Bound(2,:);
Lb=Bound(3,:);
% Data preprocessing and visualization steps can be encapsulated in a function
% Example: [TimeLabel, DataCell] = preprocessData(TimeLabel, DataCell);

%% Initialize required parameters
KFitG=[]; MeanExpFitG=[];ExpFitG=[];PiniS=[];
% FitSwitch=1; % Flag to control fitting process
% RepeatNum=3; % Number of repetitions for the fitting process
% IniModel=1; % Model initialization method
InheritModel='Label';
MaxRnaNum=100; % Maximum number of RNA
Edge=-0.5:1:MaxRnaNum+0.5; % Binning edges for histogram
Nbin=0:1:MaxRnaNum; % Binning for histogram
Pini=[1 0 0]'; % Initial state probabilities
% RnaLenth=[497+2725+3+283 2277 407+509]; % Lengths of RNA components
Tstay=sum(RnaLenth)/Vel; % Average stay time, assuming 25bp/s processing rate
Label=RnaLenth/Vel; % Labeling for RNA processing
LabelSum=cumsum(Label); % Cumulative sum of labels
LabelPart=Label/Tstay; % Proportion of labels
LabelStrength=sum([0 0.5 1].*LabelPart); % Weighted sum of labels
KannealG=[]; PiniG=[]; LldG=[]; 
KannealTime=TimeLabel'-[Tstay kfitTime];
Hanneal=figure; % Initialize a figure for annealing plots
set(gcf,'Position',[100 100 2200 900]); % Set figure size
% hInherit=figure; % Initialize a figure for annealing plots
% set(gcf,'Position',[100 100 2400 900]); % Set figure size


% Encapsulating the annealing process in a function would be beneficial
% Example: [KannealG, PiniG, LldG] = performAnnealing(...);

%% Analyze Data at each time point
FitI=zeros(1,TimeNum);
for T_i = 1:TimeNum
    Turn=0;LldPilot=inf;
    if FitI(T_i)
        continue
    end
    DataTime=DataCell{1,T_i};
    %% Setting up boundaries and initializing parameters
    % Consider encapsulating this into a function for clarity
    [Nbar,Kstart] = setupAnnealingParameters(DataTime,Edge,LabelStrength,Tstay,Model);%%Model Fix
    if T_i>1
        Kstart(UnLockIndex)=KannealG(T_i-1,:);%%K start anneal
    end
    if ~isempty(CombinPart)
        CombinNum=CombinPart(1);
        CombinPart(1)=[];
    end
    TimeCombin=T_i:(T_i+CombinNum);
    TimeCombin(TimeCombin>TimeNum)=[];
    TimeCombinSize=size(TimeCombin,2);Kl=sum(UnLockIndex,2);
    KstartAnneal=repmat(Kstart(UnLockIndex),1,TimeCombinSize);%Select parameters to anneal

    LbAnneal=repmat(Lb(UnLockIndex),1,TimeCombinSize);
    UbAnneal=repmat(Ub(UnLockIndex),1,TimeCombinSize);
    if PonLock==1
        PonLocationRep=logical(repmat(PonLocation(UnLockIndex),1,TimeCombinSize));
        LbAnneal(PonLocationRep)=PonTime(TimeCombin);
        UbAnneal(PonLocationRep)=PonTime(TimeCombin);
        KstartAnneal(PonLocationRep)=PonTime(TimeCombin);
    end
    % LabelPart=LabelWait/sum(LabelWait);
    % [Kanneal,~,~,~] = simulannealbnd(@(x) ...
    %     CombinFitFSP(x,TimeCombin,DataCell,Edge,LabelStrength,Tstay,Model, TimeLabel, PiniG, KannealG, KannealTime, InheritSwitch,InheritModel, Label, UnLockIndex),...
    %     KstartAnneal,LbAnneal,UbAnneal);
    % [lld,PDtbG,PiniG0,KannealTime]=CombinFitFSP(Kanneal,TimeCombin,DataCell,Edge,LabelStrength,Tstay,Model, TimeLabel, PiniG, KannealG, KannealTime, InheritSwitch,InheritModel, Label, UnLockIndex);
    KannealGTurn=[];
    while Turn < RepeatNum % repeat turn
        Turn = Turn + 1;
        [Kanneal,~,~,~] = simulannealbnd(@(x) ...
            CombinFitFSP(x,TimeCombin,DataCell,Edge,LabelStrength,Tstay,Model, TimeLabel, PiniG, KannealG, KannealTime, InheritSwitch,InheritModel, Label, UnLockIndex),...
            KstartAnneal,LbAnneal,UbAnneal);
        %%% Lld&Pini
        [lld,PDtbG,PiniG0,KannealTime]=CombinFitFSP(Kanneal,TimeCombin,DataCell,Edge,LabelStrength,Tstay,Model, TimeLabel, PiniG, KannealG, KannealTime, InheritSwitch,InheritModel, Label, UnLockIndex);
        if lld<=LldPilot
            LldPilot=lld;
            KannealOut=Kanneal;
        end
    end
    [lld,PDtbG,PiniG0,KannealTime]=CombinFitFSP(KannealOut,TimeCombin,DataCell,Edge,LabelStrength,Tstay,Model, TimeLabel, PiniG, KannealG, KannealTime, InheritSwitch,InheritModel, Label, UnLockIndex);
    %% Analysis and fitting process
    for Ci=1:TimeCombinSize
        Tci=TimeCombin(Ci);
        KannealG(Tci,:)=Kanneal(1,(1:Kl)+(Ci-1)*Kl);
        PiniG(Tci,:)=PiniG0(Ci,:);
        LldG(Tci,:)=lld;
        FitI(Tci)=1;
        Pdtb=PDtbG(Ci,:);
        [Nbar,~] = setupAnnealingParameters(DataCell{1,Tci},Edge,LabelStrength,Tstay,Model);%%Model Fix
        visualizeDataandfit(Hanneal, Nbin, Nbar, [Pdtb], TimeLabel, Tci);
    end
    % [~,Pdtb,~,Pinis,~] =  FSP_FreeSet(KannealOut,DataTime,PiniStart,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit);%%Model Fix
    % [~,PdtbInherit,~,~,~] =  FSP_FreeSet(KannealOut,DataTime,PiniStart,LabelPart,0,Kstart,UnLockIndex,0,DtdInherit);%%Model Fix
    % [~,Pdtb_steady,~,~,~] =  FSP_FreeSet_Steady(KannealOut,DataTime,PiniStart,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit);
    % MeanExpFitG(T_i,:)=[sum(Pdtb(1:71).*[0:1:70],2) sum(Pdtb_steady(1:71).*[0:1:70],2)];%[Mean MeanInSteady]
    % ExpFitG(T_i,:)=[sum(Pdtb(2:end),2) sum(Pdtb_steady(2:end),2)];%[Mean MeanInSteady]Logic
    %% Visualization and result storage
    % Encapsulate visualization into a separate function
    

    % ... rest of the code for the loop ...
end

%% Saving and plotting results
% These steps can also be encapsulated for better readability and reuse
% Example: saveResults(OutFolder, KannealG, PiniG, LldG, ...);
figure(Hanneal)
saveas(gcf,[OutFolder,'FitAnneal','.fig']);
saveas(gcf,[OutFolder,'FitAnneal','.png']);
% figure(hInherit)
% saveas(gcf,[OutFolder,'TimeInherit','.fig']);
% saveas(gcf,[OutFolder,'TimeInherit','.png']);
% save([OutFolder,['Results','.mat']],'KannealG','PiniG','PiniS','LldG','Lb','Ub','UnLockIndex')
%% k-transform fig
switch InheritModel
    case 'Near'
        TimeLabelPlot=TimeLabel;KannealGPlot=KannealG;
    case 'Label'
        [timeBoundaries, featureValues] = createTimelines(KannealTime, KannealG);
        [~, PstatePlot] = createTimelines(KannealTime, PiniG);
        [~, LldPlot] = createTimelines(KannealTime, LldG);
        TimeLabelPlot=mean([timeBoundaries(1:end-1);timeBoundaries(2:end)],1);
        KannealGPlot=featureValues;

        KUndefinedIndex=find(KannealGPlot(:,1)==0);
        if ~isempty(KUndefinedIndex)
            for iii=KUndefinedIndex'
                try
                    KannealGPlot(iii,:)=mean(KannealGPlot([iii-1 iii+1],:),1);
                catch
                    KannealGPlot(iii,:)=mean(KannealGPlot(iii-1,:),1);
                end
            end
        end
end
save([OutFolder,['Results','.mat']],'KannealG','KannealTime','PiniG','PiniS','LldG','Lb','Ub','UnLockIndex','TimeLabelPlot','KannealGPlot','PstatePlot')
figure%Plot Kon\Koff\Kini vs. Time
title('Anneal K')
yyaxis left
hold on
plot(TimeLabelPlot,KannealGPlot(:,1),'-.c','LineWidth',2)
plot(TimeLabelPlot,KannealGPlot(:,2),'-.m','LineWidth',2)
ylim([0 0.01]);ylabel('K');xlabel('Evolution Time');xlim([0 max(TimeLabelPlot)-1])
yyaxis right
plot(TimeLabelPlot,KannealGPlot(:,3),'LineWidth',2)
ylim([0 0.5]);ylabel('Kini');legend('Kon','Koff','Kini')
saveas(gcf,[OutFolder,'KAnneal','.fig']);
saveas(gcf,[OutFolder,'KAnneal','.png']);

figure%Plot Kon/Koff vs. Time
title('Anneal Kon/Koff')
plot(TimeLabelPlot,KannealGPlot(:,1)./(KannealGPlot(:,2)+KannealGPlot(:,1)),'-.c','LineWidth',2)
ylim([0 2]);ylabel('Kon/Koff');xlabel('Evolution Time');xlim([0 max(TimeLabelPlot)-1])
saveas(gcf,[OutFolder,'KratioAnneal','.fig']);
saveas(gcf,[OutFolder,'KratioAnneal','.png']);
end
