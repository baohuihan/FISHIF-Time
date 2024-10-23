function [KannealG,MeanExpFitG,PiniG,TimeLabelPlot,KannealGPlot,PstatePlot]=ModuleFitTime(RnaSignal,TimeLabel,RnaLenth,kfitTime,OutFolder,InheritSwitch,RepeatNum,Bound,Model,PonLock,PonTime)
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
Tstay=sum(RnaLenth)/25; % Average stay time, assuming 25bp/s processing rate
Label=RnaLenth/25; % Labeling for RNA processing
LabelSum=cumsum(Label); % Cumulative sum of labels
LabelPart=Label/Tstay; % Proportion of labels
LabelStrength=sum([0 0.5 1].*LabelPart); % Weighted sum of labels
KannealG=[]; PiniG=[]; LldG=[]; 
KannealTime=TimeLabel'-[Tstay kfitTime];
Hanneal=figure; % Initialize a figure for annealing plots
set(gcf,'Position',[100 100 2200 900]); % Set figure size
hInherit=figure; % Initialize a figure for annealing plots
set(gcf,'Position',[100 100 2400 900]); % Set figure size
PlotSwitch=1;
% Encapsulating the annealing process in a function would be beneficial
% Example: [KannealG, PiniG, LldG] = performAnnealing(...);

%% Analyze Data at each time point
for T_i = 1:TimeNum
    %% Data processing for each time point
    DataTime=DataCell{1,T_i};
    %% Setting up boundaries and initializing parameters
    % Consider encapsulating this into a function for clarity
    [Nbar,Kstart] = setupAnnealingParameters(DataTime,Edge,LabelStrength,Tstay,Model);%%Model Fix
    K0=Kstart(UnLockIndex);
    if T_i>1
        Kstart(UnLockIndex)=KannealG(T_i-1,:);%%K start anneal
        K0=KannealG(T_i-1,:);
    end
    %% Initialize parameters and handle inheritance for RNA distribution
    figure(hInherit)
    nexttile

    LabelWait=Label;
    Tstart=TimeLabel(1)-60;
    TimeElong = [TimeLabel(T_i) - Tstay, TimeLabel(T_i)];
    if TimeElong(1)<Tstart
        TimeAtK=Tstart-TimeElong(1);
        for ii =3:-1:1
            TimeAtK=TimeAtK-LabelWait(ii);
            LabelUse(ii)=LabelWait(ii)+min([TimeAtK,0]);
            if TimeAtK<0
                break
            end
        end
        LabelWait=LabelWait-LabelUse;
        TimeElong(1)=Tstart;
        KannealTime(T_i,1)=Tstart;
    end

    [PiniStart, LabelWait, Tfit, DtdInherit,timeBoundaries, featureValues,Kindex] = initializeParametersAndInherit...
        (T_i,T_i,TimeElong, TimeLabel, PiniG, KannealG, KannealTime, InheritSwitch,InheritModel, LabelWait, TimeElong(2)-TimeElong(1), Kstart, UnLockIndex,PlotSwitch);
    %% Analysis and fitting process
    % This is a core part of the code and could be a separate function
    LabelPart=LabelWait/sum(LabelWait);Turn=0;LldPilot=inf;
    while Turn < RepeatNum % repeat turn
        Turn = Turn + 1;
        [Kanneal] = performAnnealing(K0,DataTime, Kstart, PiniStart,LabelPart,UnLockIndex,Tfit,DtdInherit,Lb,Ub);%%Model Fix
        %%% Lld&Pini
        [lld,~,~,pini,~] =  FSP_FreeSet(Kanneal,K0,DataTime,PiniStart,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit,0);%%Model Fix
        if lld<LldPilot
            KannealG(T_i,:)=Kanneal;
            
            PiniG(T_i,:)=pini;
            LldPilot=lld;
            KannealOut=Kanneal;
        end
    end
    LldG(T_i,:)=LldPilot;
    PiniS(T_i,:)=PiniStart;
    [~,Pdtb,~,Pinis,~] =  FSP_FreeSet(KannealOut,K0,DataTime,PiniStart,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit,0);%%Model Fix
    [~,PdtbInherit,~,~,~] =  FSP_FreeSet(KannealOut,K0,DataTime,PiniStart,LabelPart,0,Kstart,UnLockIndex,0,DtdInherit,0);%%Model Fix
    [~,Pdtb_steady,~,~,~] =  FSP_FreeSet_Steady(KannealOut,DataTime,PiniStart,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit,0);
    MeanExpFitG(T_i,:)=[sum(Pdtb(1:71).*[0:1:70],2) sum(Pdtb_steady(1:71).*[0:1:70],2)];%[Mean MeanInSteady]
    ExpFitG(T_i,:)=[sum(Pdtb(2:end),2) sum(Pdtb_steady(2:end),2)];%[Mean MeanInSteady]Logic
    %% Visualization and result storage
    % Encapsulate visualization into a separate function
    visualizeDataandfit(Hanneal, Nbin, Nbar, [Pdtb;PdtbInherit], TimeLabel, T_i);

    % ... rest of the code for the loop ...
end

%% Saving and plotting results
% These steps can also be encapsulated for better readability and reuse
% Example: saveResults(OutFolder, KannealG, PiniG, LldG, ...);
figure(Hanneal)
saveas(gcf,[OutFolder,'FitAnneal','.fig']);
saveas(gcf,[OutFolder,'FitAnneal','.png']);
figure(hInherit)
saveas(gcf,[OutFolder,'TimeInherit','.fig']);
saveas(gcf,[OutFolder,'TimeInherit','.png']);

%% k-transform fig
switch InheritModel
    case 'Near'
        TimeLabelPlot=TimeLabel;KannealGPlot=KannealG;
    case 'Label'
        [timeBoundaries, featureValues] = createTimelines(KannealTime, KannealG);
        [~, PstatePlot] = createTimelines(KannealTime, PiniG);
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
save([OutFolder,['Results','.mat']],'KannealG','PiniG','PiniS','LldG','Lb','Ub','UnLockIndex','TimeLabelPlot','KannealGPlot','PstatePlot')

figure%Plot Kon\Koff\Kini vs. Time
title('Anneal K')
yyaxis left
hold on
plot(TimeLabelPlot,KannealGPlot(:,1),'-.c','LineWidth',2)
plot(TimeLabelPlot,KannealGPlot(:,2),'-.m','LineWidth',2)
ylim([0 0.02]);ylabel('K');xlabel('Evolution Time');xlim([0 max(TimeLabelPlot)-1])
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
