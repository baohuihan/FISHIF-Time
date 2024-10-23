function ModuleFitTimeShow(DataCellMean,TimeLabelMean,DataCell,TimeLabel,RnaLenth,kfitTime,KannealG,PiniG,OutFolder,InheritSwitch,RepeatNum,Bound,Model)
%%% A program to fit the experiment data by model

%% DataIni
TimeLabel=TimeLabel*60;%min to sec
TimeLabelMean=TimeLabelMean*60;
ShowTimeBin=1;%sec
Tmax=max(TimeLabel)+60;Tmin=0;
TimeShow=[Tmin:ShowTimeBin:Tmax,TimeLabel,TimeLabelMean];
TimeShowLabel=[zeros(1,size(Tmin:ShowTimeBin:Tmax,2)),ones(1,size(TimeLabel,2)),2*ones(1,size(TimeLabelMean,2))];
[TimeShowSort,Ilabel] = sort(TimeShow);
TimeShowLabelSort=TimeShowLabel(Ilabel);
DataIndex=find(TimeShowLabelSort>0);DataIndex=[DataIndex,DataIndex(end)];
DataIndex_i=0;

UnLockIndex=logical(Bound(1,:));
Ub=Bound(2,:);
Lb=Bound(3,:);

TimeNum=size(TimeShowSort,2);
Index=(1:TimeNum);IndexK=1:size(TimeLabelMean,2);
DataLabel={'single','muti'};DataLabelColor={'k','blue','red'};InheritText={'OFF','ON'};
DelayTime=[0.1,0.5,3];
CatLim=[89 300];
Ht=figure;
%% Initialize required parameters
KFitG=[]; MeanExpFitG=[];
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
KannealTime=TimeLabelMean'-[Tstay kfitTime];
Hanneal=figure; % Initialize a figure for annealing plots
set(gcf,'Position',[100 100 2200 900]); % Set figure size
IndexData1=0;IndexData2=0;
DataTime=0;KIStart=0;KIStart0=0;NextDistance=1;
%% Analyze Data at time i
for T_i = 1:TimeNum
    %% Data i
    switch TimeShowLabelSort(T_i)
        case 1
            IndexData1=IndexData1+1;
            DataTime=DataCell{IndexData1,1};
            DataIndex_i=DataIndex_i+1;
        case 2
            IndexData2=IndexData2+1;
            DataTime=DataCellMean{IndexData2,1};
            KIStart=KIStart+1;
            DataIndex_i=DataIndex_i+1;
    end
    Nbar = histcounts(DataTime,Edge)/size(DataTime,1);
    if DataIndex_i>0
        NextDistance=(DataIndex(DataIndex_i+1)-T_i)/(DataIndex(DataIndex_i+1)-DataIndex(DataIndex_i));
        if isnan(NextDistance)||isinf(NextDistance)
            NextDistance=1;
        end
    end
    try
        switch TimeShowLabelSort(DataIndex(DataIndex_i+1))
            case 1
                DataTimeNext=DataCell{IndexData1+1,1};
            case 2
                DataTimeNext=DataCellMean{IndexData2+1,1};
        end
    catch
        DataTimeNext=0;
    end
    NbarNext = histcounts(DataTimeNext,Edge)/size(DataTimeNext,1);

    %% Initialize required parameters
    DtdInherit=[];
    %% analysis start
    if 1
        %% calulate the initial k by poission fit & Get ini dtb
        Kstart=[2e-2,0,2e-2,0,0,0,0,0.3,0,0,1,0.5];      
            TimeElong=[TimeShowSort(T_i)-Tstay TimeShowSort(T_i)];%Elong time interval
            Kevo=min([max(Index(TimeShowSort(T_i)>[-inf TimeLabelMean inf])) length(TimeLabelMean)]);
            figure
            [PiniStart, LabelWait, Tfit, DtdInherit,timeBoundaries, featureValues,Kindex] = initializeParametersAndInherit...
        (Kevo,Kevo,TimeElong, TimeLabelMean, PiniG, KannealG, KannealTime, InheritSwitch,InheritModel, Label, Tstay, Kstart, UnLockIndex);
        %     [PiniStart, LabelWait, Tfit, DtdInherit,timeBoundaries, featureValues,Kindex] = initializeParametersAndInherit...
        % (T_i,length(TimeLabelMean),TimeElong, TimeLabelMean, PiniG, KannealG, KannealTime, InheritSwitch,InheritModel, Label, Tstay, Kstart, UnLockIndex);
            close gcf
             % [~, featureValues] = createTimelines(KannealTime(1:IndexData2+1,:), KannealG(1:IndexData2+1,:));
            % Kend=min([max([1 max(Kindex)+1]),length(featureValues)]);
            LabelPart=LabelWait/sum(LabelWait);
            [~,Pdtb,~,Pinis,~] =  FSP_FreeSet(KannealG(Kevo,:),0,PiniStart,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit);
            PdtbG(T_i,:)=Pdtb;
        %% data and fit curve figure
        figure(Ht);
        h=bar(Nbin,Nbar,'hist');h.FaceColor=[0.3010 0.7450 0.9330];
        h.FaceAlpha=NextDistance;h.EdgeAlpha=NextDistance;
        hold on
        h=bar(Nbin,NbarNext,'hist');h.FaceColor=[0.3010 0.7450 0.9330];
        h.FaceAlpha=1-NextDistance;h.EdgeAlpha=1-NextDistance;
        ylim([0 0.05]);xlim([-0.5 70])
        plot(0:70,Pdtb(1:71),'r','LineWidth',2)
        xlabel('RNA Number')
        ylabel('Frequency')
        axis square
        title([num2str(TimeShowSort(T_i)),'s '],'Color',DataLabelColor{TimeShowLabelSort(T_i)+1})
        F=getframe(gcf);
        I=frame2im(F);
        [I,map]=rgb2ind(I,256);
        if T_i == 1
            imwrite(I,map,[OutFolder,'TimeEvo',InheritText{InheritSwitch+1},'.gif'],'gif','Loopcount',inf,'DelayTime',DelayTime(TimeShowLabelSort(T_i)+1));
        else
            imwrite(I,map,[OutFolder,'TimeEvo',InheritText{InheritSwitch+1},'.gif'],'gif','WriteMode','append','DelayTime',DelayTime(TimeShowLabelSort(T_i)+1)); 
        end
        if T_i>=CatLim(1)&&T_i<=CatLim(2)
            if T_i == CatLim(1)
                imwrite(I,map,[OutFolder,'TimeCatEvo',InheritText{InheritSwitch+1},'.gif'],'gif','Loopcount',inf,'DelayTime',DelayTime(TimeShowLabelSort(T_i)+1));
            else
                imwrite(I,map,[OutFolder,'TimeCatEvo',InheritText{InheritSwitch+1},'.gif'],'gif','WriteMode','append','DelayTime',DelayTime(TimeShowLabelSort(T_i)+1)); 
            end
        end

        clf(Ht,'reset')
        
        if TimeShowLabelSort(T_i)~=0
            figure(Hanneal)
            nexttile
            h=bar(Nbin,Nbar,'hist');h.FaceColor=[0.3010 0.7450 0.9330];
            hold on
            ylim([0 0.05]);xlim([-0.5 70])
            hold on
            plot(0:70,Pdtb(1:71),'b','LineWidth',2)
    %         plot(0:70,Pdtb_Hold(1:71),'r','LineWidth',2)
            xlabel('RNA Number')
            ylabel('Frequency')
            axis square
            title([num2str(TimeShowSort(T_i)),'s ','DataEmbryo ',DataLabel{TimeShowLabelSort(T_i)}],'Color',DataLabelColor{TimeShowLabelSort(T_i)})
            DdtbG(T_i,:)=Nbin;
        end
    end       
end
sgtitle(['Inherit',InheritText{InheritSwitch+1}])
saveas(gcf,[OutFolder,'FitAnnealShowInherit',InheritText{InheritSwitch+1},'.fig']);
saveas(gcf,[OutFolder,'FitAnnealShowInherit',InheritText{InheritSwitch+1},'.png']);
%% save parameter file
save([OutFolder,['ResultsShowInherit',InheritText{InheritSwitch+1},'.mat']],'PdtbG','PiniG','DdtbG','TimeShowSort','TimeShowLabelSort')
%% k-transform fig
% figure%Plot Kon\Koff\Kini vs. Time
% title('Anneal K')
% yyaxis left
% hold on
% plot(TimeShowSort,KannealG(:,1),'-oc','LineWidth',2)
% plot(TimeShowSort,KannealG(:,2),'-om','LineWidth',2)
% plot(TimeShowSort,KannealG(:,3),'-ok','LineWidth',2)
% plot(TimeShowSort,KannealG(:,4),'-or','LineWidth',2)
% ylim([0 0.2]);ylabel('K');xlabel('Evolution Time');xlim([0 max(TimeShowSort)+60])
% yyaxis right
% plot(TimeShowSort,KannealG(:,5),'LineWidth',2)
% ylim([0 1]);ylabel('Kini');legend('ON','OFF','Kon','Koff','Kini')
% saveas(gcf,[OutFolder,'KAnneal','.fig']);
% saveas(gcf,[OutFolder,'KAnneal','.png']);
% 
% figure%Plot Kon/Koff vs. Time
% title('Anneal Kon/Koff')
% plot(TimeShowSort,KannealG(:,1)./KannealG(:,2),'-.c','LineWidth',2)
% ylabel('Kon/Koff');xlabel('Evolution Time');xlim([0 max(TimeShowSort)+60])
% saveas(gcf,[OutFolder,'KratioAnneal','.fig']);
% saveas(gcf,[OutFolder,'KratioAnneal','.png']);
% 
% figure%Plot double vs. Time
% title('Anneal Alp')
% hold on
% plot(TimeShowSort,KannealG(:,end),'LineWidth',2)
% ylim([0 1]);ylabel('Alp');xlabel('Evolution Time');xlim([0 max(TimeShowSort)+60])
% legend('Alp')
% saveas(gcf,[OutFolder,'AlpAnneal','.fig']);
% saveas(gcf,[OutFolder,'AlpAnneal','.png']);
% 
% figure%Plot Pstate vs. Time
% title('Anneal Pstate')
% hold on
% plot(TimeShowSort',PiniG,'-o','LineWidth',2)
% ylim([0 1]);ylabel('Pstate');xlabel('Evolution Time');xlim([0 max(TimeShowSort)+60])
% legend('Soff','S0','S1')
% saveas(gcf,[OutFolder,'PstateAnneal','.fig']);
% saveas(gcf,[OutFolder,'PstateAnneal','.png']);