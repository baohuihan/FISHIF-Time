function ModuleFitTime_Show(DataCellMean,TimeLabelMean,DataCell,TimeLabel,KannealG,PiniG,OutFolder,InheritSwitch)
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

TimeNum=size(TimeShowSort,2);
Index=(1:TimeNum);IndexK=1:size(TimeLabelMean,2);
DataLabel={'single','muti'};DataLabelColor={'k','blue','red'};InheritText={'OFF','ON'};
DelayTime=[0.1,0.2,0.5];
CatLim=[89 300];
%% Initialize required parameters
FitSwitch=1;
% InheritSwitch=1;
RepeatNum=1;
IniModel=1;
%%%%%%%%%%%
KFitG=[];PdtbG=[];DdtbG=[];
MaxRnaNum=100;Edge=-0.5:1:MaxRnaNum+0.5;
Nbin=0:1:MaxRnaNum;
Pini=[1 0 0]';
RnaLenth=[497+2725+3+283 2277 407+509];%CDS
% RnaLenth=[0 0 497+2725+3+283+2277+407+509];
Tstay=sum(RnaLenth)/25;% 25bp/s
Label=RnaLenth/25;
LabelSum=cumsum(Label);
LabelPart=Label/Tstay;
LabelStrength=sum([0 0.5 1].*LabelPart);
% PiniG=[];
LldG=[];MeanExpFitG=[];
Hanneal=figure;
set(gcf,'Position',[100 100 2200 900])
Ht=figure;
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
    %% Set Up\Low boundary
    %k=[k01 k02 k10 k12 k20 k21 ki0 ki1 ki2 Ton dataintensity1 alpha0]
    UnLockIndex=logical([1 0 1 1 0 1 0 0 1 0 0 1]);
%     Ub=[0.1,0,0.1,0.2,0,0.01,0,0,0.55,0,1,0.5];
%     Lb=[0,0,0,0.001,0,0.0001,0,0,0.15,0,1,0.5];

    %% Initialize required parameters
    DtdInherit=[];
    %% analysis start
    if FitSwitch
        %% calulate the initial k by poission fit & Get ini dtb
        Kstart=[2e-2,0,2e-2,5e-2,0,5e-2,0,0,0.45,0,1,0.5];
%             figure
%             plot([Tmin Tmax],[0 0],'LineWidth',2)
%             hold on
%              
            TimeElong=[TimeShowSort(T_i)-Tstay TimeShowSort(T_i)];%Elong time interval
            EvoKsIndex=TimeLabelMean>TimeElong(1)&TimeLabelMean<TimeElong(2);%Ks need to evo
            EvoStartK=IndexK(max(TimeLabelMean<TimeElong(1))+1);
            Kindex=IndexK(EvoKsIndex);
            TimeNode=sort([TimeLabelMean(EvoKsIndex),TimeElong]);%K change nodes
            switch IniModel
                case 1
                    if EvoStartK==1
                        PiniStart0=[1 0 0]';
                        TimeStart0=0;
                    else
                        PiniStart0=PiniG(EvoStartK-1,:)';
                        TimeStart0=TimeLabelMean(EvoStartK-1);
                    end
                    TevoS=max([TimeElong(1)-TimeStart0,0]);
                    [~,~,~,PiniStart,~] =  FSP_FreeSet20(KannealG(EvoStartK,:),0,PiniStart0,LabelPart,0,Kstart,UnLockIndex,TevoS,[]);%%MaybeErrorExist
                case 2
                    PiniStart=PiniG(max(Kindex),:)';
            end
            
            LabelWait=Label;Tfit=Tstay;
            if InheritSwitch==1%Inherit of Rna distribution
                EvoI=0;
                LabelUse=[];
                DtdInherit=[];
                for K_i=Kindex 
                    EvoI=EvoI+1;
                    TimeAtK=(TimeLabelMean(K_i)-TimeNode(EvoI));
                    TimeAtKev=TimeAtK;
                    for ii =3:-1:1
                        TimeAtK=TimeAtK-LabelWait(ii);
                        LabelUse(ii)=LabelWait(ii)+min([TimeAtK,0]);
                        if TimeAtK<0
                            break
                        end
                    end
                    LabelWait=LabelWait-LabelUse;
                    LabelPart=LabelUse/sum(LabelUse);
                    [~,~,~,~,DtdInherit] =  FSP_FreeSet20(KannealG(K_i,:),0,PiniStart,LabelPart,0,Kstart,UnLockIndex,TimeAtKev,DtdInherit);
                end
                Tfit=TimeElong(end)-TimeNode(end-1);
            end
            
            Kend=min([max([1 max(Kindex)+1]),IndexK(end)]);
            [~,Pdtb,~,Pinis,~] =  FSP_FreeSet20(KannealG(Kend,:),0,PiniStart,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit);
            PdtbG(T_i,:)=Pdtb;
        %while Turn < RepeatNum  % repeat turn
            % if KIStart>KIStart0
            %     PiniG(KIStart,:)=Pinis;
            %     KIStart0=KIStart;
            % end
           
%             On=KannealOut(1)/(KannealOut(1)+KannealOut(2));KannealHoldOff=KannealOut;
%             KannealHoldOff(2)=0.05;KannealHoldOff(1)=0.05*On/(1-On);
% %             [~,Pdtb_Hold,~,~,~] =  FSP_FreeSetTon(KannealHoldOff,DataTime,Pinie,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit);
%             [~,Pdtb_steady,~,~,~] =  FSP_FreeSet20_Steady(KannealOut,DataTime,Pinie,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit);
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