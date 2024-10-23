function ModuleFitTimeSimu
%%% A program to fit the experiment data by model

%% DataIni
OutFolder='Y:\TimeFit\Output\Simu\Ense-1e3 T900 Tr120-3\';
mkdir(OutFolder);
load('Y:\Model\RNAproduce_Live_Simple\Ense-1e3 T900\RnaSimu.mat');%Simulation Data
TimeResolution=120;
TimeLabel=1:TimeResolution:size(RsignalEnse,2);
TimeNum=size(TimeLabel,2);Index=(1:TimeNum);
DataCell=mat2cell(RsignalEnse(:,TimeLabel),size(RsignalEnse,1),ones(1,size(RsignalEnse(:,TimeLabel),2)));
MeanPlot(TimeLabel',DataCell','')
xlim([0 size(RsignalEnse,2)]);xlabel('Time');ylabel('MeanExpr')
saveas(gcf,[OutFolder,'Mean-Time','.png']);

%% Initialize required parameters
KFitG=[];
FitSwitch=1;
MaxRnaNum=100;Edge=-0.5:1:MaxRnaNum+0.5;
Nbin=0:1:MaxRnaNum;
Pini=[1 0 0]';
Tstay=2*60;
Label=[20 60 40];LabelSum=cumsum(Label);
LabelPart=[20 60 40]/Tstay;
LabelStrength=sum([0 0.5 1].*LabelPart);
KannealG=[];PiniG=[];LldG=[];
Hanneal=figure;
set(gcf,'Position',[100 100 2200 900])
%% Analyze Data at time i
for T_i = 1:TimeNum
    %% Data i
    DataTime=DataCell{1,T_i};
    %% Set Up\Low boundary
    %k=[k01 k02 k10 k12 k20 k21 ki0 ki1 ki2 useless0 dataintensity1 alpha0]
    UnLockIndex=logical([1 0 1 0 0 0 0 1 0 0 0 0]);
    Ub=[0.02,0,0.02,0,0,0,0,0.5,0,0,1,0];
    Lb=[0,0,0,0,0,0,0,0,0,0,1,0];

    %% Initialize required parameters
    DtdInherit=[];
    %% analysis start
    if FitSwitch
        %% poission fit
        Nbar = histcounts(DataTime,Edge)/size(DataTime,1);
        state0=Nbar(1,1);
        figure
        h=bar(Nbin,Nbar,'hist');h.FaceColor=[0.3010 0.7450 0.9330];
        hold on
        [para,poiss2_dtb_best2] = PoissonFitScan_singledata(DataTime,state0);
        poiss2_dtb_best2(1)=state0(1);
        plot(0:100,poiss2_dtb_best2,'LineWidth',1);
        legend('Data','Poisson')
        ylim([0 0.1]);xlim([0 40])
%         saveas(gcf,[OutFolder,['Pioss-Time',num2str(TimeLabel(T_i)),'.fig']]);
%         saveas(gcf,[OutFolder,['Pioss-Time',num2str(TimeLabel(T_i)),'.png']]);
        close(gcf);
        
        %% repeat the simulannealbnd to find the best one
        LldPilot=inf;
        Turn=0;
        %% calulate the initial k by poission fit & Get ini dtb
        Kstart=[5e-2*para(1),0,5e-2*(1-para(1)),0,0,0,0,para(2)/(LabelStrength.*Tstay),0,0,1,0];
        if T_i==1
            Pinie=Pini;
            LabelWait=Label;
            Tfit=Tstay;
        else
            TimeElong=[TimeLabel(T_i)-Tstay TimeLabel(T_i)];%Elong time interval
            EvoKsIndex=TimeLabel>TimeElong(1)&TimeLabel<TimeElong(2);%Ks need to evo
            Kindex=Index(EvoKsIndex);
            TimeNode=sort([TimeLabel(EvoKsIndex),TimeElong]);
            if isempty(Kindex)
                PiniStart0=PiniG(T_i-1,:)';
                TimeStart0=TimeLabel(T_i-1);
                KIStart=T_i-1;
            elseif Kindex(1)>1
                PiniStart0=PiniG(Kindex(1)-1,:)';%Pini
                TimeStart0=TimeLabel(Kindex(1)-1);
                KIStart=Kindex(1);
            else
                PiniStart0=[1 0 0]';
                TimeStart0=-Tstay;
                KIStart=Kindex(1);
            end
            [~,~,~,PiniStart,~] =  FSP_FreeSet(KannealG(KIStart,:),0,PiniStart0,LabelPart,0,Kstart,UnLockIndex,TimeElong(1)-TimeStart0,[]);
            EvoI=0;
            LabelUse=[];
            LabelWait=Label;
            for K_i=Kindex 
                EvoI=EvoI+1;
                TimeAtK=(TimeLabel(K_i)-TimeNode(EvoI));
                for ii =3:-1:1
                    TimeAtK=TimeAtK-LabelWait(ii);
                    LabelUse(ii)=LabelWait(ii)+min([TimeAtK,0]);
                    if TimeAtK<0
                        break
                    end
                end
                LabelWait=LabelWait-LabelUse;
                LabelPart=LabelUse/sum(LabelUse);
                [~,~,~,~,DtdInherit] =  FSP_FreeSet(KannealG(K_i,:),0,PiniStart,LabelPart,0,Kstart,UnLockIndex,TimeAtK,DtdInherit);
            end
            Tfit=TimeElong(end)-TimeElong(end-1);
        end
        while Turn < 1  % repeat turn
            Turn = Turn + 1;
            %% simulannealbnd and scan
            KstartAnneal=Kstart(UnLockIndex);%Select parameters to anneal
            LbAnneal=Lb(UnLockIndex);
            UbAnneal=Ub(UnLockIndex);
            LabelPart=LabelWait/sum(LabelWait);
            [Kanneal,~,~] = simulannealbnd(@(x) FSP_FreeSet(x,DataTime,Pinie,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit),...
                KstartAnneal,LbAnneal,UbAnneal);
            %% Lld&Pini
            [lld,~,~,pini,~] =  FSP_FreeSet(Kanneal,DataTime,Pinie,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit);
            if lld<LldPilot
                KannealG(T_i,:)=Kanneal;
                PiniG(T_i,:)=pini;
                LldPilot=lld;
                KannealOut=Kanneal;
            end
            LldG(T_i,:)=LldPilot;
            [~,Pdtb,~,Pinis,~] =  FSP_FreeSet(KannealOut,DataTime,Pinie,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit);%Get Distribution
        end
        %% data and fit curve figure
        figure(Hanneal)
        nexttile
        h=bar(Nbin,Nbar,'hist');h.FaceColor=[0.3010 0.7450 0.9330];
        hold on
        ylim([0 0.01]);xlim([-0.5 40])
        hold on
        plot(0:50,Pdtb(1:51),'b','LineWidth',2)
        xlabel('RNA Number')
        ylabel('Frequency')
        axis square
        title([num2str(TimeLabel(1,T_i)),'s'])
    end       
end
saveas(gcf,[OutFolder,'FitAnneal','.fig']);
saveas(gcf,[OutFolder,'FitAnneal','.png']);
%% save parameter file
save([OutFolder,['Results','.mat']],'KannealG','PiniG','LldG')
%% k-transform fig
figure%Plot Kon\Koff\Kini vs. Time
title('Anneal K')
yyaxis left
hold on
plot(TimeLabel,KannealG(:,1),'-.c','LineWidth',2)
plot(TimeLabel,KannealG(:,2),'-.m','LineWidth',2)
ylim([0 0.01]);ylabel('K');xlabel('Evolution Time');xlim([0 max(TimeLabel)-1])
yyaxis right
plot(TimeLabel,KannealG(:,3),'LineWidth',2)
ylim([0 0.5]);ylabel('Kini');legend('Kon','Koff','Kini')
saveas(gcf,[OutFolder,'KAnneal','.fig']);
saveas(gcf,[OutFolder,'KAnneal','.png']);

figure%Plot Kon/Koff vs. Time
title('Anneal Kon/Koff')
plot(TimeLabel,KannealG(:,1)./KannealG(:,2),'-.c','LineWidth',2)
ylim([0 2]);ylabel('Kon/Koff');xlabel('Evolution Time');xlim([0 max(TimeLabel)-1])
saveas(gcf,[OutFolder,'KratioAnneal','.fig']);
saveas(gcf,[OutFolder,'KratioAnneal','.png']);