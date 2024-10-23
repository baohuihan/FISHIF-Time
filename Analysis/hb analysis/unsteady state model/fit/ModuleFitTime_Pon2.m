function [KannealG,MeanExpFitG,PiniG]=ModuleFitTime_Pon2(RnaSignal,TimeLabel,OutFolder)
%%% A program to fit the experiment data by model

%% DataIni
TimeLabel=TimeLabel*60;%min to sec
TimeNum=size(TimeLabel,2);Index=(1:TimeNum);
DataCell=RnaSignal';
% MeanExp=MeanPlot(TimeLabel',DataCell','');
% DeleMax=MeanExp>45;%delet > 45
% TimeLabel=TimeLabel(~DeleMax);DataCell=DataCell(~DeleMax);
% MeanExp=MeanPlot(TimeLabel',DataCell','');
% xlabel('Time/min');ylabel('MeanExpr')
% hold on
% [~,TFrm] = rmoutliers(MeanExp(),"movmedian",30,"SamplePoints",TimeLabel);
% scatter(TimeLabel(TFrm),MeanExp(TFrm),500,'x','LineWidth',3)
% saveas(gcf,[OutFolder,'Mean-Time','.png']);
% DataCell=DataCell(~TFrm');
% TimeLabel=TimeLabel(~TFrm)';
%% Initialize required parameters
FitSwitch=1;
InheritSwitch=0;
RepeatNum=5;
IniModel=1;
%%%%%%%%%%%
KFitG=[];
MaxRnaNum=100;Edge=-0.5:1:MaxRnaNum+0.5;
Nbin=0:1:MaxRnaNum;
Pini=[1 0 0]';
RnaLenth=[3+283 2277 407+509];%CDS
Tstay=sum(RnaLenth)/25;% 25bp/s
Label=RnaLenth/25;
LabelSum=cumsum(Label);
LabelPart=Label/Tstay;
LabelStrength=sum([0 0.5 1].*LabelPart);
KannealG=[];PiniG=[];LldG=[];MeanExpFitG=[];
Hanneal=figure;
set(gcf,'Position',[100 100 2200 900])
%% Analyze Data at time i
for T_i = 1:TimeNum
    %% Data i
    DataTime=DataCell{1,T_i};
    %% Set Up\Low boundary
    %k=[k01 k02 k10 k12 k20 k21 ki0 ki1 ki2 Ton dataintensity1 alpha0]
    UnLockIndex=logical([1 0 1 0 0 0 0 1 0 1 0 1]);
    Ub=[0.02,0,0.5,0,0,0,0,0.55,0,1,1,1];
    Lb=[0.02,0,0.001,0,0,0,0,0.35,0,0,1,0];

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
%         Kstart=[5e-2*para(1)/(1-para(1)),0,5e-2,0,0,0,0,para(2)/(LabelStrength.*Tstay),0,0,1,0];
        Kstart=[5e-2,0,5e-2,0,0,0,0,para(2)/(LabelStrength.*Tstay),0,0,1,0];
        Change=[];
        if T_i==1%%First evo
            Pinie=Pini;
            LabelWait=Label;
            Tfit=Tstay;
        else
            TimeElong=[TimeLabel(T_i)-Tstay TimeLabel(T_i)];%Elong time interval
            EvoKsIndex=TimeLabel>TimeElong(1)&TimeLabel<TimeElong(2);%Ks need to evo
            Kindex=Index(EvoKsIndex);
            if ~isempty(Kindex)
                Change=[TimeLabel(Kindex)-(TimeLabel(T_i)-Tstay);...
                    KannealG(Kindex,4)'];
            end
            TimeNode=sort([TimeLabel(EvoKsIndex),TimeElong]);
            switch IniModel
                case 1
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
                        TimeStart0=0;
                        KIStart=Kindex(1);
                    end
                    if TimeElong(1)-TimeStart0<0
                        Pinie=PiniG(T_i-1,:)';
                    else
                        [~,~,~,Pinie,~] =  FSP_FreeSet(KannealG(KIStart,:),0,PiniStart0,LabelPart,0,Kstart,UnLockIndex,TimeElong(1)-TimeStart0,[]);%%MaybeErrorExist
                    end
                case 2
                    PiniStart=PiniG(max(Kindex),:)';
            end
            
            LabelWait=Label;
            if InheritSwitch==1%Inherit of Rna distribution
                EvoI=0;
                LabelUse=[];
                for K_i=Kindex 
                    EvoI=EvoI+1;
                    TimeAtK=(TimeLabel(K_i)-TimeNode(EvoI));
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
                    [~,~,~,~,DtdInherit] =  FSP_FreeSet(KannealG(K_i,:),0,PiniStart,LabelPart,0,Kstart,UnLockIndex,TimeAtKev,DtdInherit);
                end
                Tfit=TimeElong(end)-TimeNode(end-1);
            end
        end
        while Turn < RepeatNum  % repeat turn
            Turn = Turn + 1;
            %% simulannealbnd and scan
            KstartAnneal=Kstart(UnLockIndex);%Select parameters to anneal
            LbAnneal=Lb(UnLockIndex);
            UbAnneal=Ub(UnLockIndex);
            if T_i>1
                if TimeLabel(T_i)<7.5*60
                    LbAnneal(4)=KannealG(T_i-1,4);UbAnneal(4)=1;KstartAnneal(4)=KannealG(T_i-1,4);
                else
                    UbAnneal(4)=KannealG(T_i-1,4);LbAnneal(4)=0;KstartAnneal(4)=KannealG(T_i-1,4);
                end
            end
            LabelPart=LabelWait/sum(LabelWait);
            options = optimoptions('simulannealbnd','MaxFunctionEvaluations',40000);
            if UbAnneal(4)-LbAnneal(4)<=0.0002
                UbAnneal(4)=LbAnneal(4);
            end
            if ~isempty(Change)&&Change(2,1)~=sum(Pinie)
                Diff=Change(2,1)-sum(Pinie);
                Pinie=Pinie+Pinie./sum(Pinie).*Diff;
            end
            [Kanneal,~,~] = simulannealbnd(@(x) FSP_FreeSetPon2(x,DataTime,Pinie,LabelPart,0,Kstart,UnLockIndex,Tfit,Change,DtdInherit),...
                KstartAnneal,LbAnneal,UbAnneal,...
                 options);
            %% Lld&Pini
            [lld,dtb,~,pini,~] =  FSP_FreeSetPon2(Kanneal,DataTime,Pinie,LabelPart,0,Kstart,UnLockIndex,Tfit,Change,DtdInherit);
            if isnan(lld)
                continue
            end
            if lld<=LldPilot
                KannealG(T_i,:)=Kanneal;
                PiniG(T_i,:)=pini;
                LldPilot=lld;
                KannealOut=Kanneal;
            end
            LldG(T_i,:)=LldPilot;
            [~,Pdtb,~,Pinis,~] =  FSP_FreeSetPon2(KannealOut,DataTime,Pinie,LabelPart,0,Kstart,UnLockIndex,Tfit,Change,DtdInherit);%Get Distribution
            On=KannealOut(1)/(KannealOut(1)+KannealOut(2));KannealHoldOff=KannealOut;
            KannealHoldOff(2)=0.05;KannealHoldOff(1)=0.05*On/(1-On);
%             [~,Pdtb_Hold,~,~,~] =  FSP_FreeSetTon(KannealHoldOff,DataTime,Pinie,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit);
            [~,Pdtb_steady,~,~,~] =  FSP_FreeSetPon2_Steady(KannealOut,DataTime,Pinie,LabelPart,0,Kstart,UnLockIndex,Tfit,Change,DtdInherit);
            MeanExpFitG(T_i,:)=[sum(Pdtb(1:71).*[0:1:70],2) sum(Pdtb_steady(1:71).*[0:1:70],2)];%[Mean MeanInSteady]
        end
        %% data and fit curve figure
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
        title([num2str(TimeLabel(T_i)),'s'])
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
ylim([0 0.1]);ylabel('K');xlabel('Evolution Time');xlim([0 max(TimeLabel)+60])
yyaxis right
plot(TimeLabel,KannealG(:,3),'LineWidth',2)
ylim([0 1]);ylabel('Kini');legend('Kon','Koff','Kini')
saveas(gcf,[OutFolder,'KAnneal','.fig']);
saveas(gcf,[OutFolder,'KAnneal','.png']);

figure%Plot Kon/Koff vs. Time
title('Anneal Kon/Koff')
plot(TimeLabel,KannealG(:,1)./KannealG(:,2),'-.c','LineWidth',2)
ylabel('Kon/Koff');xlabel('Evolution Time');xlim([0 max(TimeLabel)+60])
saveas(gcf,[OutFolder,'KratioAnneal','.fig']);
saveas(gcf,[OutFolder,'KratioAnneal','.png']);

figure%Plot double vs. Time
title('Anneal Alp&Pon')
hold on
plot(TimeLabel,KannealG(:,4),'LineWidth',2)
plot(TimeLabel,KannealG(:,5),'LineWidth',2)
ylim([0 1]);ylabel('Alp');xlabel('Evolution Time');xlim([0 max(TimeLabel)+60])
legend('Pon','Alp')
saveas(gcf,[OutFolder,'AlpAnneal','.fig']);
saveas(gcf,[OutFolder,'AlpAnneal','.png']);