function RNAproduce_Live_Simple(Annotation)
%%% A program to simulation the rna producing.

%% SavePath Initialization
tic;%clock
fold_save=['Y:\Model\RNAproduce_Live_Simple\',...
   Annotation,'\'];
mkdir(fold_save);

%% Parameter Settings
TimeEvo=900;%Evolution time
KonF=@(t) 0.008.*exp(-(t-140).^2./(2.*5e3));%Kon function
KoffF=@(t) 0.005.*stepfun(t,0);%Koff function
syms RnaLabelF(ELt)
RnaLabelF(ELt)=piecewise(...
    ELt<=20, 0,...
    (20<ELt)<=80, (ELt-20)/60,...
    ELt>80, 1);
Kini=1.5e-1;
PlotDetail=0;%Plot switch
EnsembleNum=1e3;
TimeResolution=1;
Tstart=0;
Tend=TimeEvo;%Time of detect.
Tstay=2*60;%Stay time of RNA at allele.
DetectTime=(0:TimeResolution:Tend);

%% Overview of parameter Settings
figure%Plot Kon\Koff\Expression vs. Time
subplot(1,2,1);title('Overview K')
yyaxis left
hold on
plot(DetectTime,KonF(DetectTime),'-.c','LineWidth',2)
plot(DetectTime,KoffF(DetectTime),'-.m','LineWidth',2)
ylim([0 0.01]);ylabel('K');xlabel('Evolution Time');xlim([0 TimeEvo])
yyaxis right
plot(DetectTime,KonF(DetectTime)./(KoffF(DetectTime)+KonF(DetectTime)),'LineWidth',2)
ylim([0 1]);ylabel('Steady Expression');legend('Kon','Koff','ActiveRatio')
subplot(1,2,2);
yyaxis left
plot(0:Tstay,RnaLabelF(0:Tstay),'-m','LineWidth',2)
xlim([0 Tstay]);xlabel('Elong-Time')
ylabel('Label Signal');title('Rna Label')
yyaxis right
plot([0 Tstay],[Kini Kini],'LineWidth',2)
ylabel('Kini')
saveas(gcf,[fold_save,'Overview K setting','.fig']);
saveas(gcf,[fold_save,'Overview K setting','.png']);

%% Simulation-Ense
RsignalEnse=[];%Ini rna number ensemble
parfor_progress(EnsembleNum);
PoolCreate(20);
parfor ensmble_n=1:EnsembleNum
    P=[];MeanOFF=KoffF(0)/(KonF(0)+KoffF(0));
    P(1,1)=stepfun(rand/MeanOFF,1);
    R=[];R(1,1)=0;
    Time_detect=[];Time_detect(1,1)=0;
    time_sum=[];time_sum(1,1)=0;
    Index=1;
    %% Simulation-allele
    while time_sum(1,end)<Tend
        Index=Index+1;
        K_tot=[KonF(time_sum(1,Index-1)) KoffF(time_sum(1,Index-1))];
        if P(1,Index-1)==0, K_tot(2)=0; K_ini=0;
        elseif P(1,Index-1)==1, K_tot(1)=0; K_ini=Kini;end

        k_tot=sum([K_tot,K_ini],2);
        k_mat=[K_tot,K_ini];
        k_ratio=arrayfun(@(x) sum(k_mat(1,1:x),2),1:3);

        Time_tau=1;
        Time_detect(1,Index)=Time_tau;
        time_sum(1,Index)=time_sum(1,Index-1)+Time_tau;
        rand_num=rand;
        rand_index=find(sort([k_ratio,rand_num])==rand_num);
        if rand_index==1
            P(:,Index) = P(:,Index-1);R(:,Index) = zeros(1,1);
            P(1,Index) = 1;
        elseif rand_index==2
            P(:,Index) = P(:,Index-1);R(:,Index) = zeros(1,1);
            P(1,Index) = 0;
        elseif rand_index==3
            P(:,Index) = P(:,Index-1);R(:,Index) = zeros(1,1);
            R(1,Index) = 1;
        elseif rand_index==4
            P(:,Index) = P(:,Index-1);R(:,Index) = zeros(1,1);
        end
    end
    %% Calculate Rna Signal
    Rsignal=zeros(1,size(DetectTime,2));
    for T_i=DetectTime
        SignalStayIndex=time_sum<T_i&time_sum>=T_i-Tstay&R==1;
        ElongTime=T_i-time_sum;
        SignalStay=double(RnaLabelF(ElongTime));
        SignalNum=sum(SignalStay(SignalStayIndex),2);
        Rsignal(1,T_i+1)=SignalNum;
    end
    RsignalEnse(ensmble_n,:)=Rsignal;
    parfor_progress;
end
save([fold_save,['RnaSimu','.mat']],'RsignalEnse')

%% Plot Signal Bar at selected time
MaxRnaNum=100;Edge=-0.5:1:MaxRnaNum+0.5;
figure;
set(gcf,'Position',[200 300 2200 600])
Index=0;
for T_i=(1:3:15)*60
    Index=Index+1;
    subplot(1,5,Index)
    SignalsAtTime=RsignalEnse(:,T_i);
    Nbar = histcounts(SignalsAtTime,Edge)/EnsembleNum;
    h=bar(0:1:MaxRnaNum,Nbar,'hist');
    h.FaceColor=[0.3010 0.7450 0.9330];
    ylim([0 0.1]);xlim([-0.5 30])
    axis square
    title(['EvolutionTime: ',num2str(T_i/60),'min'])
end
saveas(gcf,[fold_save,'Overview Simulation Output','.png']);
%% clock
escape_time=toc;
disp(['Escape time: ',num2str(escape_time/3600,'%.3f'),' hours'])