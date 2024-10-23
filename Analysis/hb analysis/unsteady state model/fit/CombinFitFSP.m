function [LLd_all,PDtbG,PiniG0,KannealTime]=CombinFitFSP(KG,TimeCombin,DataCell,Edge,LabelStrength,Tstay,Model, TimeLabel, PiniG, KannealG, KannealTime, InheritSwitch,InheritModel, Label, UnLockIndex)
TimeNum=size(TimeCombin,2);
Kl=sum(UnLockIndex,2);
LLd_all=0;
for T_ic = 1:TimeNum
    T_i=TimeCombin(T_ic);
    Kx=KG(1,(1:Kl)+(T_ic-1)*Kl);
    %% Data processing for each time point
    DataTime=DataCell{1,T_i};
    %% Setting up boundaries and initializing parameters
    % Consider encapsulating this into a function for clarity
    [~,Kstart] = setupAnnealingParameters(DataTime,Edge,LabelStrength,Tstay,Model);%%Model Fix
    K0=Kstart(UnLockIndex);
    Kdiff=1;
    if T_i>1
        Kstart(UnLockIndex)=KannealG(T_i-1,:);%%K start anneal
        K0=KannealG(T_i-1,:);
        Kdiff=1;
    end
    %% Initialize parameters and handle inheritance for RNA distribution
    LabelWait=Label;
    Tstart=min([TimeLabel(1)-60 KannealTime(1,2)-1]);
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
    % TimeElong = [TimeLabel(T_i) - Tstay, TimeLabel(T_i)];
    [PiniStart, LabelWait, Tfit, DtdInherit,timeBoundaries, featureValues,Kindex] = initializeParametersAndInherit...
        (T_i,T_i,TimeElong, TimeLabel, PiniG, KannealG, KannealTime, InheritSwitch,InheritModel, LabelWait, TimeElong(2)-TimeElong(1), Kstart, UnLockIndex,0);
    %% Analysis and fitting process
    % This is a core part of the code and could be a separate function
    LabelPart=LabelWait/sum(LabelWait);
    % Turn=0;LldPilot=inf;
    % while Turn < RepeatNum % repeat turn
    %     Turn = Turn + 1;
    [lld,P_dtb_rna,~,pini,~] =  FSP_FreeSet(Kx,K0,DataTime,PiniStart,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit,Kdiff);%%Model Fix
    PDtbG(T_ic,:)=P_dtb_rna;
    PiniG0(T_ic,:)=pini;
    PiniG(T_i,:)=pini;
    KannealG(T_i,:)=Kx;
    LLd_all=LLd_all+lld;
end