function [Kanneal] = performAnnealing(K0,DataTime, Kstart, PiniStart,LabelPart,UnLockIndex,Tfit,DtdInherit,Lb,Ub)
    % Perform the annealing process
    % This function should encapsulate the core logic of your annealing process
    % Initialize variables and perform the annealing optimization
    % Return the optimized parameters and other metrics
    
    %% simulannealbnd and scan
    KstartAnneal=Kstart(UnLockIndex);%Select parameters to anneal
    LbAnneal=Lb(UnLockIndex);
    UbAnneal=Ub(UnLockIndex);
    % LabelPart=LabelWait/sum(LabelWait);
    [Kanneal,~,~] = simulannealbnd(@(x) FSP_FreeSet(x,K0,DataTime,PiniStart,LabelPart,0,Kstart,UnLockIndex,Tfit,DtdInherit,1),...
        KstartAnneal,LbAnneal,UbAnneal);
end