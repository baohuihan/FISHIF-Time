function [Nbar,Kstart] = setupAnnealingParameters(DataTime,Edge,LabelStrength,Tstay,Model)
    Nbar = histcounts(DataTime,Edge)/size(DataTime,1);
    state0=Nbar(1,1);
    [para,~] = PoissonFitScan_singledata(DataTime,state0);
    % Kstart=[1e-1*para(1),0,1e-1*(1-para(1)),0,0,0,0,para(2)/(LabelStrength.*Tstay),0,0,1,0.5];
    Kstart=[0.02,0,0.12,0,0,0,0,0.12,0,0,1,0];
end
   