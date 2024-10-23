function [KannealG,MeanExpFitG,PiniG]=ModuleFitTime(RnaSignal,TimeLabel,OutFolder,InheritSwitch,Model)
%%% A program to fit the experiment data by model
KannealG=[];MeanExpFitG=[];PiniG=[];
switch Model
    case 'TwoState'
        [KannealG,MeanExpFitG,PiniG]=ModuleFitTime_TwoState(RnaSignal,TimeLabel,OutFolder,InheritSwitch);
end