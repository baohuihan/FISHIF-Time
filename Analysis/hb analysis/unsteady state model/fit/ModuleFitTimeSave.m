function ModuleFitTimeSave(TimeLabelMean,EL,KannealELs,PiniELs,MeanExpELs,OutFolderMain,Model,Color)
%%% A program to save and show the fitting parameter
switch Model
    case 'TwoState'
        ModuleFitTimeSaveTwoState(TimeLabelMean,EL,KannealELs,PiniELs,MeanExpELs,OutFolderMain,Color);
    case 'Ton'
        ModuleFitTimeSaveTon(TimeLabelMean,EL,KannealELs,PiniELs,MeanExpELs,OutFolderMain,Color);
end