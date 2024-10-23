function PonT=FuncPon(Cycle,Time)
    load('Y:\TimeFit\OutputModel\Ton\CyclePlot\-TimeWindowMingle-PonBoole0-0726-bcd-ori\FittingResults.mat')
    % load('Y:\TimeFit\OutputModel\Ton\CyclePlot\-TimeWindowMingle-PonBoole0-0703-bcd\FittingResults.mat')
    switch Cycle
        case '11'
            FucPon=Fuc.cycle11;
        case '12'
            FucPon=Fuc.cycle12;
        case '13'
            FucPon=Fuc.cycle13;
    end
    PonT=FucPon(Time);
end