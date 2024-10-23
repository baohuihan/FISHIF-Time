function [KannealG, PiniG] = ModifiedModuleFitTime(RnaSignal, TimeLabel, OutFolder, InheritSwitch)
    %%% Modified program to fit data across multiple time points considering delay effects

    % Data initialization and parameter setting
    % ... (Initialization part of your original code) ...

    % Initialize variables related to delay effects
    delayWindow = 3; % Assume a delay window of 3 time points
    KannealG = []; % Store fitting parameters for each time point
    PiniG = []; % Store initial state estimates for each time point

    % Main loop - process data for multiple time points simultaneously
    for T_i = 1:(TimeNum - delayWindow)
        % Retrieve data for the current time point and its delay window
        DataTimeWindow = DataCell(1, T_i:(T_i + delayWindow));

        % Call the new fitting function to process data across multiple time points
        [KFit, Pini] = FitMultipleTimePoints(DataTimeWindow, ...);

        % Store fitting results
        KannealG(T_i, :) = KFit;
        PiniG(T_i, :) = Pini;

        % ... Additional logic ...
    end

    % Output results
    % ... (Output part of your original code) ...
end

function [KFitResult, PiniResult] = FitMultipleTimePoints(DataTimeWindow, ...)
    % This function is for processing data across multiple time points
    % DataTimeWindow: Contains data for multiple time points

    % Logic for parameter estimation and model fitting
    % ... (Modify according to your specific model) ...

    % Return fitting parameters for multiple time points
    KFitResult = ...;
    PiniResult = ...;
end
