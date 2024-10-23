function [PiniStart, LabelWait, Tfit, DtdInherit,timeBoundaries, featureValues,Kindex] = initializeParametersAndInherit(T_i,G_i,TimeElong, TimeLabel, PiniG, KannealG,KannealTime, InheritSwitch,InheritModel, Label, Tstay, Kstart, UnLockIndex,PlotSwitch)
    % Initialize parameters and handle inheritance for RNA distribution
    
    if T_i == 1
        PiniStart = [1 0 0]';  % Initial state probabilities for the first time point
        LabelWait = Label;     % Waiting label for RNA processing
        Tfit = Tstay;          % Fitting time
        DtdInherit=[];
        timeBoundaries=[];featureValues=[];Kindex=0;
    else
        % Calculate time elongation and determine evolution time interval
        [timeBoundaries, featureValues] = createTimelines(KannealTime(1:G_i-1,:), KannealG(1:G_i-1,:));
        % Test Assist Plot
        if PlotSwitch
            plotTimeWindows(KannealTime(1:G_i-1,:), KannealG(1:G_i-1,:),TimeElong);
        end
        % Determine the indices for evolution
        switch InheritModel
            case 'Near'
                TimeNum=size(TimeLabel,2); % Number of time points
                Index=(1:TimeNum);
                EvoKsIndex = TimeLabel > TimeElong(1) & TimeLabel < TimeElong(2);
                % EvoStartK = find(TimeLabel < TimeElong(1), 1, 'last') + 1;
                Kindex=Index(EvoKsIndex);
                EvoStartP=Index(max(TimeLabel<TimeElong(1))+1);
                EvoStartK=EvoStartP;
                TimeNode=sort([TimeLabel(EvoKsIndex),TimeElong]);
                TimeLabelE=TimeLabel;
                kUse=KannealG;
                KannealTimePini=TimeLabel;
            case 'Label'
                TimeLabelE=timeBoundaries(2:end);
                EvoKsIndex = TimeLabelE > TimeElong(1) & TimeLabelE < TimeElong(2);
                TimeNum=size(TimeLabelE,2); % Number of time points
                Index=(1:TimeNum+1);
                % EvoStartK = find(TimeLabel < TimeElong(1), 1, 'last') + 1;
                Kindex=Index(EvoKsIndex);
                EvoStartP=max([Index(TimeLabel<TimeElong(1)) 0]);
                
                TimeNode=sort([TimeLabelE(EvoKsIndex),TimeElong]);
                
                kUse=featureValues;
                KUndefinedIndex=find(kUse(:,1)==0);
                if ~isempty(KUndefinedIndex)
                    for iii=KUndefinedIndex'
                        try
                            kUse(iii,:)=mean(kUse([iii-1 iii+1],:),1);
                        catch
                            kUse(iii,:)=mean(kUse(iii-1,:),1);
                        end
                    end
                end
                KannealTimePini=[TimeLabel(1)-60 TimeLabel];
                EvoStartK=(Index(timeBoundaries<TimeElong(1)& timeBoundaries>=KannealTimePini(EvoStartP+1)));
                if isempty(EvoStartK)
                    try
                        EvoStartK=Index(timeBoundaries>=KannealTimePini(EvoStartP+1));
                        EvoStartK=EvoStartK(1);
                    catch
                        EvoStartK=length(timeBoundaries);
                    end
                end
        end
        % Initialize PiniStart based on the previous time points
        PiniStart = calculatePiniStart(EvoStartP,EvoStartK, TimeLabel,timeBoundaries,KannealTimePini, PiniG, kUse, Kstart, UnLockIndex,TimeElong,PlotSwitch);

        % Handle RNA distribution inheritance
        DtdInherit=[];
        if InheritSwitch == 1
            [LabelWait, Tfit, DtdInherit] = handleInheritance(EvoKsIndex, TimeLabelE, TimeElong, kUse, Kstart,UnLockIndex,PiniStart, Label, Kindex, TimeNode, DtdInherit,PlotSwitch);
        else
            LabelWait = Label;
            Tfit = TimeElong(end) - TimeElong(1);
        end
        % close gcf
    end
end

function PiniStart = calculatePiniStart(EvoStartP,EvoStartK, TimeLabel,timeBoundaries,KannealTimePini, PiniG, KannealG, Kstart, UnLockIndex, TimeElong,PlotSwitch)
    % Calculate initial state probabilities (PiniStart) for the current time point
    Color=['b','r'];ColorI=1;
    if EvoStartP == 0
        PiniStart0 = [1 0 0]';  % Default initial state for the first time point
        TimeStart0=TimeLabel(1)-60;
    else
        % Use the previous time point's data for initialization
        PiniStart0 = PiniG(EvoStartP, :)';
        TimeStart0=KannealTimePini(EvoStartP+1);
        % Additional logic for calculating PiniStart based on previous time points
        % This may involve evolving the state based on KannealG and other factors
    end
    if PlotSwitch
        text(TimeStart0,0,{'State0','\diamondsuit'},'color','r','HorizontalAlignment','center','FontSize',12,'FontWeigh','bold')%
    end
    KannealG0=[zeros(1,size(KannealG,2));KannealG];
    for Ki=EvoStartK
        try
            Tk=timeBoundaries(Ki+1);%deal ini evo nan k
            Kinherit=Ki;
        catch
            Tk=TimeElong(1);
            Kinherit=Ki-1;
        end
        TevoS=max([min([TimeElong(1) Tk])-TimeStart0,0]);
        if PlotSwitch
            line([min([TimeElong(1) Tk]) TimeStart0],[0 0],'linewidth',4,'color',Color(ColorI))%
            ColorI=3-ColorI;
        end
        [~,~,~,PiniStart,~] =  FSP_FreeSet(KannealG(Kinherit,:),KannealG0(Kinherit,:),0,PiniStart0,[0 0 1],0,Kstart,UnLockIndex,TevoS,[],0);%%MaybeErrorExist
        PiniStart0=PiniStart;
        TimeStart0=Tk;
    end
end

function [LabelWait, Tfit, DtdInherit] = handleInheritance(EvoKsIndex, TimeLabelE, TimeElong, KannealG,Kstart,UnLockIndex, PiniStart, Label, Kindex, TimeNode, DtdInherit,PlotSwitch)
    % Handle RNA distribution inheritance based on previous time points
    LabelWait=Label;
    Color=['b','r'];ColorI=1;
    KannealG0=[zeros(1,size(KannealG,2));KannealG];
    %Inherit of Rna distribution
        EvoI=0;
        LabelUse=[];
        for K_i=Kindex 
            EvoI=EvoI+1;
            TimeAtK=(TimeLabelE(K_i)-TimeNode(EvoI));
            if PlotSwitch
                line([TimeLabelE(K_i) TimeNode(EvoI)],[0 0],'linewidth',8,'color',Color(ColorI))%
                ColorI=3-ColorI;
            end
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
            [~,~,~,~,DtdInherit] =  FSP_FreeSet(KannealG(K_i,:),KannealG0(K_i,:),0,PiniStart,LabelPart,0,Kstart,UnLockIndex,TimeAtKev,DtdInherit,0);
        end
        Tfit=TimeElong(end)-TimeNode(end-1);
    
    % Additional logic to evolve the RNA distribution based on inheritance
    % This may involve updating LabelWait and Tfit based on KannealG and other factors
end