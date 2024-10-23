function [RnaSignal,ErrorIndex]=RnaSignalExtract(EmPath,RnaIndex,EL)
%%% a process to extract rna signal from results.
PlotSwitch=0;%figure;
UseNucleus=0;

RnaSignal=cell(size(EmPath,1),size(EL,2));
ErrorIndex=zeros(size(EmPath,1),1);
for Ei=1:size(EmPath,1)
    fname=EmPath{Ei};
    try
        clear foci_RNA_profile
        load(fname,'foci_RNA_profile')
        foci_RNA=foci_RNA_profile;
    catch
        fname=strrep(fname,'\Results\','\Results_new\');
    end
    try
        switch RnaIndex(Ei)% which rna signal to extract
                case 1
                    load(fname,'nucleus_RNA_profile','foci_RNA_profile')
                    nucleus_RNA=nucleus_RNA_profile;
                    foci_RNA=foci_RNA_profile;
                    clear nucleus_RNA_profile foci_RNA_profile
                case 2
                    load(fname,'nucleus_signal2_profile','foci_signal2_profile')
                    nucleus_RNA=nucleus_signal2_profile;
                    foci_RNA=foci_signal2_profile;
                    clear nucleus_signal2_profile foci_signal2_profile
                case 3
                    load(fname,'nucleus_signal4_profile','foci_signal4_profile')
                    nucleus_RNA=nucleus_signal4_profile;
                    foci_RNA=foci_signal4_profile;
                    clear nucleus_signal4_profile foci_signal4_profile
        end
    catch
        ErrorIndex(Ei)=1;
        disp(['Error: ',fname,' Channel ',num2str(RnaIndex(Ei))])
        continue
    end
    %% Data extract
    ind_foci = foci_RNA(:,2);
    I_foci = foci_RNA(:,3);
    N_nu = size(nucleus_RNA,1);
    [n_ind,i_ind] = hist(ind_foci,1:N_nu);
    repet_index=find(n_ind>4);
    repet_all=[];
    for i=1:length(repet_index)
        repet_foci=find(ind_foci==repet_index(i));
        repet_all=[repet_all;repet_foci];
    end
    foci_RNA(repet_all,:)=[];
    [C,ia,ic] = unique(foci_RNA(:,2));a_counts = accumarray(ic,1);
    value_counts = [C, a_counts];
    % try
    %     [nucleus_RNA_profile0,ind_foci] = foci_info(nucleus_RNA,foci_RNA);
    % catch
        [nucleus_RNA_profile0,ind_foci] = foci_info_old(nucleus_RNA,foci_RNA);
    % end
    if PlotSwitch
        nexttile
        UncleusSignalNumber = accumarray(a_counts(2:end,1),1);
        UncleusSignalNumberRatio=UncleusSignalNumber/sum(UncleusSignalNumber,1);
        bar(UncleusSignalNumberRatio);ylim([0 0.8])
        xlabel('Nucleus Signal Number');ylabel('Frequence');%title('Alignment')
    end
    signal_num=length(foci_RNA);
    for iii = 1:size(EL,2)
        reg_I = (nucleus_RNA(:,1) >= EL(1,iii)) & (nucleus_RNA(:,1) <= EL(2,iii));
        if UseNucleus==1
            f_ob = nucleus_RNA(reg_I,4);
        else
            reg_I0 = reg_I(ind_foci);
            pI0 = reg_I0;
            f_ob = nucleus_RNA_profile0(pI0,4);
        end
        % histogram(f_ob)
        RnaSignal{Ei,iii} = f_ob;
    end
end