function ModuleFit(in0f,out_folder_use,data_cycle,el,protein_concentration_all,on_off,simu,model,enrich)
%% A program to fit the experiment data by model

%% Initialize required parameters

% k=[k01 k02 k10 k12 k20 k21 ki0 ki1 ki1 useless dataintensity alpha]
k_scan=0;% scan the parameters' space
gaussfit=0;% guass fit
poiss_width=6;% the width of peak to calulate the enrichment
intr_peak=0;
EL=el;
transform = '0112';% In 3-state model: 0-1-2 transforming:0112 or 0221 
%% divide the data with cycle to analyze
for i_cycle = [1]
    %% consider [all cycle\cycle11\12\13\14]
    out_folder=out_folder_use;
    data_all = data_cycle{i_cycle};
    fig_title='allcycles';
    protein_concentration=protein_concentration_all{1};
    hb1k10=3;hb7k10=1;
    hb1k21=1.202;hb7k21=0.15;
    switch i_cycle
        case 2
            out_folder=[out_folder,'\','cycle11'];
            mkdir(out_folder);
            fig_title='cycle11';
            hb1k10=3;hb7k10=1;
            hb1k21=1.202;hb7k21=0.15;
            protein_concentration=protein_concentration_all{2};
        case 3
            out_folder=[out_folder,'\','cycle12'];
            mkdir(out_folder);
            fig_title='cycle12';
            hb1k10=0.3548;hb7k10=0.8414;
            hb1k21=0.6026;hb7k21=0.79433;
            protein_concentration=protein_concentration_all{3};
        case 4
            out_folder=[out_folder,'\','cycle13'];
            mkdir(out_folder);
            fig_title='cycle13';
            hb1k10=0.0631;hb7k10=0.1496;
            hb1k21=0.3981;hb7k21=1.8197;
            protein_concentration=protein_concentration_all{4};
        case 5
            out_folder=[out_folder,'\','cycle14'];
            mkdir(out_folder);
            fig_title='cycle14';
    end  

    %% choose model to fit the data(model 2\3 state)
    switch model
        case 2
            ksim_group1=zeros(length(data_all),5);
            ksim_group2=zeros(length(data_all),5);
            pini_group1=zeros(length(data_all),2);
            pini_group2=zeros(length(data_all),2);
            lb=0*ones(1,5);
            ub=500*ones(1,5);
            ub(3)=0;lb([1,2])=0.1;ub(5)=1;
            ub([1,2])=100;
        case 3
            ksim_group1=zeros(length(data_all),12);
            ksim_group2=zeros(length(data_all),12);
            pini_group1=zeros(length(data_all),3);
            pini_group2=zeros(length(data_all),3);
            lb=[zeros(1,10),1];
            lb([1:6])=0;
            lb([3,6])=[0.2,0.5];
            lb([1,4])=[0.001,0.0001];
            ub=[1*ones(1,6),1000*ones(1,3),0,1];
            ub([3,6])=[0.3,0.7];
            ub([1,4])=[0.1,0.3];
%             ub([2,5])=[20,20];
            if strcmp(transform,'0112')
                ub([2,5,7])=0;lb([2,5,7])=0;
            else
                ub ([4,6,7])=0;lb([4,6,7])=0;
            end
            lb_2=lb;
            lb_1=lb;
            ub_1=ub;
%             lb_1([3,6])=[0.2,0.5];
            lb_1([3,6])=[0.1,0.1];
            lb_1([1,4])=[0.001,0.0001];
            lb_1([8,9])=[1,8];
%             ub_1([3,6])=[0.5,0.7];
            ub_1([3,6])=[5,5];
%             ub_1([1,4])=[0.1,0.3];
            ub_1([1,4])=[5,5];
            ub_1([8,9])=[5,18];
            ub_2=ub;
%             lb_2([3,6])=[0.8,0.75];
            lb_2([3,6])=[0.1,0.1];
            lb_2([1,4])=[0.001,0.0001];
            lb_2([8,9])=[45,85];
%             ub_2([3,6])=[0.9,0.85];
            ub_2([3,6])=[5,5];
            ub_2([1,4])=[5,2];
            ub_2([8,9])=[55,95];
%             lb_1(9)=13.5;ub_1(9)=13.5;lb_1(8)=4;ub_1(8)=4;
%             ub_1(8)=0;ub_2(8)=0;
%             lb_1(8)=0;lb_2(8)=0;
%             ub_1(9)=15;ub_2(9)=30;
%             lb_1(9)=0;lb_2(9)=0;
    end
    alpha_double_up=1;alpha_double_dn=0;alpha_double=0.5;
    %% Initialize required parameters
    [lld1_cell,lld2_cell,lld_scan,hb1_scan_g,hb7_scan_g,hb1_scan_k,hb7_scan_k]=deal...
        (cell(length(data_all),1),cell(length(data_all),1),zeros(length(data_all),2),zeros(length(data_all),2),zeros(length(data_all),2),zeros(length(data_all),4),zeros(length(data_all),4));
    lambda_group=zeros(length(data_all),4);
    guass_group=zeros(length(data_all),4);
    lld_group=zeros(length(data_all),2);
    n_group=zeros(length(data_all),8);
    hb_all=zeros(length(data_all),2);
    enrich_peak_gr=zeros(length(data_all),4);
%     enrich_peak_origin_gr=cell(length(data_all),4);
    nb=0:100;
    T = [0.0037,0.0722,0.4211,0.48,0.0096,0.0619,0.0071];
    T1=sum(T(1:4));
    T2=sum(T(end-3:end));
    Nbin1=2000;Nbin2=Nbin1*T2/T1;
    MeanPlot(mean(el)',data_all','')
    legend('hb1','hb7')
    saveas(gcf,[[out_folder,'\'],'Mean_el','.fig']);
    saveas(gcf,[[out_folder,'\'],'Mean_el','.png']);
    close(gcf);
    
    MeanPlot(mean(protein_concentration)',data_all','')
    legend('hb1','hb7')
    saveas(gcf,[[out_folder,'\'],'Mean_protein','.fig']);
    saveas(gcf,[[out_folder,'\'],'Mean_protein','.png']);
    close(gcf);
    
    MeanPlot(mean(el)',data_all','nozero')
    legend('hb1','hb7')
    saveas(gcf,[[out_folder,'\'],'Mean_nozero','.fig']);
    saveas(gcf,[[out_folder,'\'],'Mean_nozero','.png']);
    close(gcf);
    
    figure
    plot(mean(el)',mean(protein_concentration)')
    saveas(gcf,[[out_folder,'\'],'EL_protein','.fig']);
    saveas(gcf,[[out_folder,'\'],'EL_protein','.png']);
    close(gcf);
    %% parallel fitting of different el ranges
    PoolCreate(length(data_all));
%     parpool(length(data_all))
    parfor fig_order=1:length(data_all)
        LockList=[0 0 0 0 0 0 0 0 0 0 0 0];
        LockNum1=[0 0 hb1k10 0 0 hb1k21 0 0 0 0 0 0];
        LockNum2=[0 0 hb7k10 0 0 hb7k21 0 0 0 0 0 0];
        %% data extraction
        if isempty(data_all{fig_order})
            data_all{fig_order}=zeros(1,4);
        end
        f_ob=data_all{fig_order}(:,1);
        f_ob_nostandard=data_all{fig_order}(:,1);
        f_ob(f_ob>100)=100;
        hb1_all=sum(f_ob);
        f_ob_2=data_all{fig_order}(:,2);
        f_ob_2(f_ob_2>100)=100;
        f_ob_2_nz=f_ob_2(f_ob_2>0);
        f_ob_2_fano=var(f_ob_2_nz)/mean(f_ob_2_nz);
%         [mean(f_ob_2),var(f_ob_2),var(f_ob_2)/(mean(f_ob_2).^2)]
        hb7_all=sum(f_ob_2);
        hb_all(fig_order,:)=[hb1_all,hb7_all];
        Nbin_data=0:1:100;
        nf = histc(reshape(f_ob,1,length(f_ob)),Nbin_data);nb=Nbin_data;
        nf_2 = histc(reshape(f_ob_2,1,length(f_ob_2)),Nbin_data);
        nf_all=[nf/sum(nf);nf_2/sum(nf_2)]';
    %             state0=nf_all(1,:);
        %% plot the original data
        bar(nb,nf_all,'hist');
        ylim([0 0.08])
        xlim([0 60])
        saveas(gcf,[[out_folder,'\'],['hb17_',num2str(fig_order),'.fig']]);
        saveas(gcf,[[out_folder,'\'],['hb17_',num2str(fig_order),'.png']]);
        close(gcf);
        
        subplot(1,2,1)    
            bar(nb,nf'/sum(nf),'y');
            ylim([0 0.08])
            xlim([0 60])
            xlabel('# hb1')
            ylabel('Frequency')
        subplot(1,2,2)
            bar(nb,nf_2'/sum(nf_2),'g');
            ylim([0 0.08])
            xlim([0 60])
            xlabel('# hb7')
            ylabel('Frequency') 
            text(48,0.07,{num2str(f_ob_2_fano)})
        saveas(gcf,[[out_folder,'\'],['hb17_div_',num2str(fig_order),'.fig']]);
        saveas(gcf,[[out_folder,'\'],['hb17_div_',num2str(fig_order),'.png']]);
        close(gcf);
        
        Nbin_data_2=0:2:100;
        nf_bin2 = histc(reshape(f_ob_2,1,length(f_ob_2)),Nbin_data_2);
        nf_bin2=nf_bin2'/sum(nf_bin2);
        bar(Nbin_data_2,nf_bin2,'hist');
        ylim([0 0.14])
        xlim([0 60])
        xlabel('# RNA2')
        ylabel('Frequency') 
        text(48,0.1,{num2str(f_ob_2_fano)})
        saveas(gcf,[[out_folder,'\'],['RNA2_BIN2_',num2str(fig_order),'.fig']]);
        saveas(gcf,[[out_folder,'\'],['RNA2_BIN2_',num2str(fig_order),'.png']]);
        close(gcf);
        %% Trying to find cyclicity from data  
        if intr_peak==1
            f_ob_nostandard=f_ob_nostandard(f_ob_nostandard<300);
            for intr_bin=0.1
                for i0 = 2:0.1:10
                    f_ob_nostandard_use=f_ob_nostandard/i0;
                    intron_bin=intr_bin;
                    Nbin_data=0:intron_bin:300;nb=Nbin_data;
                    nf_2 = histc(reshape(f_ob_nostandard_use,1,length(f_ob_nostandard_use)),Nbin_data);
                    nf_2(1)=0;
                    nf_2=[nf_2/sum(nf_2)]';
                    nfnz_index=find(nf_2~=0);
%                     nf_2_nz=nf_2(nfnz_index(1):end);
                    in_cycle=5+round(nfnz_index(1)/10);

                    i0_binx=1/intr_bin;
                    nf_2_nz_i00=nf_2(1:in_cycle*i0_binx);
                    nf_2_nz_i0=reshape(nf_2_nz_i00,[i0_binx,in_cycle]);
                    i0_stack=sum(nf_2_nz_i0,2);
                    plot(intr_bin:intr_bin:1,i0_stack,'-ks','LineWidth',1)
                    hold on
                    plot([0.55,0.55],ylim,'m--');
                    xlabel('I_0')
                    ylabel('Pre')
                    title([fig_title,'-- El:',num2str(el(1,fig_order)),'~',num2str(el(2,fig_order)),'bin',num2str(intron_bin),'I0',num2str(i0)])
                    saveas(gcf,[[out_folder,'\'],['intronnz',num2str(fig_order),'bin',num2str(intr_bin),'棗i0',num2str(i0),'.fig']]);
                    saveas(gcf,[[out_folder,'\'],['intronnz',num2str(fig_order),'bin',num2str(intr_bin),'棗i0',num2str(i0),'.png']]);
                    close(gcf);
               end
            end
        end
        %% analysis start
        if on_off==1
            %% data extraction
            f_ob=data_all{fig_order}(:,1);
            f_ob(f_ob>100)=100;
            f_ob_2=data_all{fig_order}(:,2);
            f_ob_2(f_ob_2>100)=100;
            Nbin_data=0:1:100;
            nf = histc(reshape(f_ob,1,length(f_ob)),Nbin_data);nb=Nbin_data;
            nf_2 = histc(reshape(f_ob_2,1,length(f_ob_2)),Nbin_data);
            nf_all=[nf/sum(nf);nf_2/sum(nf_2)]';
            nf_guass=[[0,0];nf_all(2:end,:)];
            %% gauss fit
            if gaussfit==1
                %% gauss2 fit
                fit_1=fit(Nbin_data(1:end)',nf_guass(1:end,1),'gauss2');
                var_1=[(fit_1.c1^2/2),(fit_1.c2^2/2)];
                mean_1=[fit_1.b1,fit_1.b2];
                guass_intensity_1=[fit_1.a1,fit_1.a2];
                lambda_1=var_1./mean_1;

                subplot(1,2,1)
                    scatter(Nbin_data(1:end),nf_guass(1:end,1),15,'filled','k')
                    hold on
                    plot(fit_1,'r')
                    xlim([0 50])
                    xlabel('# rna1')
                    ylabel('# frequence')
                    title('Gauss2 Fitting')
                    text(30,0.8*max(ylim),{['\lambda_{peak1} = ',num2str(lambda_1(1))],['\lambda_{peak2} = ',num2str(lambda_1(2))]})
                fit_2=fit(Nbin_data(1:end)',nf_guass(1:end,2),'gauss2');
                var_2=[(fit_2.c1^2/2),(fit_2.c2^2/2)];
                mean_2=[fit_2.b1,fit_2.b2];
                lambda_2=var_2./mean_2;
                lambda_group(fig_order,:)=[lambda_1,lambda_2];
                guass_intensity_2=[fit_2.a1,fit_2.a2];
                guass_group(fig_order,:)=[guass_intensity_1,guass_intensity_2];
                subplot(1,2,2)
                    scatter(Nbin_data(1:end),nf_guass(1:end,2),15,'filled','k')
                    hold on
                    plot(fit_2,'g')
                    xlim([0 50])
                    xlabel('# rna2')
                    ylabel('# frequence')
                    title('Gauss2 Fitting')
                    text(30,0.8*max(ylim),{['\lambda_{peak1} = ',num2str(lambda_2(1))],['\lambda_{peak2} = ',num2str(lambda_2(2))]})
                saveas(gcf,[[out_folder,'\'],'Guass2_',num2str(fig_order),'.fig']);
                saveas(gcf,[[out_folder,'\'],'Guass2_',num2str(fig_order),'.png']);
                close(gcf);
                %% gauss1 fit
                fit_1=fit(Nbin_data(1:end)',nf_guass(1:end,1),'gauss1');
                var_1=(fit_1.c1^2/2);
                mean_1=fit_1.b1;
                lambda_1=var_1./mean_1;
                subplot(1,2,1)
                    scatter(Nbin_data(1:end),nf_guass(1:end,1),15,'filled','k')
                    hold on
                    plot(fit_1,'r')
                    xlim([0 50])
                    xlabel('# rna1')
                    ylabel('# frequence')
                    title('Gauss Fitting')
                    text(30,0.8*max(ylim),{['var = ',num2str(var_1)],['mean = ',num2str(mean_1)]})
                fit_2=fit(Nbin_data(1:end)',nf_guass(1:end,2),'gauss1');
                var_2=(fit_2.c1^2/2);
                mean_2=fit_2.b1;
                lambda_2=var_2./mean_2;
    %             lambda_group(fig_order,:)=[lambda_1,lambda_2];
                subplot(1,2,2)
                    scatter(Nbin_data(1:end),nf_guass(1:end,2),15,'filled','k')
                    hold on
                    plot(fit_2,'g')
                    xlim([0 50])
                    xlabel('# rna2')
                    ylabel('# frequence')
                    title('Gauss Fitting')
                    text(30,0.8*max(ylim),{['var = ',num2str(var_2)],['mean = ',num2str(mean_2)]})
                saveas(gcf,[[out_folder,'\'],'Guass_',num2str(fig_order),'.fig']);
                saveas(gcf,[[out_folder,'\'],'Guass_',num2str(fig_order),'.png']);
                close(gcf);
            end
            %% poission fit
            nf = histc(reshape(f_ob,1,length(f_ob)),Nbin_data);nb=Nbin_data;
            nf_2 = histc(reshape(f_ob_2,1,length(f_ob_2)),Nbin_data);
            nf_all=[nf/sum(nf);nf_2/sum(nf_2)]';
            state0=nf_all(1,:);
            bar(nb,nf_all,'hist');
            ylim([0 0.08])
            xlim([0 60])
            hold on
            % double poisson %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [para1,poiss2_dtb_best1,para2,poiss2_dtb_best2] = PoissonFitScan(f_ob,state0(1),f_ob_2,state0(2));
            poiss2_dtb_best1(1)=state0(1);
            poiss2_dtb_best2(1)=state0(2);
            plot(0:60,poiss2_dtb_best1);
            plot(0:60,poiss2_dtb_best2);
            legend('RNA1','RNA2','poiss1','poiss2')
            saveas(gcf,[[out_folder,'\'],['dualpioss',num2str(fig_order),'.fig']]);
            saveas(gcf,[[out_folder,'\'],['dualpioss',num2str(fig_order),'.png']]);
            close(gcf);
%             para_fig=[para1([2,1])',para2([2,1])'];
%             para_fig=cat(1,1-[sum(para1(1:2)),sum(para2(1:2))],para_fig);
%             x_fig=categorical({'Ground State','Strong Activation','Weak Activation'});
%             bar(x_fig,para_fig);
%             legend('RNA 1','RNA 2')
%             xlabel('State')
%             ylabel('Proportion')
%             saveas(gcf,[[out_folder,'\'],'Propor_state',num2str(fig_order),'.fig']);
%             saveas(gcf,[[out_folder,'\'],'Propor_state',num2str(fig_order),'.png']);
%             close(gcf);
            %% enrichment calulate
            if enrich == 1
                peak11 = f_ob(:,1)<=para1(3) & f_ob(:,1)>=para1(3)-2*poiss_width &f_ob(:,1)>0;
                peak_enrich11_or = data_all{fig_order}(peak11,3);            
                peak_enrich11_or = peak_enrich11_or(peak_enrich11_or~=0);
                peak_enrich11 = mean(peak_enrich11_or);
                peak12 = f_ob(:,1)<=para1(4)+2*poiss_width & f_ob(:,1)>=para1(4)-poiss_width;
                peak_enrich12_or = data_all{fig_order}(peak12,3);
                peak_enrich12_or = peak_enrich12_or(peak_enrich12_or~=0);
                peak_enrich12 = mean(peak_enrich12_or);

                peak21 = f_ob_2(:,1)<=para2(3) & f_ob_2(:,1)>=para2(3)-2*poiss_width & f_ob(:,1)>0;
                peak_enrich21_or = data_all{fig_order}(peak21,4);
                peak_enrich21_or = peak_enrich21_or(peak_enrich21_or~=0);
                peak_enrich21 = mean(peak_enrich21_or);
                peak22 = f_ob_2(:,1)<=para2(4)+2*poiss_width & f_ob_2(:,1)>=para2(4)-poiss_width;
                peak_enrich22_or = data_all{fig_order}(peak22,4);
                peak_enrich22_or = peak_enrich22_or(peak_enrich22_or~=0);
                peak_enrich22 = mean(peak_enrich22_or);
                enrich_peak=[[peak_enrich12,peak_enrich11];[peak_enrich22,peak_enrich21]];
                enrich_peak_gr(fig_order,:)=[peak_enrich11,peak_enrich12,peak_enrich21,peak_enrich22];
%                 enrich_fig=categorical({'Strong Activation','Weak Activation'});
%                 bar(enrich_fig,enrich_peak);
%                 legend('RNA 1','RNA 2')
%                 xlabel('State')
%                 ylabel('Enrichment Value')
%                 saveas(gcf,[[out_folder,'\'],'Enrichment_rna',num2str(fig_order),'.fig']);
%                 saveas(gcf,[[out_folder,'\'],'Enrichment_rna',num2str(fig_order),'.png']);
%                 close(gcf);
% 
%                 enrich_fig2=categorical({'RNA 1','RNA 2'});
%                 bar(enrich_fig2,enrich_peak');
%                 legend('Strong Activation','Weak Activation')
%                 xlabel('RNA')
%                 ylabel('Enrichment Value')
%                 saveas(gcf,[[out_folder,'\'],'Enrichment_state',num2str(fig_order),'.fig']);
%                 saveas(gcf,[[out_folder,'\'],'Enrichment_state',num2str(fig_order),'.png']);
%                 close(gcf);
            end
            %% model2\3 fit
            if simu == 1
                %% repeat the simulannealbnd to find the best one
                lld1_max=inf;lld2_max=inf;
                turn=0;
                while turn < 3  % repeat turn
                    turn = turn + 1;
                    %% pick the model
                    switch model
                        case 2
                        %% calulate the initial k by poission fit
                            k_new=rand(1,20);
                            k=k_new;
                            n11=para1(1);
                            n10=para1(2);
                            n01=para2(1);
                            n00=state0(2)-n10;
                            n1=2*(0.5*73.9+431+491.3);n2=2*(0.5*110.8+12.7);
                            k_new=[0.1*ones(1,6),1000*rand(1,3),0,1];
                            k_new([2,5,7])=0;
                            k_new1=k_new;
                            k_new2=k_new;
                            k_new1(3)=k_new1(1)*(1-(para1(1)+para1(2)))/para1(1);
                            k_new1(6)=k_new1(4)*para1(1)/para1(2);
                            k_new2(3)=k_new2(1)*(1-(para2(1)+para2(2)))/para2(1);
                            k_new2(6)=k_new2(4)*para2(1)/para2(2);
                            k_new1([9,8])=[Nbin1*para1(4)/n1,Nbin1*para1(3)/n1];
                            k_new2([8,9])=[Nbin1*para2(3)/n2,Nbin1*para2(4)/n2];

                            k_new1=[0.01*(para1(1)+para1(2)),0.01*(1-(para1(1)+para1(2))),0,Nbin1*para1(3)/n1,0.3];
                            k_new2=[0.01*(para2(1)+para2(2)),0.01*(1-(para2(1)+para2(2))),0,Nbin1*para2(3)/n2,0.3];

                            p_ini=[0.5,0.3,0.2]';
                        %     p_ini=[n00,n01,n10,n11]';
                        %% simulannealbnd
                            [kk_simu1,~,~] = simulannealbnd(@(x) FSP_NSX_pc(x,2,f_ob,sum(T(3:4))/sum([T(1:2),T(3)+T(4)]),T(1)/sum([T(1:2),T(3)+T(4)])),k_new1,lb,ub);
                            [kk_simu2,~,~] = simulannealbnd(@(x) FSP_NSX_pc(x,2,f_ob_2,T(end)/sum([T(end-3)+T(end-2),T(end-1:end)]),sum(T(end-3:end-2))/sum([T(end-3)+T(end-2),T(end-1:end)])),k_new2,lb,ub);
                            ksim_group1(fig_order,:)=kk_simu1;
                            ksim_group2(fig_order,:)=kk_simu2;
%                             [~,P_dtb1] =  FSP_NSX_pc(kk_simu1,2,f_ob,sum(T(3:4))/sum([T(1:2),T(3)+T(4)]),T(1)/sum([T(1:2),T(3)+T(4)]));
%                             [~,P_dtb2] =  FSP_NSX_pc(kk_simu2,2,f_ob_2,T(end)/sum([T(end-3)+T(end-2),T(end-1:end)]),sum(T(end-3:end-2))/sum([T(end-3)+T(end-2),T(end-1:end)]));
%                             pini_group1(fig_order,:)=[(1-(para1(1)+para1(2))),(para1(1)+para1(2))];
%                             pini_group2(fig_order,:)=[(1-(para2(1)+para2(2))),(para2(1)+para2(2))];
%                             P_dtb1=P_dtb1(:,2);
%                             P_dtb2=P_dtb2(:,2);
                        case 3
                        %% calulate the initial k by poission fit
                            k_new=rand(1,20);
                            k=k_new;
                            n11=para1(1);
                            n10=para1(2);
                            n01=para2(1);
                            n00=state0(2)-n10;
                            n1=2*(0.5*73.9+431+491.3);n2=2*(0.5*110.8+12.7);
                        %% model 3 -0112 or 0221
                            if strcmp(transform,'0112')
                                k_new=[0.01*ones(1,6),1000*rand(1,3),0,1];
                                k_new([2,5,7])=0;
                                k_new1=k_new;
                                k_new2=k_new;
                                k_new1(3)=k_new1(1)*(1-(para1(1)+para1(2)))/para1(1);
                                k_new1(6)=k_new1(4)*para1(1)/para1(2);
                                k_new2(3)=k_new2(1)*(1-(para2(1)+para2(2)))/para2(1);
                                k_new2(6)=k_new2(4)*para2(1)/para2(2);
%                                 hold_time_1=0.4753;hold_time=0.4218;
                                hold_time_1=0.4753;hold_time=0.1528;
                                Nbin1_2=Nbin1+Nbin1*hold_time;
                                n1=2*(0.5*73.9+431+491.3+Nbin1_2/2*hold_time);n2=2*(0.5*110.8+12.7+Nbin1_2/2*hold_time);
                                Nbin1_1=2000+2000*hold_time_1;
                                n1_1=2*(0.5*73.9+431+491.3+Nbin1_1/2*hold_time_1);n2_1=2*(0.5*110.8+12.7+Nbin1_1/2*hold_time_1);
                                k_new1([9,8])=[Nbin1*para1(4)/n1_1,Nbin1*para1(3)/n1_1];
                                k_new2([8,9])=[Nbin1*para2(3)/n2,Nbin1*para2(4)/n2];
                            else
                                k_new=[0.01*ones(1,6),1000*rand(1,3),0,1];
    %                             k_new([2,5])=1*ones(1,2);
                                k_new([4,6,7])=0;
                                k_new1=k_new;
                                k_new2=k_new;
                                k_new1(3)=k_new1(1)*(1-(para1(1)+para1(2)))/para1(1);
                                k_new1(5)=k_new1(2)*(1-(para1(1)+para1(2)))/para1(2);
                                k_new2(3)=k_new2(1)*(1-(para2(1)+para2(2)))/para2(1);
                                k_new2(5)=k_new2(2)*(1-(para2(1)+para2(2)))/para2(2);
    %                             hold_time=0.94;hold_time_1=0.39;
    %                             hold_time=0.1224;hold_time_1=0.5375;
                                hold_time_1=0.4753;hold_time=0.4218;
                                Nbin1_2=Nbin1+Nbin1*hold_time;
                                n1=2*(0.5*73.9+431+491.3+Nbin1_2/2*hold_time);n2=2*(0.5*110.8+12.7+Nbin1_2/2*hold_time);
                                Nbin1_1=2000+2000*hold_time_1;
                                n1_1=2*(0.5*73.9+431+491.3+Nbin1_1/2*hold_time_1);n2_1=2*(0.5*110.8+12.7+Nbin1_1/2*hold_time_1);
                                k_new1([8,9])=[Nbin1_1*para1(3)/n1_1,Nbin1_1*para1(4)/n1_1];
                                % kk_simu1(8)=Nbin1_1*(para11+5)/n1_1;
                                k_new2([8,9])=[Nbin1_2*para2(3)/n2,Nbin1_2*para2(4)/n2];          
    %                             k_new1([9,8])=[Nbin1*para1(4)/n1,Nbin1*para1(3)/n1];
    %                             k_new2([8,9])=[Nbin1*para2(3)/n2,Nbin1*para2(4)/n2];
                            end
                            p_ini=[0.5,0.3,0.2]';
                        %     p_ini=[n00,n01,n10,n11]';
                            errorornot1=1;errorornot2=1;
                        %% double copy ratio
                            
                        %% simulannealbnd and scan
                            if strcmp(transform,'0112')
                                while errorornot1
%                                     k_new1([1 3])=k_new1([1 3])*5;
                                    k_new1(9)=13.5;k_new1(8)=3.5;
                                    LockList=logical(LockList);
                                    Start=[k_new1,alpha_double];
                                    Lb=[lb_1,alpha_double_dn];
                                    Ub=[ub_1,alpha_double_up];
                                    Start(LockList)=LockNum1(LockList);
                                    Lb(LockList)=LockNum1(LockList);
                                    Ub(LockList)=LockNum1(LockList);
%                                      k_new1(8)=0;k_new1(9)=4;
%                                         k_new1(8:9)=k_new1([9,8]);
                                    [kk_simu1_double,~,~] = simulannealbnd(@(x) FSP_threestate(x,f_ob,p_ini,[T(1:2),T(3)+T(4)],hold_time_1),...
                                        Start,Lb,Ub);
%                                     kk_simu1_double=[k_new1,alpha_double];kk_simu1_best=[k_new1,alpha_double];
                                    kk_simu1=kk_simu1_double;
                                    if kk_simu1(8)<kk_simu1(9)
                                        errorornot1=0;          
                                    end
                                end
                                while errorornot2
                                    k_new2(9)=60;k_new2(8)=25;
                                    LockList=logical(LockList);
                                    Start=[k_new2,alpha_double];
                                    Lb=[lb_2,alpha_double_dn];
                                    Ub=[ub_2,alpha_double_up];
                                    Start(LockList)=LockNum2(LockList);
                                    Lb(LockList)=LockNum2(LockList);
                                    Ub(LockList)=LockNum2(LockList);
%                                         k_new2(8)=0;k_new2(9)=13;
%                                         k_new2(8:9)=k_new2([9,8]);
%                                     [kk_simu2_double,~,~] = simulannealbnd(@(x) FSP_threestate(x,f_ob_2,p_ini,[T(end-3)+T(end-2),T(end-1:end)],hold_time),...
%                                         [k_new2,alpha_double],[lb_2,alpha_double_dn],[ub_2,alpha_double_up]);
%                                     T_cds_2=[500,4325,930+1145]/sum([500,4325,930+1145]); 
                                    T_cds_2=[T(end-3)+T(end-2),T(end-1:end)];
                                    [kk_simu2_double,~,~] = simulannealbnd(@(x) FSP_threestate(x,f_ob_2,p_ini,T_cds_2,hold_time),...
                                        Start,Lb,Ub);
%                                     kk_simu2_double=[k_new2,alpha_double];
                                    kk_simu2=kk_simu2_double;
                                    if kk_simu2(8)<kk_simu2(9)
                                        errorornot2=0;  
                                    end
                                end
                            else
                                [kk_simu1_double,~,~] = simulannealbnd(@(x) FSP_threestate(x,f_ob,p_ini,[T(1:2),T(3)+T(4)],hold_time_1),[k_new1,0.12],[lb,0.1],[ub_1,0.1]);
                                [kk_simu2_double,~,~] = simulannealbnd(@(x) FSP_threestate(x,f_ob_2,p_ini,[T(end-3)+T(end-2),T(end-1:end)],hold_time),[k_new2,0.1],[lb,0.1],[ub_2,0.1]);
                                if kk_simu1_double(9)<kk_simu1_double(8)
                                    kk_simu1_double([8,9])=kk_simu1_double([9,8]);
                                    kk_simu1_double([1,2,3,5])=kk_simu1_double([2,1,5,3]);
                                end
                                if kk_simu2_double(9)<kk_simu2_double(8)
                                    kk_simu2_double([8,9])=kk_simu2_double([9,8]);
                                    kk_simu2_double([1,2,3,5])=kk_simu2_double([2,1,5,3]);
                                end
                                kk_simu1=kk_simu1_double;
                                kk_simu2=kk_simu2_double;
                            end 
                            [lld1,~,~,p_ini1] =  FSP_threestate(kk_simu1_double,f_ob,p_ini,[T(1:2),T(3)+T(4)],hold_time_1);
                            [lld2,~,~,p_ini2] =  FSP_threestate(kk_simu2_double,f_ob_2,p_ini,[T(end-3)+T(end-2),T(end-1:end)],hold_time);
                            if lld1<lld1_max
                                ksim_group1(fig_order,:)=kk_simu1;
                                pini_group1(fig_order,:)=p_ini1;
                                lld1_max=lld1;
                                kk_simu1_best=kk_simu1;
                            end
                            if lld2<lld2_max
                                ksim_group2(fig_order,:)=kk_simu2;
                                pini_group2(fig_order,:)=p_ini2;
                                lld2_max=lld2;
                                kk_simu2_best=kk_simu2;
                            end
                            lld_group(fig_order,:)=[lld1_max,lld2_max];
                            [~,P_dtb1,~,p_ini1] =  FSP_threestate(kk_simu1_best,f_ob,p_ini,[T(1:2),T(3)+T(4)],hold_time_1);
%                             [~,P_dtb2,~,p_ini2] =  FSP_threestate(kk_simu2_double,f_ob_2,p_ini,[T(end-3)+T(end-2),T(end-1:end)],hold_time);
                            [~,P_dtb2,~,p_ini2] =  FSP_threestate(kk_simu2_best,f_ob_2,p_ini,T_cds_2,hold_time);
                        %% scan the parameters' space
                            if k_scan==1
%                                     [lld1_group,lld2_group,hb1_k_best,hb7_k_best,P1_scan,P2_scan]=...
%                                         FitScan_ratio(kk_simu1_best,kk_simu2_best,f_ob,f_ob_2,p_ini,T,hold_time_1,hold_time);
%                                     hb1_scan_g(fig_order,:)=hb1_k_best;hb7_scan_g(fig_order,:)=hb7_k_best;
                                [lld1_group,lld2_group,hb1_k_best,hb7_k_best,P1_scan_k,P2_scan_k]=...
                                    FitScan(kk_simu1_double,kk_simu2_double,f_ob,f_ob_2,p_ini,T,hold_time_1,hold_time,i_cycle);
%                                 [lld1_group,lld2_group,hb1_k_best,hb7_k_best,P1_scan_k,P2_scan_k]=...
%                                     FitScan_remove(kk_simu1_double,kk_simu2_double,f_ob,f_ob_2,p_ini,T,hold_time_1,hold_time,i_cycle,in0f);
%                                 [lld1_group,lld2_group,hb1_k_best,hb7_k_best,P1_scan_k,P2_scan_k]=...
%                                     FitScan_OffFix(kk_simu1_best,kk_simu2_best,f_ob,f_ob_2,p_ini,T,hold_time_1,hold_time);
%                                 [lld1_group,lld2_group,hb1_k_best,hb7_k_best,P1_scan_k,P2_scan_k]=...
%                                     FitScan_FixRatio(kk_simu1_best,kk_simu2_best,f_ob,f_ob_2,p_ini,T,hold_time_1,hold_time,fig_title,fig_order);
                                %% scan fit curve
                                [lld1_best,P1_scan,~,~] =  FSP_threestate(P1_scan_k,f_ob,p_ini,[T(1:2),T(3)+T(4)],hold_time_1);
                                [lld2_best,P2_scan,~,~] =  FSP_threestate(P2_scan_k,f_ob_2,p_ini,T_cds_2,hold_time);
                                lld1_cell{fig_order}=lld1_group;lld2_cell{fig_order}=lld2_group;
                                hb1_scan_g(fig_order,:)=hb1_k_best(1,[1,3])./hb1_k_best(1,[2,4]);
                                hb7_scan_g(fig_order,:)=hb7_k_best(1,[1,3])./hb7_k_best(1,[2,4]);
                                lld_scan(fig_order,:)=[lld1_best,lld2_best];
                                
                                hb1_scan_k(fig_order,:)=hb1_k_best;
                                hb7_scan_k(fig_order,:)=hb7_k_best;   
                            end
                    end
                end
                %% Two-dimension data and fit figure 
%                 P_dtb=P_dtb1'*P_dtb2;
%                 imagesc(0:length(P_dtb)-1,0:length(P_dtb)-1,log10((P_dtb)'));
%                 xlabel('# HB 1')
%                 ylabel('# HB 7')
%                 hold on
%                 scatter(f_ob,f_ob_2,10,'filled','k')
%                 cb = colorbar;
%                 cb.Label.String = 'Frequence';
%                 caxis([1e-5,1])
%                 axis square
%                 xlim([0,30])
%                 ylim([0,30])
%                 title('FIT-DATA-2D')
%                 legend('Data')
%                 saveas(gcf,[[out_folder,'\'],'three_state_2D-',num2str(fig_order),'-',num2str(EL(1,fig_order)),'-',num2str(EL(2,fig_order)),'.fig']);
%                 saveas(gcf,[[out_folder,'\'],'three_state_2D-',num2str(fig_order),'-',num2str(EL(1,fig_order)),'-',num2str(EL(2,fig_order)),'.png']);
%                 close(gcf);
                %% data and fit curve figure
                Nbin_data=0:1:100;nb=Nbin_data;
                nf_2 = histc(reshape(f_ob_2,1,length(f_ob_2)),Nbin_data);
                nf_all_2=nf_2'/sum(nf_2);
                nf_1 = histc(reshape(f_ob,1,length(f_ob)),Nbin_data);
                nf_all_1=nf_1'/sum(nf_1);

                subplot(1,2,1)
                    bar(nb,nf_all_1,'g');
                    ymax1_2=sort(nf_all_1);ymax1_2=ymax1_2(end-1);
                    ylim([0 ymax1_2+0.01])
                    xlim([-0.5 40])
                    hold on
                    plot(0:50,P_dtb1(1:51),'b','LineWidth',1)
                    xlabel('# RNA1')
                    ylabel('frequency')
                    axis square
                    text(30,0.01,['\alpha = ',num2str(kk_simu1_double(end))])
                subplot(1,2,2)
                    bar(nb,nf_all_2,'y');
                    ymax2_2=sort(nf_all_2);ymax2_2=ymax2_2(end-1);
                    ylim([0 ymax2_2+0.01])
                    xlim([-0.5 40])
                    hold on
                    plot(0:50,P_dtb2(1:51),'r','LineWidth',1)
                    xlabel('# RNA2')
                    ylabel('frequency')
                    axis square
                    text(30,0.01,['\alpha = ',num2str(kk_simu1_double(end))])
                saveas(gcf,[[out_folder,'\'],'three_state_double-',num2str(fig_order),'-',num2str(EL(1,fig_order)),'-',num2str(EL(2,fig_order)),'.fig']);
                saveas(gcf,[[out_folder,'\'],'three_state_double-',num2str(fig_order),'-',num2str(EL(1,fig_order)),'-',num2str(EL(2,fig_order)),'.png']);
                close(gcf);
                if k_scan==1
                    subplot(1,2,1)
                        bar(nb,nf_all_1,'g');
                        ymax1_2=sort(nf_all_1);ymax1_2=ymax1_2(end-1);
                        ylim([0 ymax1_2+0.01])
                        xlim([-0.5 40])
                        hold on
                        plot(0:50,P1_scan(1:51),'b','LineWidth',1)
                        xlabel('# RNA1')
                        ylabel('frequency')
                        axis square
                        text(30,0.01,['\alpha = ',num2str(kk_simu1_double(end))])
                    subplot(1,2,2)
                        bar(nb,nf_all_2,'y');
                        ymax2_2=sort(nf_all_2);ymax2_2=ymax2_2(end-1);
                        ylim([0 ymax2_2+0.01])
                        xlim([-0.5 40])
                        hold on
                        plot(0:50,P2_scan(1:51),'r','LineWidth',1)
                        xlabel('# RNA2')
                        ylabel('frequency')
                        axis square
                        text(30,0.01,['\alpha = ',num2str(kk_simu1_double(end))])
                    saveas(gcf,[[out_folder,'\'],'three_state_double-scan',num2str(fig_order),'-',num2str(EL(1,fig_order)),'-',num2str(EL(2,fig_order)),'.fig']);
                    saveas(gcf,[[out_folder,'\'],'three_state_double-scan',num2str(fig_order),'-',num2str(EL(1,fig_order)),'-',num2str(EL(2,fig_order)),'.png']);
                    close(gcf);
                end
  
                Nbin_data_2=0:2:100;
                nf_2_bin2 = histc(reshape(f_ob_2,1,length(f_ob_2)),Nbin_data_2);
                nf_all_2=nf_2_bin2'/sum(nf_2_bin2);
                nf_1_bin2 = histc(reshape(f_ob,1,length(f_ob)),Nbin_data_2);
                nf_all_1=nf_1_bin2'/sum(nf_1_bin2);
                P_dtb1_bin2=sum([[P_dtb1,0];[0,P_dtb1]]);
                P_dtb1_bin2=P_dtb1_bin2(2:2:end); 
                P_dtb2_bin2=sum([[P_dtb2,0];[0,P_dtb2]]);
                P_dtb2_bin2=P_dtb2_bin2(2:2:end); 
                subplot(1,2,1)
                    bar(Nbin_data_2,nf_all_1,'g');
                    ylim([0 0.1])
                    xlim([-0.5 40])
                    hold on
                    plot(0:2:40,P_dtb1_bin2(1:21),'b','LineWidth',1)
                    xlabel('# RNA1')
                    ylabel('frequency')
                subplot(1,2,2)
                    bar(Nbin_data_2,nf_all_2,'y');
                    ylim([0 0.1])
                    xlim([-0.5 40])
                    hold on
                    plot(0:2:40,P_dtb2_bin2(1:21),'r','LineWidth',1)
                    xlabel('# RNA2')
                    ylabel('frequency')
                saveas(gcf,[[out_folder,'\'],'three_state_bin2-',num2str(EL(1,fig_order)),'-',num2str(EL(2,fig_order)),'.fig']);
                saveas(gcf,[[out_folder,'\'],'three_state_bin2-',num2str(EL(1,fig_order)),'-',num2str(EL(2,fig_order)),'.png']);
                close(gcf);
            end   
        end
    end
    
    %% k_scan fix analysis all el
%     lld1_cell{fig_order}=lld1_group;lld2_cell{fig_order}=lld2_group;
    hb1_scan_fix=zeros(length(data_all),4);hb7_scan_fix=zeros(length(data_all),4);
%     lld_allel_max=inf;
%     for index_k10 = 1:size(lld1_cell{1},2) 
%         for index_k21 = 1:size(lld1_cell{1},4) 
%             lld_allel=0;
%             for el_tr = 1:length(data_all)
%                 lld1_group=lld1_cell{el_tr};
%                 lld_allel=lld_allel+min(min(lld1_group(:,index_k10,:,index_k21)));
%             end
%             if lld_allel<lld_allel_max
%                 index_best=[index_k10,index_k21];
%                 lld_allel_max=lld_allel;
%             end
%         end
%     end
%     hb1_k01 = 10.^[-3:0.075:0];
%     hb1_k10 = 10.^[-1.8:0.075:0];
%     hb1_k12 = 10.^[-3:0.075:0];
%     hb1_k21 = 10.^[-1:0.06:0.2];
%     for el_tr = 1:length(data_all) 
%         lld1_group=lld1_cell{el_tr};
%         [index_k01,index_k12]=find(lld1_group(:,index_best(1),:,index_best(2))==min(min(lld1_group(:,index_best(1),:,index_best(2)))));
%         hb1_scan_fix(el_tr,:)=[hb1_k01(index_k01),hb1_k10(index_best(1)),hb1_k12(index_k12),hb1_k21(index_best(2))];        
%     end
%             
%     lld_allel_max=inf;
%     for index_k10 = 1:size(lld2_cell{1},2) 
%         for index_k21 = 1:size(lld2_cell{1},4) 
%             lld_allel=0;
%             for el_tr = 1:length(data_all)
%                 lld2_group=lld2_cell{el_tr};
%                 lld_allel=lld_allel+min(min(lld2_group(:,index_k10,:,index_k21)));
%             end
%             if lld_allel<lld_allel_max
%                 index_best=[index_k10,index_k21];
%                 lld_allel_max=lld_allel;
%             end
%         end
%     end
%     hb7_k01 = 10.^[-3:0.075:0];
%     hb7_k10 = 10.^[-1.8:0.075:0];
%     hb7_k12 = 10.^[-3:0.075:0];
%     hb7_k21 = 10.^[-1:0.06:0.5];
%     for el_tr = 1:length(data_all) 
%         lld2_group=lld2_cell{el_tr};
%         [index_k01,index_k12]=find(lld2_group(:,index_best(1),:,index_best(2))==min(min(lld2_group(:,index_best(1),:,index_best(2)))));
%         hb7_scan_fix(el_tr,:)=[hb7_k01(index_k01(1)),hb7_k10(index_best(1)),hb7_k12(index_k12(1)),hb7_k21(index_best(2))];        
%     end
    
    %% Hill fit of k
%     if simu==1
%         k_fit1=ksim_group1(:,[1,2])'./ksim_group1(:,[3,5])';
%         k_fit2=ksim_group2(:,[1,2])'./ksim_group2(:,[3,5])';
%         hill_qun = fittype( @(h,c,a,d,x) a*x.^h./(c^h+x.^h)+d);
%         Hill_fit11 = fit(mean(protein_concentration,1)',k_fit1(1,:)',hill_qun,'Lower',[0,0,0,0],'Upper',[inf,10e-6,inf,0.1]);
%         Hill_fit12 = fit(mean(protein_concentration,1)',k_fit1(2,:)',hill_qun,'Lower',[0,0,0,0],'Upper',[inf,10e-6,inf,0.1]);
%         Hill_fit21 = fit(mean(protein_concentration,1)',k_fit2(1,:)',hill_qun,'Lower',[0,0,0,0],'Upper',[inf,10e-6,inf,0.1]);
%         Hill_fit22 = fit(mean(protein_concentration,1)',k_fit2(2,:)',hill_qun,'Lower',[0,0,0,0],'Upper',[inf,10e-6,inf,0.1]);
%         subplot(2,2,1)
%             scatter(mean(protein_concentration,1)',k_fit1(1,:)',15,'filled','k')
%             hold on
%             plot(Hill_fit11,'r');
%             legend('k01/k10','Hill fit')
%             xlabel('protein concentration')
%             ylabel('ratio')
%             title('HB1')
%             text(0.6*max(xlim),0.2*max(ylim),{['h = ',num2str(Hill_fit11.h)],['Ka = ',num2str(Hill_fit11.c)],...
%                 ['a = ',num2str(Hill_fit11.a)],['d = ',num2str(Hill_fit11.d)]})
%         subplot(2,2,2)
%             scatter(mean(protein_concentration,1)',k_fit1(2,:)',15,'filled','k')
%             hold on
%             plot(Hill_fit12,'r');
%             legend('k02/k20','Hill fit')
%             xlabel('protein concentration')
%             ylabel('ratio')
%             title('HB1')
%             text(0.6*max(xlim),0.2*max(ylim),{['h = ',num2str(Hill_fit12.h)],['Ka = ',num2str(Hill_fit12.c)],...
%                 ['a = ',num2str(Hill_fit12.a)],['d = ',num2str(Hill_fit12.d)]})
%         subplot(2,2,3)
%             scatter(mean(protein_concentration,1)',k_fit2(1,:)',15,'filled','k')
%             hold on
%             plot(Hill_fit21,'r');
%             legend('k01/k10','Hill fit')
%             xlabel('protein concentration')
%             ylabel('ratio')
%             title('HB7')
%             text(0.6*max(xlim),0.2*max(ylim),{['h = ',num2str(Hill_fit21.h)],['Ka = ',num2str(Hill_fit21.c)],...
%                 ['a = ',num2str(Hill_fit21.a)],['d = ',num2str(Hill_fit21.d)]})
%         subplot(2,2,4)
%             scatter(mean(protein_concentration,1)',k_fit2(2,:)',15,'filled','k')
%             hold on
%             plot(Hill_fit22,'r');
%             legend('k02/k20','Hill fit')
%             xlabel('protein concentration')
%             ylabel('ratio')
%             title('HB7')
%             text(0.6*max(xlim),0.2*max(ylim),{['h = ',num2str(Hill_fit22.h)],['Ka = ',num2str(Hill_fit22.c)],...
%                 ['a = ',num2str(Hill_fit22.a)],['d = ',num2str(Hill_fit22.d)]})
%         saveas(gcf,[[out_folder,'\'],'Hill-fit','.fig']);
%         saveas(gcf,[[out_folder,'\'],'Hill-fit','.png']);
%         close(gcf);
%     end
    %% fig x-axis el\protein
    for fig_x={'el','protein'}
%     for fig_x={'el'}
        %% x-axis
        if strcmp(fig_x{1},'protein')
            el=protein_concentration;
        else
            el=EL;
        end
        %% total signal Value fig
        plot(mean(el,1),hb_all(:,1:2)','-d','MarkerSize',5)
        legend('hb1','hb7')
        xlabel(fig_x{1})
        ylabel('total signal Value')
        saveas(gcf,[[out_folder,'\'],'total signal Value','-',fig_x{1},'.fig']);
        saveas(gcf,[[out_folder,'\'],'total signal Value','-',fig_x{1},'.png']);
        close(gcf);
        %% Guass fit fig
        if gaussfit==1
            save([out_folder,'\',['Guass2 var-divide-mean ','.mat']],'lambda_group')
            subplot(1,2,1)
                plot(mean(el,1),lambda_group(:,1:2)','-d','MarkerSize',5)
                legend('guass peak weak','guass peak strong')
                xlabel(fig_x{1})
                ylabel('var-divide-mean Value')
                title('Hb1')
                ylim([0 5])
            subplot(1,2,2)
                plot(mean(el,1),lambda_group(:,3:4)','-d','MarkerSize',5)
                legend('guass peak weak','guass peak strong')
                xlabel(fig_x{1})
                ylabel('var-divide-mean Value')
                title('Hb7')
                ylim([0 5])
            saveas(gcf,[[out_folder,'\'],'var-divide-mean_hb17','-',fig_x{1},'.fig']);
            saveas(gcf,[[out_folder,'\'],'var-divide-mean_hb17','-',fig_x{1},'.png']);
            close(gcf);

            subplot(1,2,1)
                plot(mean(el,1),guass_group(:,1:2)','-d','MarkerSize',5)
                legend('guass peak weak','guass peak strong')
                xlabel(fig_x{1})
                ylabel('Intensity Value')
                title('Hb1')
            subplot(1,2,2)
                plot(mean(el,1),guass_group(:,3:4)','-d','MarkerSize',5)
                legend('guass peak weak','guass peak strong')
                xlabel(fig_x{1})
                ylabel('Intensity Value')
                title('Hb7')
            saveas(gcf,[[out_folder,'\'],'guassIntensity_hb17','-',fig_x{1},'.fig']);
            saveas(gcf,[[out_folder,'\'],'guassIntensity_hb17','-',fig_x{1},'.png']);
            close(gcf);
        end
        %% Enrichment fig
        if enrich==1
            save([out_folder,'\',['Enrichment','.mat']],'enrich_peak_gr','peak_enrich11_or','peak_enrich12_or','peak_enrich21_or','peak_enrich22_or')
            subplot(1,2,1)
                plot(mean(el,1),enrich_peak_gr(:,1:2)','-d','MarkerSize',5)
                legend('hb1 weak','hb1 strong')
                xlabel(fig_x{1})
                ylabel('Enrichment Value')
                title('Hb1')
            subplot(1,2,2)
                plot(mean(el,1),enrich_peak_gr(:,3:4)','-d','MarkerSize',5)
                legend('hb7 weak','hb7 strong')
                xlabel(fig_x{1})
                ylabel('Enrichment Value')
                title('Hb7')
            saveas(gcf,[[out_folder,'\'],'Enrichment_state_hb17','-',fig_x{1},'.fig']);
            saveas(gcf,[[out_folder,'\'],'Enrichment_state_hb17','-',fig_x{1},'.png']);
            close(gcf);

            subplot(1,2,1)
                scatter(mean(el,1),enrich_peak_gr(:,1),'b')
                legend('hb1 weak')
                hold on
                scatter(mean(el,1),enrich_peak_gr(:,2)','r')
                legend('hb1 strong')
                xlabel(fig_x{1})
                ylabel('Enrichment Value')
                title('Hb1')
            subplot(1,2,2)
                scatter(mean(el,1),enrich_peak_gr(:,3)','b')
                legend('hb7 weak')
                hold on
                scatter(mean(el,1),enrich_peak_gr(:,4)','r')
                legend('hb7 strong')
                xlabel(fig_x{1})
                ylabel('Enrichment Value')
                title('Hb7')
            saveas(gcf,[[out_folder,'\'],'Enrichment_state_hb17sc','-',fig_x{1},'.fig']);
            saveas(gcf,[[out_folder,'\'],'Enrichment_state_hb17sc','-',fig_x{1},'.png']);
            close(gcf);

            plot(mean(el,1),enrich_peak_gr(:,[1,3])','-d','MarkerSize',5)
            legend('hb1 weak','hb7 weak')
            xlabel(fig_x{1})
            ylabel('Enrichment Value')
            saveas(gcf,[[out_folder,'\'],'Enrichment_rna_weak','-',fig_x{1},'.fig']);
            saveas(gcf,[[out_folder,'\'],'Enrichment_rna_weak','-',fig_x{1},'.png']);
            close(gcf);

            plot(mean(el,1),enrich_peak_gr(:,[2,4])','-d','MarkerSize',5)
            legend('hb1 strong','hb7 strong')
            xlabel(fig_x{1})
            ylabel('Enrichment Value')
            saveas(gcf,[[out_folder,'\'],'Enrichment_rna_strong','-',fig_x{1},'.fig']);
            saveas(gcf,[[out_folder,'\'],'Enrichment_rna_strong','-',fig_x{1},'.png']);
            close(gcf);
        end
        %% model k parameters fig
        if simu==1
            switch model
                case 2
                    %% k-transform fig
                    save([out_folder,'\',['para_all','.mat']],...
                        'lld1_cell','lld2_cell','ksim_group1','ksim_group2','pini_group1','pini_group2')
                    subplot(1,2,1)
                        plot(mean(el,1),ksim_group1(:,[1,2])','-d','MarkerSize',5)
                        legend('k01','k10')
                        xlabel(fig_x{1})
                        ylabel('HB1 Knf value')
                    subplot(1,2,2)
                        plot(mean(el,1),ksim_group2(:,[1,2])','-d','MarkerSize',5)
                        legend('k01','k10')
                        xlabel(fig_x{1})
                        ylabel('HB7 Kon value')
                    saveas(gcf,[[out_folder,'\'],'k_tra','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'k_tra','-',fig_x{1},'.png']);
                    close(gcf);
                    %% KI value fig
                    subplot(2,2,[1,2])
                        plot(mean(el,1),ksim_group1(:,[3,4])','-d','MarkerSize',5)
                        legend('KI_{0}','KI_{1}')
                        xlabel(fig_x{1})
                        ylabel('HB1 KI value')
                    subplot(2,2,[3,4])
                        plot(mean(el,1),ksim_group2(:,[3,4])','-d','MarkerSize',5)
                        legend('KI_{0}','KI_{1}')
                        xlabel(fig_x{1})
                        ylabel('HB7 Kon value')
                    saveas(gcf,[[out_folder,'\'],'KI_tra','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'KI_tra','-',fig_x{1},'.png']);
                    close(gcf);
                    %% k-transform ratio fig
                    subplot(2,2,[1,2])
                        plot(mean(el,1),ksim_group1(:,2)'./ksim_group1(:,1)','-d','MarkerSize',5)
                        legend('k_{10}/k_{01}')
                        xlabel(fig_x{1})
                        ylabel('HB1 K_{ratio} value')
                    subplot(2,2,[3,4])
                        plot(mean(el,1),ksim_group2(:,2)'./ksim_group2(:,1)','-d','MarkerSize',5)
                        legend('k_{10}/k_{01}')
                        xlabel(fig_x{1})
                        ylabel('HB7 K_{ratio} value')
                    saveas(gcf,[[out_folder,'\'],'k_ratio','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'k_ratio','-',fig_x{1},'.png']);
                    close(gcf);

                    subplot(1,2,1)
                        plot(mean(el,1),[pini_group1(:,1)./pini_group1(:,2)]','-d','MarkerSize',5)
                        legend('n_{0}/n_{1}')
                        xlabel(fig_x{1})
                        ylabel('HB1 state ratio')
                    subplot(1,2,2)
                        plot(mean(el,1),[pini_group2(:,1)./pini_group2(:,2)]','-d','MarkerSize',5)
                        legend('n_{0}/n_{1}')
                        xlabel(fig_x{1})
                        ylabel('HB7 state ratio')
                    saveas(gcf,[[out_folder,'\'],'p_state_ratio','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'p_state_ratio','-',fig_x{1},'.png']);
                    close(gcf);
                    %% double copy ratio
                    subplot(1,2,1)
                        plot(mean(el,1),1-ksim_group1(:,5)','-d','MarkerSize',5)
                        legend('HB1')
                        xlabel(fig_x{1})
                        ylabel('P_{double copy}')
                    subplot(1,2,2)
                        plot(mean(el,1),1-ksim_group2(:,5)','-d','MarkerSize',5)
                        legend('HB7')
                        xlabel(fig_x{1})
                        ylabel('P_{double copy}')
                    saveas(gcf,[[out_folder,'\'],'k_pc','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'k_pc','-',fig_x{1},'.png']);
                    close(gcf);
                    %% state fig
                    subplot(1,2,1)
                        plot(mean(el,1),pini_group1','-d','MarkerSize',5)
                        legend('n0','n1')
                        xlabel(fig_x{1})
                        ylabel('HB1 state')
                    subplot(1,2,2)
                        plot(mean(el,1),pini_group2','-d','MarkerSize',5)
                        legend('n0','n1')
                        xlabel(fig_x{1})
                        ylabel('HB7 state')
                    saveas(gcf,[[out_folder,'\'],'p_state','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'p_state','-',fig_x{1},'.png']);
                    close(gcf);
                case 3
                    %% save parameter file
                    save([out_folder,'\',['para_all','.mat']],'el','protein_concentration'...
                        ,'lld1_cell','lld2_cell','lld_scan','lld_group','ksim_group1','ksim_group2','pini_group1','pini_group2')
                    %% k-transform fig
                    subplot(2,2,1)
                        plot(mean(el,1),ksim_group1(:,[1,2,4])','-d','MarkerSize',5)
                        legend('k01','k02','k12')
                        xlabel(fig_x{1})
                        ylabel('HB1 Kon value')
                    subplot(2,2,2)
                        plot(mean(el,1),ksim_group2(:,[1,2,4])','-d','MarkerSize',5)
                        legend('k01','k02','k12')
                        xlabel(fig_x{1})
                        ylabel('HB7 Kon value')
                    subplot(2,2,3)
                        plot(mean(el,1),ksim_group1(:,[3,5,6])','-d','MarkerSize',5)
                        legend('k10','k20','k21')
                        xlabel(fig_x{1})
                        ylabel('HB1 Koff value')
                    subplot(2,2,4)
                        plot(mean(el,1),ksim_group2(:,[3,5,6])','-d','MarkerSize',5)
                        legend('k10','k20','k21')
                        xlabel(fig_x{1})
                        ylabel('HB7 Koff value')
                    saveas(gcf,[[out_folder,'\'],'k_tra','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'k_tra','-',fig_x{1},'.png']);
                    close(gcf);
                    %% scan k-transform fig
                    subplot(2,2,1)
                        plot(mean(el,1),hb1_scan_k(:,[1,3])','-d','MarkerSize',5)
                        legend('k01','k12')
                        xlabel(fig_x{1})
                        ylabel('HB1 Kon value')
                    subplot(2,2,2)
                        plot(mean(el,1),hb7_scan_k(:,[1,3])','-d','MarkerSize',5)
                        legend('k01','k12')
                        xlabel(fig_x{1})
                        ylabel('HB7 Kon value')
                    subplot(2,2,3)
                        plot(mean(el,1),hb1_scan_k(:,[2,4])','-d','MarkerSize',5)
                        legend('k10','k21')
                        xlabel(fig_x{1})
                        ylabel('HB1 Koff value')
                    subplot(2,2,4)
                        plot(mean(el,1),hb7_scan_k(:,[2,4])','-d','MarkerSize',5)
                        legend('k10','k21')
                        xlabel(fig_x{1})
                        ylabel('HB7 Koff value')
                    saveas(gcf,[[out_folder,'\'],'k_tra_scan','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'k_tra_scan','-',fig_x{1},'.png']);
                    close(gcf); 
                    %% scan k-transform fix fig
                    subplot(2,2,1)
                        plot(mean(el,1),hb1_scan_fix(:,[1,3])','-d','MarkerSize',5)
                        legend('k01','k12')
                        xlabel(fig_x{1})
                        ylabel('HB1 Kon value')
                    subplot(2,2,2)
                        plot(mean(el,1),hb7_scan_fix(:,[1,3])','-d','MarkerSize',5)
                        legend('k01','k12')
                        xlabel(fig_x{1})
                        ylabel('HB7 Kon value')
                    subplot(2,2,3)
                        plot(mean(el,1),hb1_scan_fix(:,[2,4])','-d','MarkerSize',5)
                        legend('k10','k21')
                        xlabel(fig_x{1})
                        ylabel('HB1 Koff value')
                    subplot(2,2,4)
                        plot(mean(el,1),hb7_scan_fix(:,[2,4])','-d','MarkerSize',5)
                        legend('k10','k21')
                        xlabel(fig_x{1})
                        ylabel('HB7 Koff value')
                    saveas(gcf,[[out_folder,'\'],'k_tra_scan_fix','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'k_tra_scan_fix','-',fig_x{1},'.png']);
                    close(gcf); 
                    %% KI value fig
                    subplot(2,2,[1,2])
                        plot(mean(el,1),ksim_group1(:,[7,8,9])','-d','MarkerSize',5)
                        legend('KI_{0}','KI_{1}','KI_{2}')
                        xlabel(fig_x{1})
                        ylabel('HB1 KI value')
                    subplot(2,2,[3,4])
                        plot(mean(el,1),ksim_group2(:,[7,8,9])','-d','MarkerSize',5)
                        legend('KI_{0}','KI_{1}','KI_{2}')
                        xlabel(fig_x{1})
                        ylabel('HB7 KI value')
                    saveas(gcf,[[out_folder,'\'],'KI_tra','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'KI_tra','-',fig_x{1},'.png']);
                    close(gcf);
                    %% k-transform ratio fig
                    subplot(2,2,[1,2])
                        plot(mean(el,1),ksim_group1(:,[3,6])'./ksim_group1(:,[1,4])','-d','MarkerSize',5)
                        legend('k_{10}/k_{01}','k_{21}/k_{12}')
                        xlabel(fig_x{1})
                        ylabel('HB1 K_{ratio} value')
                    subplot(2,2,[3,4])
                        plot(mean(el,1),ksim_group2(:,[3,6])'./ksim_group2(:,[1,4])','-d','MarkerSize',5)
                        legend('k_{10}/k_{01}','k_{21}/k_{12}')
                        xlabel(fig_x{1})
                        ylabel('HB7 K_{ratio} value')
                    saveas(gcf,[[out_folder,'\'],'k_ratio','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'k_ratio','-',fig_x{1},'.png']);
                    close(gcf);                 
                    %% k-transform ratio fig reciprocal
                    subplot(2,2,[1,2])
                        plot(mean(el,1),ksim_group1(:,[1,4])'./ksim_group1(:,[3,6])','-d','MarkerSize',5)
                        legend('k_{01}/k_{10}','k_{12}/k_{21}')
                        xlabel(fig_x{1})
                        ylabel('HB1 K_{ratio} value')
                    subplot(2,2,[3,4])
                        plot(mean(el,1),ksim_group2(:,[1,4])'./ksim_group2(:,[3,6])','-d','MarkerSize',5)
                        legend('k_{01}/k_{10}','k_{12}/k_{21}')
                        xlabel(fig_x{1})
                        ylabel('HB7 K_{ratio} value')
                    saveas(gcf,[[out_folder,'\'],'k_ratio_re','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'k_ratio_re','-',fig_x{1},'.png']);
                    close(gcf);
                    %% scan k-transform ratio fig
                    subplot(2,2,[1,2])
                        plot(mean(el,1),hb1_scan_g(:,[1,2])','-d','MarkerSize',5)
                        legend('k_{01}/k_{10}','k_{12}/k_{21}')
                        xlabel(fig_x{1})
                        ylabel('HB1 K_{ratio} value')
                    subplot(2,2,[3,4])
                        plot(mean(el,1),hb7_scan_g(:,[1,2])','-d','MarkerSize',5)
                        legend('k_{01}/k_{10}','k_{12}/k_{21}')
                        xlabel(fig_x{1})
                        ylabel('HB7 K_{ratio} value')
                    saveas(gcf,[[out_folder,'\'],'k_ratio_scan','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'k_ratio_scan','-',fig_x{1},'.png']);
                    close(gcf);
                    %% k-transform ratio 0102 fig
                    subplot(2,2,[1,2])
                        plot(mean(el,1),ksim_group1(:,[1,2])'./ksim_group1(:,[3,5])','-d','MarkerSize',5)
                        legend('k_{01}/k_{10}','k_{02}/k_{20}')
                        xlabel(fig_x{1})
                        ylabel('HB1 K_{ratio} value')
                    subplot(2,2,[3,4])
                        plot(mean(el,1),ksim_group2(:,[1,2])'./ksim_group2(:,[3,5])','-d','MarkerSize',5)
                        legend('k_{01}/k_{10}','k_{02}/k_{20}')
                        xlabel(fig_x{1})
                        ylabel('HB7 K_{ratio} value')
                    saveas(gcf,[[out_folder,'\'],'k_ratio_0102','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'k_ratio_0102','-',fig_x{1},'.png']);
                    close(gcf);
                    %% double copy ratio
                    plot(mean(el,1),[ksim_group1(:,12),ksim_group2(:,12)]','-d','MarkerSize',5)
                    legend('\alpha_{hb1}','\alpha_{hb7}')
                    xlabel(fig_x{1})
                    ylabel('\alpha value')
                    saveas(gcf,[[out_folder,'\'],'alp','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'alp','-',fig_x{1},'.png']);
                    close(gcf);
                    %% state ratio fig
                    subplot(1,2,1)
                        plot(mean(el,1),[pini_group1(:,1)./pini_group1(:,2),pini_group1(:,2)./pini_group1(:,3)]','-d','MarkerSize',5)
                        legend('n_{0}/n_{1}','n_{1}/n_{2}')
                        xlabel(fig_x{1})
                        ylabel('HB1 state ratio')
                    subplot(1,2,2)
                        plot(mean(el,1),[pini_group2(:,1)./pini_group2(:,2),pini_group2(:,2)./pini_group2(:,3)]','-d','MarkerSize',5)
                        legend('n_{0}/n_{1}','n_{1}/n_{2}')
                        xlabel(fig_x{1})
                        ylabel('HB7 state ratio')
                    saveas(gcf,[[out_folder,'\'],'p_state_ratio','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'p_state_ratio','-',fig_x{1},'.png']);
                    close(gcf);
                    %% state fig
                    subplot(1,2,1)
                        plot(mean(el,1),pini_group1','-d','MarkerSize',5)
                        legend('n0','n1','n2')
                        xlabel(fig_x{1})
                        ylabel('HB1 state')
                    subplot(1,2,2)
                        plot(mean(el,1),pini_group2','-d','MarkerSize',5)
                        legend('n0','n1','n2')
                        xlabel(fig_x{1})
                        ylabel('HB7 state')
                    saveas(gcf,[[out_folder,'\'],'p_state','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'p_state','-',fig_x{1},'.png']);
                    close(gcf);
                    %% scan likelihood fig
                    plot(mean(el,1),lld_scan(:,[1,2])','-d','MarkerSize',5)
                    legend('RNA1','RNA2')
                    xlabel(fig_x{1})
                    ylabel('LIKELIHOOD')
                    saveas(gcf,[[out_folder,'\'],'lld_scan','-',fig_x{1},'.fig']);
                    saveas(gcf,[[out_folder,'\'],'lld_scan','-',fig_x{1},'.png']);
                    close(gcf);
            end
        end
    end
end