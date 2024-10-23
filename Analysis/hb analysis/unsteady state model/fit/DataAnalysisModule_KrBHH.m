function DataAnalysisModule_KrBHH(in0f,elreg,elbin,elmove,embryo,sort_ornot,on_off,simu,model,enrich,suffix)
%% A module to analysis the data

%% input parameter
%in0f: data path %i.e. 'S1\'
%elreg: embryo length range %i.e. [0 1]
%elbin: analytic embryo length range ones %i.e. 0.1
%elmove: the EL-move between two analytic ranges  %i.e. 0.05
%embryo: single embryo, sort_ornot: sort the loci, on_off:piossion, simu:FSP,model, enrich:enrichment
%suffix: output fold name
%% Initialize required parameters
EmbryoName_all=[];
data_del=0;data_num=[];
big_information=0;
close_distance=0.71;% judge loci close or not
foci_normal_all=[];
foci_RNA_profile_all=[];
message_distr_allembryo=[];
corr_colsesignal_allembryo=[];
protein_concentration_el_allembryo=[];
protein_concentration_el_10=[];
protein_concentration_el_11=[];
protein_concentration_el_12=[];
protein_concentration_el_13=[];
protein_concentration_el_14=[];
elmin=elreg(1):elmove:elreg(2)-elbin;
elmax=elmin+elbin;
if elmax(end)<elreg(2)
    elmax(end)=elreg(2);
end
el=[elmin;elmax];
Npool=length(el);
% parpool('local',Npool)
DataCorrelation_all=[];
cycle11_enrich=[];
cycle12_enrich=[];
cycle13_enrich=[];
cycle14_enrich=[];
p22_cycle11=[];
p22_cycle12=[];
p22_cycle13=[];
p22_cycle14=[];
p22_allcycle=[];
signal_num_cycle11=[];
signal_num_cycle12=[];
signal_num_cycle13=[];
signal_num_cycle14=[];
closesignal1_allembryo=[];
closesignal2_allembryo=[];
contribution_me_allembryo=[];
histogram1_close_intensity_allem=[];
histogram2_close_intensity_allem=[];
[histogram2_close_intensity_cycle11,histogram2_close_intensity_cycle12,...
    histogram2_close_intensity_cycle13,histogram2_close_intensity_cycle14]=deal([],[],[],[]);
[histogram1_close_intensity_cycle11,histogram1_close_intensity_cycle12,...
    histogram1_close_intensity_cycle13,histogram1_close_intensity_cycle14]=deal([],[],[],[]);
[DataCorrelation_11,DataCorrelation_12,...
    DataCorrelation_13,DataCorrelation_14]=deal([],[],[],[]);
data_oneembryo_allem=[];
p22_cycle11_name='';
p22_cycle12_name='';
p22_cycle13_name='';
p22_cycle14_name='';
p22_allcycle_name='';
data_all=cell(1,length(elmin));
data_11=cell(1,length(elmin));
data_12=cell(1,length(elmin));
data_13=cell(1,length(elmin));
data_14=cell(1,length(elmin));
data_all_noreset=cell(1,length(elmin));
data_11_noreset=cell(1,length(elmin));
data_12_noreset=cell(1,length(elmin));
data_13_noreset=cell(1,length(elmin));
data_14_noreset=cell(1,length(elmin));
data_oneembryo=cell(1,length(elmin));
ExpressionLoci=struct('allcycle',[],'cycle11',[],...
    'cycle12',[],'cycle13',[],'cycle14',[]);
ExpressionQuantity=struct('allcycle',[],'cycle11',[],...
    'cycle12',[],'cycle13',[],'cycle14',[]);
%% processing directory
% fff_name='Z:\FISHIF_WANG\confocal\';
fff_name='X:\kr-enhancer\';
% fff_name='Z:\FISHIF_WANG\confocal\others\';
mkfff_name='Y:\OtherWorks\Kr\';
Excel_deal='Duallist.xlsx';
% Excel_deal='S1.xls';
[~,~,Duallist]=xlsread([mkfff_name,Excel_deal]);
ff_name=[fff_name];
mkff_name=[mkfff_name,in0f];
Duallist=Duallist(cellfun(@(x)strcmp(x,in0f(1:end-1)),Duallist(:,1))&cellfun(@(x)strcmp(x,'T'),Duallist(:,7)),2);
%% output path
out0f=mkff_name;
out_folder = [out0f,'Output','\',mat2str(elreg),'_',num2str(elbin),'_',num2str(elmove),'_',suffix];
if simu == 1
    out_folder = [out0f,'Output','\',mat2str(elreg),'_',num2str(elbin),'_',num2str(elmove),'_','simu-model',num2str(model),suffix];
end
% mkdir(out_folder);
% mkdir([out_folder,'\scatter_dist'])
%% embryos to analysis pick in Duallist
for fold = 1:length(Duallist)
    %% file deal
    index=strfind(Duallist{fold},'/');
    in_folder_name0=Duallist{fold};
    in_folder_name=in_folder_name0(index(end-1)+1:end-1);
    mkdir([out_folder,'\',in_folder_name]);
    in_folder = [ff_name,in_folder_name,'\','Results_new\'];
    in_folder_st = [ff_name,in_folder_name,'\','stacks\'];
%     in_folder = [ff_name,in_folder_name,'\','Results2\'];
    in_folder_histogram=[ff_name,in_folder_name,'\','Histogram\'];
    in_folder_histogramA=[ff_name,in_folder_name,'\','Histogram_A\'];
    in_folder_histogram2=[ff_name,in_folder_name,'\','Histogram_RNA2\'];
    in_folder_histogram2A=[ff_name,in_folder_name,'\','Histogram_A_RNA2\'];
    in_folder_mask = [ff_name,in_folder_name,'\','masks\'];
    try
        try
            [~,~,matchlist]=xlsread([in_folder_st,'matchlist.xls']);    
        catch
            ff_name=strrep(ff_name,'\S3\','\S3_new\');
            in_folder = [ff_name,in_folder_name,'\','Results\'];
            in_folder_st = [ff_name,in_folder_name,'\','stacks\'];
        %     in_folder = [ff_name,in_folder_name,'\','Results2\'];
            in_folder_histogram=[ff_name,in_folder_name,'\','Histogram\'];
            in_folder_histogramA=[ff_name,in_folder_name,'\','Histogram_A\'];
            in_folder_histogram2=[ff_name,in_folder_name,'\','Histogram_RNA2\'];
            in_folder_histogram2A=[ff_name,in_folder_name,'\','Histogram_A_RNA2\'];
            in_folder_mask = [ff_name,in_folder_name,'\','masks\'];
            [~,~,matchlist]=xlsread([in_folder_st,'matchlist.xls']);  
        end
    catch
        try
            [~,~,matchlist]=xlsread([Duallist{fold},'stacks/matchlist.xls']);
        catch
            Duallist_new=cellfun(@(x) replace(x,'Z:','\\192.168.1.108\wangjingyaodata'),Duallist(:,1),'UniformOutput',false);
            Duallist(:,1)=Duallist_new;
            [~,~,matchlist]=xlsread([Duallist{fold},'stacks/matchlist.xls']);
        end
    end
    %% type a or b meaning cycle early or late
    %matchlist=matchlist(strcmp(matchlist(:,17),'b'),:);  
    %% embryo to analysis pick in matchlist
    for I_fit=1:size(matchlist,1)
        %% file deal
        EmbryoName=matchlist{I_fit,3}(1:end-1);
        EmbryoName_label=strrep(EmbryoName,'_','-');
        xy_unit=matchlist{I_fit,12};
        z_unit=matchlist{I_fit,14};
        try
            fname = [in_folder,[matchlist{I_fit,3}(1:end-1),'.mat']];
            fname_mask=[in_folder_mask,[matchlist{I_fit,3}(1:end-1),'\mask.mat']];
            data=cell(1,length(elmin));
            load(fname_mask,'mask_stack')
            try
                load(fname,'nucleus_RNA_profile','foci_RNA_profile')
                load(fname,'nucleus_signal2_profile','foci_signal2_profile')
                load(fname,'nucleus_protein_profile','nucleus_protein_profile_ab')
                load(fname,'foci_data','foci_data2')
                load(fname,'h','h2','r_size','r_size2')
            catch
                continue
            end
        catch
            fname = [Duallist{fold},'Results\',[matchlist{I_fit,3}(1:end-1),'.mat']];
            fname_mask=[Duallist{fold},'masks\',[matchlist{I_fit,3}(1:end-1),'\mask.mat']];
            data=cell(1,length(elmin));
            load(fname_mask,'mask_stack')
            load(fname,'nucleus_RNA_profile','foci_RNA_profile')
            load(fname,'nucleus_signal2_profile','foci_signal2_profile')
            load(fname,'nucleus_protein_profile','nucleus_protein_profile_ab')
            load(fname,'foci_data','foci_data2')
            load(fname,'h','h2','r_size','r_size2')
        end
        try
            N_cycle=matchlist{I_fit,18}; 
        catch
            N_cycle=matchlist{I_fit,16};
        end
        if N_cycle>14
            N_cycle=14;
        end
        if N_cycle<11
            N_cycle=11;
        end
      
        %% enrichment deal twice
%         enrich_fold='Z:\shihe zhang\S1\embryo_data_hb17\[0.2 0.5]_0.3_0.2__enrichembryo2\';
%         enrich_fo=[enrich_fold,in_folder_name,'\',EmbryoName,'\Enrichment.mat'];
%         load(enrich_fo)
%         eval(['cycle',num2str(N_cycle),'_enrich = [cycle',num2str(N_cycle),'_enrich;enrich_peak_gr]']);       
        %% hb-cds deal use
%         road_s3=[in_folder_name,'\',EmbryoName];     
% %         road_s3='cycle12';
%         p22_best_em=[];
%         parfor i=1:length(el)
%             A_Hb17_Cds_simu2(i,0.5375,0.1224,road_s3,N_cycle);
%             [p22_best]=A_Hb17_Cds_simu(i,0.5375,0.1224,road_s3,N_cycle);
%             p22_best_em=[p22_best_em;p22_best];
%         end
%         comm_p22=['p22_cycle',num2str(N_cycle),'=[p22_cycle',num2str(N_cycle),',p22_best_em];'];
%         comm_p22_name=['p22_cycle',num2str(N_cycle),'_name=char(p22_cycle',num2str(N_cycle),'_name,EmbryoName_label);'];
%         eval(comm_p22);
%         eval(comm_p22_name);
%         p22_allcycle=[p22_allcycle,p22_best_em];
%         p22_allcycle_name=char(p22_allcycle_name,EmbryoName_label);
%         save([out_folder,'\',['p22_data','.mat']],'p22_allcycle','p22_cycle11','p22_cycle12','p22_cycle13','p22_cycle14');
        %% close foci from histogram
        %% load foci rna profile
%         data=cell(1,length(elmin));
%         load(fname_mask,'mask_stack')
%         load(fname,'nucleus_RNA_profile','foci_RNA_profile')
%         load(fname,'nucleus_signal2_profile','foci_signal2_profile')
%         load(fname,'nucleus_protein_profile','nucleus_protein_profile_ab')
%         load(fname,'foci_data','foci_data2')
%         load(fname,'h','h2','r_size','r_size2')
        Histogram=0;
        if Histogram==1
            histogram_RNA = dir([in_folder_histogram,'*_raw.xls']);
            histogram_RNA_A = dir([in_folder_histogramA,'*fit.mat']);
            [~,~,histogram1]=xlsread([in_folder_histogram,histogram_RNA(I_fit).name]);
    %         if strcmp(in0f,'S43\')||strcmp(in0f,'S44\')||strcmp(in0f,'Control_Intron_yell_CDS_60X\')...
    %                 ||strcmp(in0f,'PE_remove_CDS_Intron_yellow_60X\')||strcmp(in0f,'S3_new\')
            if ~(strcmp(in0f,'S1\')||strcmp(in0f,'S3\'))
                try 
                    load([in_folder_histogramA,histogram_RNA_A(I_fit).name],'b');
                catch error
                    b=1;
                end
            else
                load([in_folder_histogramA,histogram_RNA_A(I_fit).name],'b');
            end
            histogram1=cell2mat(histogram1);
            histogram1_original=histogram1;
            histogram1(:,[6,7])=xy_unit*histogram1(:,[6,7]);histogram1(:,8)=z_unit*histogram1(:,8);

            % scatter of foci p1
    %         histogram1_intensity=histogram1(:,1).*histogram1(:,2).*histogram1(:,3)*2*pi/b;
    %         histogram1_scatter=[histogram1_intensity,histogram1(:,[6,7,8])];
    %         histogram1_scatter=histogram1_scatter(histogram1_scatter(:,1)<60,:);
    %         scatter(histogram1_scatter(:,2),histogram1_scatter(:,3),25,histogram1_scatter(:,1),'filled')
    %         xlabel('x')
    %         ylabel('y')
    %         title('hb1')
    %         grid on
    %         h = colorbar;
    %         colormap('cool')
    %         set(get(h,'label'),'string','intensity');
    %         saveas(gcf,[[out_folder,'\scatter_dist\'],[EmbryoName,'hb1','.fig']])
    %         saveas(gcf,[[out_folder,'\scatter_dist\'],[EmbryoName,'hb1','.png']])
    %         close(1)

            Distance_Loc=[];
            for i = 1:size(histogram1,1)
                Loc=histogram1(i,[6,7,8]);
                histogram1_rloc=(histogram1(:,[6,7,8])-Loc).^2;
                Distance_Loc=[Distance_Loc,sqrt(sum(histogram1_rloc,2))];
            end
            [Close_Index_x,Close_Index_y]=find(Distance_Loc<close_distance&Distance_Loc~=0);
            Close_Index=[Close_Index_x,Close_Index_y];
            Close_Index=sort(Close_Index,2);
            Close_Index=unique(Close_Index,'rows');
            Distance_close=Distance_Loc(Close_Index(:,1),Close_Index(:,2));
            Distance_close=diag(Distance_close);
            Distance_close=[Distance_close;Distance_close];
            histogram1_close=histogram1_original(Close_Index(:),:);
            histogram1_unpair=histogram1_original;histogram1_unpair(Close_Index(:),:)=[];
            histogram1_close_intensity=histogram1_close(:,1).*histogram1_close(:,2).*histogram1_close(:,3)*2*pi/b;
            histogram1_unpair_intensity=histogram1_unpair(:,1).*histogram1_unpair(:,2).*histogram1_unpair(:,3)*2*pi/b;

            histogram1_close_intensity=[histogram1_close_intensity,Distance_close];
            %%p2
            histogram_RNA2 = dir([in_folder_histogram2,'*_raw.xls']);
            histogram_RNA2_A = dir([in_folder_histogram2A,'*fit.mat']);
            [~,~,histogram2]=xlsread([in_folder_histogram2,histogram_RNA2(I_fit).name]);

            histogram2=cell2mat(histogram2);
            histogram2_original=histogram2;
            histogram2(:,[6,7])=xy_unit*histogram2(:,[6,7]);histogram2(:,8)=z_unit*histogram2(:,8);

            % scatter of foci p2
    %         histogram2_intensity=histogram2(:,1).*histogram2(:,2).*histogram2(:,3)*2*pi/b;
    %         histogram2_scatter=[histogram2_intensity,histogram2(:,[6,7,8])];
    %         histogram2_scatter=histogram2_scatter(histogram2_scatter(:,1)<60,:);
    %         scatter(histogram2_scatter(:,2),histogram2_scatter(:,3),25,histogram2_scatter(:,1),'filled')
    %         xlabel('x')
    %         ylabel('y')
    %         title('hb7')
    %         grid on
    %         h = colorbar;
    %         colormap('cool')
    %         set(get(h,'label'),'string','intensity');
    %         saveas(gcf,[[out_folder,'\scatter_dist\'],[EmbryoName,'hb7','.fig']])
    %         saveas(gcf,[[out_folder,'\scatter_dist\'],[EmbryoName,'hb7','.png']])
    %         close(1)

            Distance_Loc=[];
            for i = 1:size(histogram2,1)
                Loc=histogram2(i,[6,7,8]);
                histogram2_rloc=(histogram2(:,[6,7,8])-Loc).^2;
                Distance_Loc=[Distance_Loc,sqrt(sum(histogram2_rloc,2))];
            end
            [Close_Index_x,Close_Index_y]=find(Distance_Loc<close_distance&Distance_Loc~=0);
            Close_Index=[Close_Index_x,Close_Index_y];
            Close_Index=sort(Close_Index,2);
            Close_Index=unique(Close_Index,'rows');
            Distance_close2=Distance_Loc(Close_Index(:,1),Close_Index(:,2));
            Distance_close2=diag(Distance_close2);
            Distance_close2=[Distance_close2;Distance_close2];
            histogram2_close=histogram2_original(Close_Index(:),:);
            histogram2_unpair=histogram2_original;histogram2_unpair(Close_Index(:),:)=[];

    %         if strcmp(in0f,'S43\')||strcmp(in0f,'S44\')||strcmp(in0f,'Control_Intron_yell_CDS_60X\')...
    %                 ||strcmp(in0f,'PE_remove_CDS_Intron_yellow_60X\')||strcmp(in0f,'S3_new\')
            if ~(strcmp(in0f,'S1\')||strcmp(in0f,'S3\'))
                try 
                    load([in_folder_histogram2A,histogram_RNA2_A(I_fit).name],'b');
                catch error
                    b=1;
                end
            else
                load([in_folder_histogram2A,histogram_RNA2_A(I_fit).name],'b');
            end

            histogram2_close_intensity=histogram2_close(:,1).*histogram2_close(:,2).*histogram2_close(:,3)*2*pi/b;
            histogram2_unpair_intensity=histogram2_unpair(:,1).*histogram2_unpair(:,2).*histogram2_unpair(:,3)*2*pi/b;
            histogram2_close_intensity=[histogram2_close_intensity,Distance_close2];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%close signal with el
            index_1D_1=sub2ind(size(mask_stack),ceil(histogram1_close(:,6)),ceil(histogram1_close(:,7)),ceil(histogram1_close(:,8)));
            index_1D_1_unpair=sub2ind(size(mask_stack),ceil(histogram1_unpair(:,6)),ceil(histogram1_unpair(:,7)),ceil(histogram1_unpair(:,8)));
            cell_num1=mask_stack(index_1D_1);
            cell_num1(cell_num1==0)=1;
            cell_num1=[cell_num1(1:length(cell_num1)/2);cell_num1(1:length(cell_num1)/2)];
            close_el1=nucleus_RNA_profile(cell_num1,1);

            cell_num1=mask_stack(index_1D_1_unpair);
            cell_num1(cell_num1==0)=1;
            unpair_el1=nucleus_RNA_profile(cell_num1,1);

            histogram1_unpair_in_el=[histogram1_unpair_intensity,unpair_el1];
            histogram1_close_intensity=[histogram1_close_intensity,close_el1];

            index_1D_2=sub2ind(size(mask_stack),ceil(histogram2_close(:,6)),ceil(histogram2_close(:,7)),ceil(histogram2_close(:,8)));
            index_1D_2_unpair=sub2ind(size(mask_stack),ceil(histogram2_unpair(:,6)),ceil(histogram2_unpair(:,7)),ceil(histogram2_unpair(:,8)));

            cell_num2=mask_stack(index_1D_2);
            cell_num2(cell_num2==0)=1;
            cell_num2=[cell_num2(1:length(cell_num2)/2);cell_num2(1:length(cell_num2)/2)];
            close_el2=nucleus_signal2_profile(cell_num2,1);

            cell_num2=mask_stack(index_1D_2_unpair);
            cell_num2(cell_num2==0)=1;
            unpair_el2=nucleus_signal2_profile(cell_num2,1);
    %         cell_num2(cell_num2>length(nucleus_signal2_profile(:,1)))=1;

            histogram2_close_intensity=[histogram2_close_intensity,close_el2];
            histogram2_unpair_in_el=[histogram2_unpair_intensity,unpair_el2];
            histogram1_close_intensity_el=cell(1,length(elmin));
            histogram2_close_intensity_el=cell(1,length(elmin));
            histogram1_unpair_intensity_el=cell(1,length(elmin));
            histogram2_unpair_intensity_el=cell(1,length(elmin));
            protein_concentration_el=[];
            for iii = 1:length(elmin)
                el_I = (histogram1_close_intensity(:,end) >= elmin(iii)) & (histogram1_close_intensity(:,end) <= elmax(iii));
                el_I2 = (histogram2_close_intensity(:,end) >= elmin(iii)) & (histogram2_close_intensity(:,end) <= elmax(iii));
                histogram1_close_intensity_el{1,iii}=histogram1_close_intensity(el_I,:);
                histogram2_close_intensity_el{1,iii}=histogram2_close_intensity(el_I2,:);
                el_I = (histogram1_unpair_in_el(:,end) >= elmin(iii)) & (histogram1_unpair_in_el(:,end) <= elmax(iii));
                el_I2 = (histogram2_unpair_in_el(:,end) >= elmin(iii)) & (histogram2_unpair_in_el(:,end) <= elmax(iii));
                histogram1_unpair_intensity_el{1,iii}=histogram1_unpair_in_el(el_I,:);
                histogram2_unpair_intensity_el{1,iii}=histogram2_unpair_in_el(el_I2,:);
                %%%%%%%%%%%%%protein concentration
                if strcmp(in0f,'S43\')||strcmp(in0f,'S44\')||strcmp(in0f,'Control_Intron_yell_CDS_60X\')||strcmp(in0f,'PE_remove_CDS_Intron_yellow_60X\')||strcmp(in0f,'OreR_hb1hb7CDS_60X\')
                    protein_concentration_el = ones(1,length(elmin));
                elseif (strcmp(in0f,'S1\')||strcmp(in0f,'S3\'))
                    el_protein = (nucleus_protein_profile_ab(:,1) >= elmin(iii)) & (nucleus_protein_profile_ab(:,1) <= elmax(iii));
                    protein_concentration_el = [protein_concentration_el,mean(nucleus_protein_profile_ab(el_protein,5))];
                end
            end
            protein_concentration_el_allembryo=[protein_concentration_el_allembryo;protein_concentration_el];
        end
        protein_concentration_el=[];
        for iii = 1:length(elmin)
            if strcmp(in0f,'S43\')||strcmp(in0f,'B1\')
                protein_concentration_el = ones(1,length(elmin));
            else
                el_protein = (nucleus_protein_profile_ab(:,1) >= elmin(iii)) & (nucleus_protein_profile_ab(:,1) <= elmax(iii));
                protein_concentration_el = [protein_concentration_el,mean(nucleus_protein_profile_ab(el_protein,5))];
            end
        end
        protein_concentration_el_allembryo=[protein_concentration_el_allembryo;protein_concentration_el];
        EmbryoName_all=cat(1,EmbryoName_all,string(EmbryoName));
        eval(['protein_concentration_el_',num2str(N_cycle),'=[protein_concentration_el_',num2str(N_cycle),';protein_concentration_el];'])
        %% combine embryo
        if enrich == 1
            ind_foci = foci_RNA_profile(:,2);
            I_foci = foci_RNA_profile(:,3);
            N_nu = size(nucleus_RNA_profile,1);
            [n_ind,i_ind] = hist(ind_foci,1:N_nu);
            repet_index=find(n_ind>2);
            repet_all=[];
            for i=1:length(repet_index)
                repet_foci=find(ind_foci==repet_index(i));
                repet_all=[repet_all;repet_foci];
            end
            foci_RNA_profile(repet_all,:)=[];

            ind_foci = foci_signal2_profile(:,2);
            I_foci = foci_signal2_profile(:,3);
            N_nu = size(nucleus_signal2_profile,1);
            [n_ind,i_ind] = hist(ind_foci,1:N_nu);
            repet_index=find(n_ind>2);
            repet_all=[];
            for i=1:length(repet_index)
                repet_foci=find(ind_foci==repet_index(i));
                repet_all=[repet_all;repet_foci];
            end
            foci_signal2_profile(repet_all,:)=[];
            
            [~,eh_loc] = ismember(foci_RNA_profile(:,3),foci_data{1,5}(:,1));
            foci_rna = foci_RNA_profile(:,[1,3,5,4]);
%             foci_rna = [zeros(1,4);foci_rna];
%             eh_loc_logic=eh_loc;
%             eh_loc_logic(eh_loc_logic>0)=1;
%             eh_loc_logic=logical(eh_loc_logic);
            enrich_foci_value=foci_data{1,5}(:,[3,4,11,15,16,9,13]);
            enrich_foci_value=[zeros(1,size(enrich_foci_value,2));enrich_foci_value];
            enrich_foci_value=enrich_foci_value(eh_loc+1,:);       
%             enrich_foci_value=enrich_foci_value(eh_loc_logic,:);
            enrich_foci_rna = foci_rna;
%             enrich_foci_rna = enrich_foci_rna(eh_loc_logic,:);

            enrich_foci_value(:,1)=h_area*enrich_foci_value(:,6).*enrich_foci_value(:,7);
            enrich_foci_value_abs=((enrich_foci_value(:,2)-enrich_foci_value(:,1))./enrich_foci_value(:,4))./enrich_foci_value(:,5);
            enrich_foci_value_abs(isnan(enrich_foci_value_abs))=0;
    %         enrich_foci_value_abs(enrich_foci_value_abs<0)=0;
            %         enrich_foci_value_abs=enrich_foci_value(:,2);
    %         enrich_foci_value_abs(enrich_foci_value_abs<0)=0;
    %         foci_rna_profile_enrich=[enrich_foci_rna(:,1),enrich_foci_value(:,2)-enrich_foci_value(:,1),enrich_foci_rna(:,2),enrich_foci_rna(:,3)];
            foci_RNA_profile_enrich=[enrich_foci_rna(:,1),enrich_foci_value_abs,enrich_foci_rna(:,2),enrich_foci_rna(:,4),enrich_foci_rna(:,3)];        
            nucleus_protein_profile_ab(:,3) = 12;

            %%% signal 2
            [~,eh_loc2] = ismember(foci_signal2_profile(:,3),foci_data2{1,5}(:,1));
            foci_rna2 = foci_signal2_profile(:,[1,3,5,4]);
%             foci_rna2 = [zeros(1,4);foci_rna2];
%             eh_loc_logic=eh_loc2;
%             eh_loc_logic(eh_loc_logic>0)=1;
%             eh_loc_logic=logical(eh_loc_logic);
            enrich_foci_value2=foci_data2{1,5}(:,[3,4,11,15,16,9,13]);
            enrich_foci_value2=[zeros(1,size(enrich_foci_value2,2));enrich_foci_value2];
            enrich_foci_value2=enrich_foci_value2(eh_loc2+1,:);   
%             enrich_foci_value2=enrich_foci_value2(eh_loc_logic,:);
            enrich_foci_rna2 = foci_rna2;
%             enrich_foci_rna2 = enrich_foci_rna2(eh_loc_logic,:);

            enrich_foci_value2(:,1)=h_area2*enrich_foci_value2(:,6).*enrich_foci_value2(:,7);
            enrich_foci_value_abs2=((enrich_foci_value2(:,2)-enrich_foci_value2(:,1))./enrich_foci_value2(:,4))./enrich_foci_value2(:,5);
            enrich_foci_value_abs2(isnan(enrich_foci_value_abs2))=0;
    %         enrich_foci_value_abs2(enrich_foci_value_abs2<0)=0;
    %         enrich_foci_value_abs2=enrich_foci_value2(:,2);
    %         enrich_foci_value_abs2(enrich_foci_value_abs2<0)=0;
    %         foci_rna_profile_enrich=[enrich_foci_rna(:,1),enrich_foci_value(:,2)-enrich_foci_value(:,1),enrich_foci_rna(:,2),enrich_foci_rna(:,3)];
            foci_signal2_profile_enrich=[enrich_foci_rna2(:,1),enrich_foci_value_abs2,enrich_foci_rna2(:,2),enrich_foci_rna2(:,4),enrich_foci_rna2(:,3)];        
            [nucleus_RNA_profile0,ind_foci] = foci_info(nucleus_RNA_profile,foci_RNA_profile);%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [nucleus_signal2_profile0,ind_foci_2] = foci_info(nucleus_signal2_profile,foci_signal2_profile);
%             [nucleus_RNA_profile0,ind_foci] = foci_info_old(nucleus_RNA_profile,foci_RNA_profile);%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             [nucleus_signal2_profile0,ind_foci_2] = foci_info_old(nucleus_signal2_profile,foci_signal2_profile);
            
            for iii = 1:length(elmin)
                reg_I = (nucleus_RNA_profile(:,1) >= elmin(iii)) & (nucleus_RNA_profile(:,1) <= elmax(iii));
                reg_I2 = (nucleus_signal2_profile(:,1) >= elmin(iii)) & (nucleus_signal2_profile(:,1) <= elmax(iii));
                reg_I0 = reg_I(ind_foci);
                reg_I0_2 = reg_I2(ind_foci_2);
                pI0 = reg_I0;
                pI0_2 = reg_I0_2;

                f_ob = nucleus_RNA_profile0(pI0,4);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
                [~,loc]=ismember(f_ob,foci_RNA_profile_enrich(:,3));
                el1=foci_RNA_profile_enrich(:,[1,2]);
                el1=[zeros(1,2);el1];
%                 %not0=(f_ob~=0);
%                 is00=(f_ob==0);
%                 loc(is00)=0;
                f_ob_el=el1(loc+1,:);
                f_ob=[f_ob,f_ob_el];

                f_ob_2=nucleus_signal2_profile0(pI0_2,4);
                [is02,loc2]=ismember(f_ob_2,foci_signal2_profile_enrich(:,3));
%                 is002=(f_ob_2==0);
%                 loc2(is002)=0;
                el2=foci_signal2_profile_enrich(:,[1,2]);
                el2=[zeros(1,2);el2];
                f_ob_el2=el2(loc2+1,:);
                f_ob_2=[f_ob_2,f_ob_el2];
                f_enrich=f_ob(:,3);
                f_enrich2=f_ob_2(:,3);
                f_ob=f_ob(:,1);
                f_ob_2=f_ob_2(:,1);
                data_single=[f_ob,f_ob_2,f_enrich,f_enrich2];     
                data_mean = mean(data_single);
                if data_mean~=0
                    data{iii}=data_single;
                    switch N_cycle
                        case 11
                            if isempty(data_11{iii})
                                data_11{iii}=cat(1,data_11{iii},data{iii});
                            else
                                data_11_mean = mean(data_11{iii});
                                data_11{iii}=cat(1,data_11{iii},data{iii}./data_mean.*data_11_mean);
                            end
                            data_11_noreset{iii}=cat(1,data_11_noreset{iii},data{iii});
                        case 12
                            if isempty(data_12{iii})
                                data_12{iii}=cat(1,data_12{iii},data{iii});
                            else
                                data_12_mean = mean(data_12{iii});
                                data_12{iii}=cat(1,data_12{iii},data{iii}./data_mean.*data_12_mean);
                            end
                            data_12_noreset{iii}=cat(1,data_12_noreset{iii},data{iii});
                        case 13
                            if isempty(data_13{iii})
                                data_13{iii}=cat(1,data_13{iii},data{iii});
                            else
                                data_13_mean = mean(data_13{iii});
                                data_13{iii}=cat(1,data_13{iii},data{iii}./data_mean.*data_13_mean);
                            end
                            data_13_noreset{iii}=cat(1,data_13_noreset{iii},data{iii});
                        case 14
                            if isempty(data_14{iii})
                                data_14{iii}=cat(1,data_14{iii},data{iii});
                            else
                                data_14_mean = mean(data_14{iii});
                                data_14{iii}=cat(1,data_14{iii},data{iii}./data_mean.*data_14_mean);
                            end
                            data_14_noreset{iii}=cat(1,data_14_noreset{iii},data{iii});
                    end
                    if isempty(data_all{iii})
                        data_all{iii}=cat(1,data_all{iii},data{iii});
                    else
                        data_all_mean = mean(data_all{iii});
                        data_all{iii}=cat(1,data_all{iii},data{iii}./data_mean.*data_all_mean);
                    end
                    data_all_noreset{iii}=cat(1,data_all_noreset{iii},data{iii});
                end
            end
        else
            ind_foci = foci_RNA_profile(:,2);
            if N_cycle==13
                data_num=[data_num,size(foci_RNA_profile,1)];
            end
            I_foci = foci_RNA_profile(:,3);
            N_nu = size(nucleus_RNA_profile,1);
            [n_ind,i_ind] = hist(ind_foci,1:N_nu);
            repet_index=find(n_ind>2);
            repet_all=[];
            for i=1:length(repet_index)
                repet_foci=find(ind_foci==repet_index(i));
                repet_all=[repet_all;repet_foci];
            end
            foci_RNA_profile(repet_all,:)=[];
            data_del=data_del+length(repet_all);

            ind_foci = foci_signal2_profile(:,2);
            I_foci = foci_signal2_profile(:,3);
            N_nu = size(nucleus_signal2_profile,1);
            [n_ind,i_ind] = hist(ind_foci,1:N_nu);
            repet_index=find(n_ind>2);
            repet_all=[];
            for i=1:length(repet_index)
                repet_foci=find(ind_foci==repet_index(i));
                repet_all=[repet_all;repet_foci];
            end
            foci_signal2_profile(repet_all,:)=[];
            
%             [nucleus_RNA_profile0,ind_foci] = foci_info_old(nucleus_RNA_profile,foci_RNA_profile);%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             [nucleus_signal2_profile0,ind_foci_2] = foci_info_old(nucleus_signal2_profile,foci_signal2_profile);
            try
                [nucleus_RNA_profile0,ind_foci] = foci_info(nucleus_RNA_profile,foci_RNA_profile);%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [nucleus_signal2_profile0,ind_foci_2] = foci_info(nucleus_signal2_profile,foci_signal2_profile);
            catch error
                [nucleus_RNA_profile0,ind_foci] = foci_info_old(nucleus_RNA_profile,foci_RNA_profile);%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [nucleus_signal2_profile0,ind_foci_2] = foci_info_old(nucleus_signal2_profile,foci_signal2_profile);
            end
            signal_num=[length(foci_RNA_profile),length(foci_signal2_profile)];
            for iii = 1:length(elmin)
                reg_I = (nucleus_RNA_profile(:,1) >= elmin(iii)) & (nucleus_RNA_profile(:,1) <= elmax(iii));
                reg_I2 = (nucleus_signal2_profile(:,1) >= elmin(iii)) & (nucleus_signal2_profile(:,1) <= elmax(iii));
                reg_I0 = reg_I(ind_foci);
                reg_I0_2 = reg_I2(ind_foci_2);
                pI0 = reg_I0;
                pI0_2 = reg_I0_2;

                f_ob = nucleus_RNA_profile0(pI0,4);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
                f_ob_2=nucleus_signal2_profile0(pI0_2,4);

                if sort_ornot==1
                    [is0,loc]=ismember(f_ob,foci_RNA_profile(:,3));
                    el1=foci_RNA_profile(:,1);
                    el1=[0;el1];
%                     not0=(f_ob~=0);
                    is00=(f_ob==0);
                    loc(is00)=0;
                    f_ob_el=el1(loc+1);
                    f_ob=[f_ob,f_ob_el];
                    [is02,loc2]=ismember(f_ob_2,foci_signal2_profile(:,3));
                    is002=(f_ob_2==0);
                    loc2(is002)=0;
                    el2=foci_signal2_profile(:,1);
                    el2=[0;el2];
                    f_ob_el2=el2(loc2+1);
                    f_ob_2=[f_ob_2,f_ob_el2];
                    for ii=1:length(f_ob(:,1))
                        if rem(ii,2)==1
                            if f_ob(ii,1)~=0
                                if f_ob_2(ii,1)~=0
                                    if abs(f_ob(ii,2)-f_ob_2(ii,2))>0.001
                                        xx=f_ob_2(ii:ii+1,:);
                                        f_ob_2(ii,:)=xx(2,:);
                                        f_ob_2(ii+1,:)=xx(1,:);
                                    end
                                elseif f_ob_2(ii,1)==0 && f_ob_2(ii+1,1)~=0  
                                    if abs(f_ob(ii,2)-f_ob_2(ii+1,2))<0.001
                                        xx=f_ob_2(ii:ii+1,:);
                                        f_ob_2(ii,:)=xx(2,:);
                                        f_ob_2(ii+1,:)=xx(1,:);
                                    end
                                end
                            end 
                        else
                            if f_ob(ii,1)~=0
                                if f_ob_2(ii,1)~=0
                                    if abs(f_ob(ii,2)-f_ob_2(ii,2))>0.001
                                        xx=f_ob_2(ii-1:ii,:);
                                        f_ob_2(ii-1,:)=xx(2,:);
                                        f_ob_2(ii,:)=xx(1,:);
                                    end
                                elseif f_ob_2(ii,1)==0 && f_ob_2(ii-1,1)~=0  
                                    if abs(f_ob(ii,2)-f_ob_2(ii-1,2))<0.001
                                        xx=f_ob_2(ii-1:ii,:);
                                        f_ob_2(ii-1,:)=xx(2,:);
                                        f_ob_2(ii,:)=xx(1,:);
                                    end
                                end
                            end
                        end
                    end 
                end
                
                f_ob=f_ob(:,1);
                f_ob_2=f_ob_2(:,1);
                data_single=[f_ob,f_ob_2];
                data_mean = mean(data_single);
                if data_mean~=0
                    data{iii}=data_single;
                    switch N_cycle
                        case 11
                            if isempty(data_11{iii})
                                data_11{iii}=cat(1,data_11{iii},data{iii});
                            else
                                data_11_mean = mean(data_11{iii});
                                data_11{iii}=cat(1,data_11{iii},data{iii}./data_mean.*data_11_mean);
                            end
                            data_11_noreset{iii}=cat(1,data_11_noreset{iii},data{iii});
                        case 12
                            if isempty(data_12{iii})
                                data_12{iii}=cat(1,data_12{iii},data{iii});
                            else
                                data_12_mean = mean(data_12{iii});
                                data_12{iii}=cat(1,data_12{iii},data{iii}./data_mean.*data_12_mean);
                            end
                            data_12_noreset{iii}=cat(1,data_12_noreset{iii},data{iii});
                        case 13
                            if isempty(data_13{iii})
                                data_13{iii}=cat(1,data_13{iii},data{iii});
                            else
                                data_13_mean = mean(data_13{iii});
                                data_13{iii}=cat(1,data_13{iii},data{iii}./data_mean.*data_13_mean);
                            end
                            data_13_noreset{iii}=cat(1,data_13_noreset{iii},data{iii});
                        case 14
                            if isempty(data_14{iii})
                                data_14{iii}=cat(1,data_14{iii},data{iii});
                            else
                                data_14_mean = mean(data_14{iii});
                                data_14{iii}=cat(1,data_14{iii},data{iii}./data_mean.*data_14_mean);
                            end
                            data_14_noreset{iii}=cat(1,data_14_noreset{iii},data{iii});
                    end
                    if isempty(data_all{iii})
                        data_all{iii}=cat(1,data_all{iii},data{iii});
                    else
                        data_all_mean = mean(data_all{iii});
                        data_all{iii}=cat(1,data_all{iii},data{iii}./data_mean.*data_all_mean);
                    end
                    data_all_noreset{iii}=cat(1,data_all_noreset{iii},data{iii});
                end
            end
        end
        
        out_fold_oneembryo=[out_folder,'\',in_folder_name,'\',EmbryoName];
        mkdir(out_fold_oneembryo);
        data_oneembryo={data};
        data_corr=data;
        data_exp=data;
        data_empty=cell2mat(cellfun(@isempty,data_corr,'UniformOutput',false));
        for ii=1:length(data_empty)
            if data_empty(ii);data_exp{1,ii}=[0 0];end
        end
        ExpLoci=cell2mat(cellfun(@(x) sum(sum(x~=0))/(2*size(x,1)),data_corr,'UniformOutput',false));
        ExpressionLoci.allcycle=[ExpressionLoci.allcycle;ExpLoci];
        ExpQuan=cell2mat(cellfun(@(x) mean(mean(x(:,1))),data_exp,'UniformOutput',false));
        ExpressionQuantity.allcycle=[ExpressionQuantity.allcycle;ExpQuan];
        eval(['ExpressionLoci.','cycle',num2str(N_cycle)...
            ,'=[ExpressionLoci.','cycle',num2str(N_cycle),';ExpLoci];'])
        eval(['ExpressionQuantity.','cycle',num2str(N_cycle)...
            ,'=[ExpressionQuantity.','cycle',num2str(N_cycle),';ExpQuan];'])
        %% set ~=0
%         data_empty=cell2mat(cellfun(@isempty,data_corr,'UniformOutput',false));
%         for ii=1:length(data_empty)
%             if data_empty(ii);data_corr{1,ii}=[nan nan];end
%         end
%         data_corr=cellfun(@(x) x(x(:,1)~=0&x(:,2)~=0,:),data_corr,'UniformOutput',false);
        for ii=1:length(data_empty)
            if data_empty(ii);data_corr{1,ii}=[nan nan];end
        end
        DataCorrelation=cell2mat(cellfun(@(x) corr(x(:,1),x(:,2)),data_corr,'UniformOutput',false));
        DataCorrelation_all=[DataCorrelation_all;DataCorrelation];
        eval(['DataCorrelation_',num2str(N_cycle),'=[DataCorrelation_',num2str(N_cycle),';DataCorrelation];'])
        %% close signal analysis embryo
%         [message_distr,corr_closesignal,closesignal1_allel_re,closesignal2_allel_re,contribution_me]=...
%             ModuleClosesignal(out_fold_oneembryo,data_oneembryo,histogram1_unpair_intensity_el,histogram2_unpair_intensity_el...
%             ,el,0,histogram1_close_intensity_el,histogram2_close_intensity_el);
%         message_distr_allembryo=[message_distr_allembryo;message_distr];
%         corr_colsesignal_allembryo=[corr_colsesignal_allembryo;corr_closesignal];
%         closesignal1_allembryo=[closesignal1_allembryo;closesignal1_allel_re];
%         closesignal2_allembryo=[closesignal2_allembryo;closesignal2_allel_re];
%         contribution_me_allembryo=[contribution_me_allembryo;contribution_me];
%         
%         histogram1_close_intensity_allem=[histogram1_close_intensity_allem;histogram1_close_intensity_el{1}];
%         histogram2_close_intensity_allem=[histogram2_close_intensity_allem;histogram2_close_intensity_el{1}];
%         data_oneembryo_allem=[data_oneembryo_allem;data_oneembryo{1}{1}];
%         eval(['histogram2_close_intensity_cycle',num2str(N_cycle),'=[histogram2_close_intensity_cycle',num2str(N_cycle),';histogram2_close_intensity_el{1}(:,1)];'])
%         eval(['histogram1_close_intensity_cycle',num2str(N_cycle),'=[histogram1_close_intensity_cycle',num2str(N_cycle),';histogram1_close_intensity_el{1}(:,1)];'])
%         eval(['signal_num_cycle',num2str(N_cycle),'=[signal_num_cycle',num2str(N_cycle),';signal_num];'])
        %% single embryo model fit
        if embryo==1 
%             ModuleFit_ini(out_fold_oneembryo,data_oneembryo,el,protein_concentration_el,on_off,simu,model,enrich,N_cycle);
%             ModuleFit('',out_fold_oneembryo,data_oneembryo,el,protein_concentration_el,on_off,simu,model,enrich);
%             if N_cycle==12
%                 ModuleFit_signalem_CDS(out_fold_oneembryo,data_oneembryo,el,protein_concentration_el,on_off,simu,model,enrich,N_cycle);
%             end
%             ModuleFit_Signalem_simple(out_fold_oneembryo,data_oneembryo,el,protein_concentration_el,on_off,simu,model,enrich,N_cycle);
            if N_cycle            
                ModuleFit_signalem_Kr(out_fold_oneembryo,data_oneembryo,el,protein_concentration_el,on_off,simu,model,enrich,N_cycle);
            end
        end
        %% store embryo-data
        save([[out_folder,'\',in_folder_name],'\',[EmbryoName,'_data.mat']],'data','el')
    end
end
%% Expression Loci/Quantity
% for Exp_mode={'Loci','Quantity'}
%     Exp_mode=Exp_mode{1};
%     index_c=0;
%     Box_c=[];
%     eval(['Expression_cell=struct2cell(Expression',Exp_mode,');'])
%     N_max=max(cell2mat(cellfun(@(x) length(x(:)),Expression_cell,'UniformOutput',false)));
%     for Cycle={'allcycle','cycle11','cycle12','cycle13'}
%         index_c=index_c+1;
%         Cycle=Cycle{1};
%         eval(['Exp_Cy=Expression',Exp_mode,'.',Cycle,';'])
%         Box_c(:,index_c)=[Exp_Cy(:);nan(N_max-size(Exp_Cy(:),1),1)];
%         StdCorr=nanstd(Exp_Cy,1);
%         fill_x=[mean(el),fliplr(mean(el))];
%         fill_y1=[min(Exp_Cy,[],1),fliplr(max(Exp_Cy,[],1))];
%         fill_x=fill_x(~isnan(fill_y1));fill_y1=fill_y1(~isnan(fill_y1));
%         fill_y2=[nanmean(Exp_Cy,1)-StdCorr,fliplr(nanmean(Exp_Cy,1)+StdCorr)];
%         fill_y2=fill_y2(~isnan(fill_y2));
%         figure
%         hold on
%         plot(mean(el),nanmean(Exp_Cy),'r','LineWidth',2)
%         fill(fill_x,...
%             fill_y1,'b','facealpha',.1,'edgecolor','none')
%         fill(fill_x,...
%             fill_y2,'c','facealpha',.1,'edgecolor','none')
%         xlabel('Embryo Length')
%         ylabel(['Expression',Exp_mode])
%         title(['-',num2str(size(Exp_Cy,1)),'- Count Cycle',Cycle,' embryos'])
%         legend('Mean','MinMax','Std')
% %         ylim([-0.1 1])
%         xlim([0 1])
%         saveas(gcf,[[out_folder,'\'],'Exp',Exp_mode,'_Cy-p1',Cycle,'.fig'])
%         saveas(gcf,[[out_folder,'\'],'Exp',Exp_mode,'_Cy-p1',Cycle,'.png'])
%         close(gcf)
% 
%         figure
%         hold on
%         Exp_Cy_T=Exp_Cy';
%         plot(mean(el),nanmean(Exp_Cy),'LineWidth',2)
%         %     scatter(repmat(mean(el),1,size(DataCorrelation_embryo,1)),...
%         %         DataCorrelation_embryo_T(:))
%         for i=1:size(Exp_Cy,1)
%             scatter(mean(el),...
%                 Exp_Cy_T(:,i))
%         end
%         xlabel('Embryo Length')
%         ylabel(['Expression',Exp_mode])
%         title(['-',num2str(size(Exp_Cy,1)),'- Count Cycle',Cycle,' embryos'])
%         legend('Mean')
% %         ylim([-0.1 1])
%         xlim([0 1])
%         saveas(gcf,[[out_folder,'\'],'Exp',Exp_mode,'_Cy-p2',Cycle,'.fig'])
%         saveas(gcf,[[out_folder,'\'],'Exp',Exp_mode,'_Cy-p2',Cycle,'.png'])
%         close(gcf)
%     end
%     figure
%     label_box=[{'All Cycle'},{'Cycle 11'},{'Cycle 12'},{'Cycle 13'}];
%     boxplot(Box_c,label_box)
%     ylabel(['Expression',Exp_mode])
%     hold on
% %     ylim([0 1])
%     title(['Expression ',Exp_mode,' BOX'])
%     set(gca,'FontSize',15,...
%            'FontWeight','bold','FontName','times new Roman')
%     saveas(gcf,[[out_folder,'\'],'Exp',Exp_mode,'_Cy-BOX',Cycle,'.fig'])
%     saveas(gcf,[[out_folder,'\'],'Exp',Exp_mode,'_Cy-BOX',Cycle,'.png'])
% end
% Correlation fig
% index_c=0;
% figure
% for Cycle={'all','11','12','13'}
%     index_c=index_c+1;
%     Cycle=Cycle{1};
%     eval(['DataCorrelation_embryo=DataCorrelation_',Cycle,';'])
%     StdCorr=nanstd(DataCorrelation_embryo,1);
%     fill_x=[mean(el),fliplr(mean(el))];
%     fill_y1=[min(DataCorrelation_embryo,[],1),fliplr(max(DataCorrelation_embryo,[],1))];
%     fill_x=fill_x(~isnan(fill_y1));fill_y1=fill_y1(~isnan(fill_y1));
%     fill_y2=[nanmean(DataCorrelation_embryo,1)-StdCorr,fliplr(nanmean(DataCorrelation_embryo,1)+StdCorr)];
%     fill_y2=fill_y2(~isnan(fill_y2));
%     figure
%     hold on
%     plot(mean(el),nanmean(DataCorrelation_embryo),'LineWidth',2)
%     fill(fill_x,...
%         fill_y1,'m','facealpha',.1,'edgecolor','none')
%     fill(fill_x,...
%         fill_y2,'r','facealpha',.1,'edgecolor','none')
%     xlabel('Embryo Length')
%     ylabel('Correlation')
%     title(['-',num2str(size(DataCorrelation_embryo,1)),'- Count Cycle',Cycle,' embryos'])
%     legend('Mean','MinMax','Std')
%     ylim([-0.1 1])
%     xlim([0 1])
%     saveas(gcf,[[out_folder,'\'],'CorrelationStatistic-p1',Cycle,'.fig'])
%     saveas(gcf,[[out_folder,'\'],'CorrelationStatistic-p1',Cycle,'.png'])
%     close(gcf)
%     
%     figure
%     hold on
%     DataCorrelation_embryo_T=DataCorrelation_embryo';
%     plot(mean(el),nanmean(DataCorrelation_embryo),'LineWidth',2)
% %     scatter(repmat(mean(el),1,size(DataCorrelation_embryo,1)),...
% %         DataCorrelation_embryo_T(:))
%     for i=1:size(DataCorrelation_embryo,1)
%         scatter(mean(el),...
%             DataCorrelation_embryo_T(:,i))
%     end
%     xlabel('Embryo Length')
%     ylabel('Correlation')
%     title(['-',num2str(size(DataCorrelation_embryo,1)),'- Count Cycle',Cycle,' embryos'])
%     legend('Mean')
%     ylim([-0.1 1])
%     xlim([0 1])
%     saveas(gcf,[[out_folder,'\'],'CorrelationStatistic-p2',Cycle,'.fig'])
%     saveas(gcf,[[out_folder,'\'],'CorrelationStatistic-p2',Cycle,'.png'])
%     close(gcf)
%     
%     subplot(2,4,index_c)
%     hold on
%     plot(mean(el),nanmean(DataCorrelation_embryo),'LineWidth',2)
%     fill(fill_x,...
%         fill_y1,'m','facealpha',.1,'edgecolor','none')
%     fill(fill_x,...
%         fill_y2,'r','facealpha',.1,'edgecolor','none')
%     xlabel('Embryo Length')
%     ylabel('Correlation')
%     title(['-',num2str(size(DataCorrelation_embryo,1)),'- Count Cycle',Cycle,' embryos'])
%     legend('Mean','MinMax','Std')
%     ylim([-0.1 1])
%     xlim([0 1])
%     subplot(2,4,index_c+4)
%     hold on
%     DataCorrelation_embryo_T=DataCorrelation_embryo';
%     plot(mean(el),nanmean(DataCorrelation_embryo),'LineWidth',2)
% %     scatter(repmat(mean(el),1,size(DataCorrelation_embryo,1)),...
% %         DataCorrelation_embryo_T(:))
%     for i=1:size(DataCorrelation_embryo,1)
%         scatter(mean(el),...
%             DataCorrelation_embryo_T(:,i))
%     end
%     xlabel('Embryo Length')
%     ylabel('Correlation')
%     title(['-',num2str(size(DataCorrelation_embryo,1)),'- Count Cycle',Cycle,' embryos'])
%     legend('Mean')
%     ylim([-0.1 1])
%     xlim([0 1])
% end
% saveas(gcf,[[out_folder,'\'],'CorrelationStatistic','.fig'])
% saveas(gcf,[[out_folder,'\'],'CorrelationStatistic','.png'])
%% signal num
% figure
% bar([mean(signal_num_cycle11(1,:)),mean(signal_num_cycle12(1,:)),mean(signal_num_cycle13(1,:))],0.3);hold on 
% errorbar([mean(signal_num_cycle11(1,:)),mean(signal_num_cycle12(1,:)),mean(signal_num_cycle13(1,:))],...
%     [std0(signal_num_cycle11(1,:)),std0(signal_num_cycle12(1,:)),std0(signal_num_cycle13(1,:))],'r','Linestyle','None','LineWidth',1);
% ylabel('signal number')
% set(gca,'xtick',1:3);
% set(gca,'xticklabel',{'cycle11','cycle12','cycle13'});
% title(in0f(1:end-1))
%% closesignal fig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%compareallsignal
% lambda_all=contribution_me_allembryo(message_distr_allembryo(:,7)>10&message_distr_allembryo(:,10)>10,1:2:end);
% lambda_all(1,:)=[];
% label_box=[{'HB1 OVERALL'},{'HB7 OVERALL'},{'HB1 COLSESIGNAL'},{'HB7 COLSESIGNAL'}];
% boxplot(lambda_all,label_box)
% hold on
% scatter(0.9+0.2*rand(size(lambda_all,1),1),lambda_all(:,1),12,'filled','k')
% scatter(1.9+0.2*rand(size(lambda_all,1),1),lambda_all(:,2),12,'filled','k')
% scatter(2.9+0.2*rand(size(lambda_all,1),1),lambda_all(:,3),12,'filled','k')
% scatter(3.9+0.2*rand(size(lambda_all,1),1),lambda_all(:,4),12,'filled','k')
% ylabel('Peak2 contribution','FontSize',15,...
%        'FontWeight','bold')
% hold on
% ylim([0 1])
% title('Peak Contribution BOX')
% 
% lambda_all=contribution_me_allembryo(message_distr_allembryo(:,7)>10&message_distr_allembryo(:,10)>10,:);
% lambda_all=contribution_me_allembryo(message_distr_allembryo(:,10)>10,:);
% lambda_all_mean=mean(lambda_all);
% lambda_all_std0=std0(lambda_all);
% figure
% subplot(1,2,1)
% bar(lambda_all_mean([1,5]),0.3);hold on 
% errorbar(lambda_all_mean([1,5]),lambda_all_std0([1,5]),'r','Linestyle','None','LineWidth',1);
% ylim([0 1])
% ylabel('Peak2 Contribution')
% set(gca,'xtick',1:2);
% set(gca,'xticklabel',{'HB1 unpaired','HB1 paired'});
% title('HB1')
% subplot(1,2,2)
% bar(lambda_all_mean([3,7]),0.3);hold on 
% errorbar(lambda_all_mean([3,7]),lambda_all_std0([3,7]),'r','Linestyle','None','LineWidth',1);
% ylim([0 1])
% ylabel('Peak2 Contribution')
% set(gca,'xtick',1:2);
% set(gca,'xticklabel',{'HB7 unpaired','HB7 paired'});
% title('HB7')
% barweb(lambda_all_mean,...
%     lambda_all_std0);
% ylabel('Peak2 Contribution')
% legend('HB1 OVERALL','HB7 OVERALL','HB1 COLSESIGNAL','HB7 COLSESIGNAL')
% ylim([0 1])
% set(gca,'xtick',1);
% set(gca,'xticklabel',{'data'});
% 
% message_distr_allembryo_num=message_distr_allembryo(message_distr_allembryo(:,7)>10&message_distr_allembryo(:,10)>10,:);
% lambda_all=message_distr_allembryo_num(:,2:3:end)./message_distr_allembryo_num(:,3:3:end);
% label_box=[{'HB1 OVERALL'},{'HB1 COLSESIGNAL'},{'HB7 OVERALL'},{'HB7 COLSESIGNAL'}];
% boxplot(lambda_all,label_box)
% hold on
% scatter(0.9+0.2*rand(size(lambda_all,1),1),lambda_all(:,1),12,'filled','k')
% scatter(1.9+0.2*rand(size(lambda_all,1),1),lambda_all(:,2),12,'filled','k')
% scatter(2.9+0.2*rand(size(lambda_all,1),1),lambda_all(:,3),12,'filled','k')
% scatter(3.9+0.2*rand(size(lambda_all,1),1),lambda_all(:,4),12,'filled','k')
% ylabel('\sigma^2/\mu','FontSize',15,...
%        'FontWeight','bold')
% hold on
% line([0 5],[1 1],'Color','r','LineWidth',1)
% ylim([0 max(ylim)])
% text(0.7,1.5,'\sigma^2/\mu = 1')
% title('\sigma^2/\mu BOX')
% % 
% figure
% corr_colsesignal_allembryo_use=corr_colsesignal_allembryo(corr_colsesignal_allembryo(:,3)>=6&...
%     corr_colsesignal_allembryo(:,4)>=6,:);
% label_box_corr=[{'HB1 COLSESIGNAL CORR'},{'HB7 COLSESIGNAL CORR'}];
% boxplot(abs(corr_colsesignal_allembryo_use(:,[1,2])),label_box_corr)
% hold on
% scatter(0.9+0.2*rand(size(corr_colsesignal_allembryo_use,1),1),abs(corr_colsesignal_allembryo_use(:,1)),12,'filled','k')
% scatter(1.9+0.2*rand(size(corr_colsesignal_allembryo_use,1),1),abs(corr_colsesignal_allembryo_use(:,2)),12,'filled','k')
% ylim([0 1])
% 
% figure
% data1_allel=data_oneembryo_allem(:,1);data1_allelnz=data1_allel(data1_allel~=0);
% n_alldata1=length(data1_allelnz);
% var_1=var(data1_allelnz);mean_1=mean(data1_allelnz);
% 
% data2_allel=data_oneembryo_allem(:,2);data2_allelnz=data2_allel(data2_allel~=0);
% n_alldata2=length(data2_allelnz);
% var_2=var(data2_allelnz);mean_2=mean(data2_allelnz);
% 
% n_closesignal1_allel=length(histogram1_close_intensity_allem);
% var_1close=var(histogram1_close_intensity_allem);mean_1close=mean(histogram1_close_intensity_allem);
% n_closesignal2_allel=length(histogram2_close_intensity_allem);
% var_2close=var(histogram2_close_intensity_allem);mean_2close=mean(histogram2_close_intensity_allem);
% 
% f_oball=data_oneembryo_allem(:,1);
% f_oball(f_oball>100)=100;
% f_ob_2all=data_oneembryo_allem(:,2);
% f_ob_2all(f_ob_2all>100)=100;
% Nbin_data=0:1:100;
% nfall = histc(reshape(f_oball,1,length(f_oball)),Nbin_data);nb=Nbin_data;
% nf_2all = histc(reshape(f_ob_2all,1,length(f_ob_2all)),Nbin_data);
% nfall=nfall'/sum(nfall);nf_2all=nf_2all'/sum(nf_2all);
% nfclose1 = histc(reshape(histogram1_close_intensity_allem(:,1),1,length(histogram1_close_intensity_allem)),Nbin_data);
% nfclose2 = histc(reshape(histogram2_close_intensity_allem(:,1),1,length(histogram2_close_intensity_allem)),Nbin_data);
% nfclose1=nfclose1'/sum(nfclose1);nfclose2=nfclose2'/sum(nfclose2);
% nfclose1(isnan(nfclose1))=0;nfclose2(isnan(nfclose2))=0;
% 
% subplot(2,2,1)    
%     bar(nb,nfall,'r');
%     ylim([0 0.08])
%     xlim([0 60])
%     text(35,0.5*max(ylim),{['\lambda = ',num2str(var_1/mean_1)],['n = ',num2str(n_alldata1)]})
%     xlabel('# hb1')
%     ylabel('Frequency')
% subplot(2,2,3)
%     bar(nb,nfclose1,'y');
%     xlim([0 60])
%     text(35,0.5*max(ylim),{['\lambda = ',num2str(var_1close/mean_1close)],['n = ',num2str(n_closesignal1_allel)]})
%     xlabel('# hb1 close signal')
%     ylabel('Frequency') 
% subplot(2,2,2)
%     bar(nb,nf_2all,'c');
%     ylim([0 0.08])
%     xlim([0 60])
%     text(35,0.5*max(ylim),{['\lambda = ',num2str(var_2/mean_2)],['n = ',num2str(n_alldata2)]})
%     xlabel('# hb7')
%     ylabel('Frequency')
% subplot(2,2,4)
%     bar(nb,nfclose2,'g');
%     xlim([0 60])
%     text(35,0.5*max(ylim),{['\lambda = ',num2str(var_2close/mean_2close)],['n = ',num2str(n_closesignal2_allel)]})
%     xlabel('# hb7 close signal')
%     ylabel('Frequency') 
% figure
% try
%     corr_signal1=corr(closesignal1_allel_re(:,1),closesignal1_allel_re(:,2));
% catch err
%     corr_signal1=nan;
% end
% try
%     corr_signal2=corr(closesignal2_allel_re(:,1),closesignal2_allel_re(:,2));
% catch err
%     corr_signal2=nan;
% end
% subplot(1,2,1)
%     scatter(closesignal1_allembryo(:,1),closesignal1_allembryo(:,2),14,'filled','k')
%     xlabel('intensity_1')
%     ylabel('intensity_2')
%     title(['Hb1 ','\rho=',num2str(corr_signal1)])
%     axis square
%     axis([0 30 0 30])
%     line([5 5],[0 30],'Color','green','LineStyle','--')
%     line([0 30],[5 5],'Color','green','LineStyle','--')
%     line([15 15],[0 30],'Color','red')
%     line([0 30],[15 15],'Color','red')
%     line([20 20],[0 30],'Color','green','LineStyle','--')
%     line([0 30],[20 20],'Color','green','LineStyle','--')
% subplot(1,2,2)
%     scatter(closesignal2_allembryo(:,1),closesignal2_allembryo(:,2),14,'filled','k')
%     xlabel('intensity_1')
%     ylabel('intensity_2')
%     title(['Hb7 ','\rho=',num2str(corr_signal2)])
%     axis square
%     axis([0 30 0 30])
%     line([5 5],[0 30],'Color','green','LineStyle','--')
%     line([0 30],[5 5],'Color','green','LineStyle','--')
%     line([15 15],[0 30],'Color','green','LineStyle','--')
%     line([0 30],[15 15],'Color','green','LineStyle','--')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%hb17-cds
% p22_allcycle_mean=mean(p22_allcycle);
% p22_allcycle_mshape=reshape(p22_allcycle_mean,3,[])';
% p22_allcycle_error=std0(p22_allcycle);
% p22_allcycle_eshape=reshape(p22_allcycle_error,3,[])';
% subplot(2,2,[1,2])
%     barweb(p22_allcycle_mshape(:,1)',p22_allcycle_eshape(:,1)');
%     ylabel('hold ratio')
%     legend(p22_allcycle_name(2:end,:))
%     title('all cycle')
% subplot(2,2,[3,4])
%     barweb(p22_allcycle_mshape(:,2)',p22_allcycle_eshape(:,2)');
%     ylabel('cat ratio')
% saveas(gcf,[['Z:\shihe zhang\S3\embryo_data_hb17\','[0.2 0.5]_0.1_0.02__0107\'],'all-cycle.fig']);
% saveas(gcf,[['Z:\shihe zhang\S3\embryo_data_hb17\','[0.2 0.5]_0.1_0.02__0107\'],'all-cycle.png']);
% for cycle=11:13
%     eval(['p22_allcycle=p22_cycle',num2str(cycle),';'])
%     eval(['p22_allcycle_name=p22_cycle',num2str(cycle),'_name;'])
%     p22_allcycle_mean=mean(p22_allcycle);
%     p22_allcycle_mshape=reshape(p22_allcycle_mean,3,[])';
%     p22_allcycle_error=std0(p22_allcycle);
%     p22_allcycle_eshape=reshape(p22_allcycle_error,3,[])';
%     subplot(2,2,[1,2])
%         barweb(p22_allcycle_mshape(:,1)',p22_allcycle_eshape(:,1)');
%         ylabel('hold ratio')
%         legend(p22_allcycle_name(2:end,:))
%         title(['cycle ',num2str(cycle)])
%     subplot(2,2,[3,4])
%         barweb(p22_allcycle_mshape(:,2)',p22_allcycle_eshape(:,2)');
%         ylabel('cat ratio')
%     saveas(gcf,[['Z:\shihe zhang\S3\embryo_data_hb17\','[0.2 0.5]_0.1_0.02__0107\'],'all-cycle',num2str(cycle),'.fig']);
%     saveas(gcf,[['Z:\shihe zhang\S3\embryo_data_hb17\','[0.2 0.5]_0.1_0.02__0107\'],'all-cycle',num2str(cycle),'.png']);
% end   
%% store combine-data
protein_concentration_el_all={protein_concentration_el_allembryo,protein_concentration_el_11,protein_concentration_el_12,protein_concentration_el_13,protein_concentration_el_14};
for iii = 1:length(elmin)
    data_11{iii}=data_11{iii}./mean(data_11{iii}).*mean(data_11_noreset{iii});
    data_12{iii}=data_12{iii}./mean(data_12{iii}).*mean(data_12_noreset{iii});
    data_13{iii}=data_13{iii}./mean(data_13{iii}).*mean(data_13_noreset{iii});
    data_14{iii}=data_14{iii}./mean(data_14{iii}).*mean(data_14_noreset{iii});
    data_all{iii}=data_all{iii}./mean(data_all{iii}).*mean(data_all_noreset{iii});
end
save([out_folder,'\',['data_all','.mat']],'data_11','data_11_noreset','data_12','data_12_noreset','data_13','data_13_noreset','data_14','data_14_noreset','data_all','data_all_noreset','el','protein_concentration_el_all','EmbryoName_all')
%% combine-data analysis
data_cycle={data_all,data_11,data_12,data_13,data_14};
data_cycle_noreset={data_all_noreset,data_11_noreset,data_12_noreset,data_13_noreset,data_14_noreset};
% ModuleFit_ini(out_folder,data_cycle_noreset,el,protein_concentration_el_allembryo,on_off,simu,model,enrich,11);
%% Treatment

ModuleFit_Kr(in0f,out_folder,data_cycle_noreset,el,protein_concentration_el_all,on_off,simu,model,enrich);

%% send email to me when this procedure has finished
% mailTome('MATLAB Reminder','DataAnalysis has finished')