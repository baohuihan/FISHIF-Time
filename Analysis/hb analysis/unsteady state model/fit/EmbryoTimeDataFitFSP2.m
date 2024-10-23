%% 
clc;
clear;
close all;
foci_folder = 'Histogram\';
result_folder = 'Results\';
masks_folder = 'masks\';
mask_name = 'mask';
mat_tail = '.mat';
%% 

list_name = 'number_hb_13.xlsx';
[hb_13_num,hb_13_list] = xlsread(list_name,'all');
hb_13_num = hb_13_num(2:end,:); hb_13_list = hb_13_list(2:end,:);
T = hb_13_num(:,6)*60;
Embryo0 = find(T>=300 & T<=350);
o=0;
figure
foci_inten_all = [];
for list_I = Embryo0(2:end)'
o=o+1;
folder_list = [hb_13_list{list_I,1},hb_13_list{list_I,2},'\'];
load([folder_list,result_folder,hb_13_list{list_I,3},'_new',mat_tail],'foci_RNA_profile','nucleus_RNA_profile');
% load([folder_list,masks_folder,hb_13_list{list_I,3},'_new\',mask_name,mat_tail],'mask_stack');
% nucleus_cen = regionprops(mask_stack,'Centroid');
% nucleus_num = max(mask_stack(:));
% EL_info = get_EL(em_mask);
% dis_embryo=get_mRNA_embryo_dis(foci_spot,EL_info);
foci_el_hist = hist(foci_RNA_profile(:,1),[0:0.1:1]);
Itrue = mean(foci_el_hist(1:5)) >= mean(foci_el_hist(6:10));

nucleusNum=1:size(nucleus_RNA_profile,1);
nmax = 120;
bin = 0:1:nmax;
if  Itrue
   foci_RNA_el = foci_RNA_profile(foci_RNA_profile(:,1)>=0.2&foci_RNA_profile(:,1)<=0.4,:);
   nucleus_RNA_el = nucleus_RNA_profile(nucleus_RNA_profile(:,1)>=0.2&nucleus_RNA_profile(:,1)<=0.4,:);
   NucleusELIndex=nucleusNum(nucleus_RNA_profile(:,1)>=0.2&nucleus_RNA_profile(:,1)<=0.4);
   Lia=ismember(foci_RNA_profile(:,2),NucleusELIndex);
   foci_RNA_el2=foci_RNA_profile(Lia,:);
else
    
     foci_RNA_el = foci_RNA_profile(foci_RNA_profile(:,1)>=0.6&foci_RNA_profile(:,1)<=0.8,:);
     nucleus_RNA_el = nucleus_RNA_profile(nucleus_RNA_profile(:,1)>=0.6&nucleus_RNA_profile(:,1)<=0.8,:);
end

[nucleus_RNA_profile0,ind_foci] = foci_info_old(nucleus_RNA_profile,foci_RNA_profile);

reg_I = (nucleus_RNA_profile(:,1) >= 0.2) & (nucleus_RNA_profile(:,1) <= 0.4);

    reg_I0 = reg_I(ind_foci);
    pI0 = reg_I0;
    f_ob = nucleus_RNA_profile0(pI0,4);
    % foci_nucleus = hist(nucleus_RNA_el(:,2),[1:1:size(nucleus_RNA_el,1)]);
    nd_mRNA = hist(foci_RNA_el2(:,3),bin);
    foci_nucleus0 = length(find(nucleus_RNA_el(:,3)==0));
foci_nucleus1 = length(find(nucleus_RNA_el(:,3)==1));
foci_nucleus3 = length(find(nucleus_RNA_el(:,3)==3));%+length(find(nucleus_RNA_el(:,3)==2));
nd_mRNA(1) = 2*foci_nucleus0+foci_nucleus1+foci_nucleus3;
r1max = 120;
foci_inten_all = [foci_inten_all;foci_RNA_el2(:,3);ones(nd_mRNA(1),1)];
kon =2;koff =8;ktx =15;
tel = 1;tel2 = 0;
kbd = 35;kdc = 0.5;
pc=1;f_RNAP1 = 1;
k_on = 1:4;
k_off = 3:7;
k_tx = 6:12;
k_bd = 30:35;
k_dc = 0.1:0.1:0.3;
lld_All = zeros(length(k_on)*length(k_off)*length(k_tx)*length(k_bd)*length(k_dc),1);
n=0;
for m = 1:length(k_dc)
    kdc = k_dc(m);
    for l = 1:length(k_bd)
        kbd = k_bd(l);
        for k = 1:length(k_tx)
            ktx = k_tx(k);
            for j = 1:length(k_off)
                koff = k_off(j);
                for i = 1:length(k_on)
                    n = n+1;
                    kon = k_on(i);
[~,Pr1,PPr1] = FSP_NSX_RNAP1_twocopy([kon koff 0 ktx kbd kdc pc],2,f_RNAP1,r1max);
[lld,P_joint,P_m,P_r] = FSP_NSX_joint_r_m_twocopy([kon koff 0 ktx kbd kdc pc],2,foci_RNA_el(:,3),PPr1,r1max,r1max,tel,tel2);
lld_All(n) =lld; 
                end
            end
        end
    end
end


best_fit = find(lld_All==min(lld_All(:)));
[i,j,k,l,m]=ind2sub([length(k_on),length(k_off),length(k_tx),length(k_bd),length(k_dc)],best_fit);
kon = k_on(i);koff = k_off(j);ktx = k_tx(k);kbd = k_bd(l);kdc = k_dc(m);
[~,Pr1,PPr1] = FSP_NSX_RNAP1_twocopy([kon koff 0 ktx kbd kdc pc],2,f_RNAP1,r1max);
kfit{o} = [kon koff 0 ktx kbd kdc pc];
[lld,P_joint,P_m,P_r] = FSP_NSX_joint_r_m_twocopy([kon koff 0 ktx kbd kdc pc],2,foci_RNA_el(:,3),PPr1,r1max,r1max,tel,tel2);

subplot(3,5,o)
bar(bin,nd_mRNA/sum(nd_mRNA));
hold on
plot([0:1:r1max],P_m,LineWidth=2);
% title([hb_13_list{list_I,3},',T=',num2str(T(list_I))]);
legend(['Experiment,T=',num2str(T(list_I))],['FSP,lld=',num2str(lld)])
ylabel('Probability');
xlabel('# of mRNA')
title(['k_{ON}=',num2str(kon),',k_{OFF}=',num2str(koff),',k_{TX} =',num2str(ktx),',k_{DC} =',num2str(kdc),',k_{BD}=',num2str(kbd),',P_{c}=',num2str(pc),',t_{el}=',num2str(tel)])

end
figure;
subplot(1,2,1)
nd_mRNA = hist(foci_inten_all,bin);
bar(bin,nd_mRNA/sum(nd_mRNA));
subplot(1,2,2)
nd_mRNA = hist(f_ob,bin);
bar(bin,nd_mRNA/sum(nd_mRNA));