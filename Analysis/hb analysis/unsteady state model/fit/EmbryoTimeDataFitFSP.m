%% 
clc;
clear;
close all;
foci_folder = 'Histogram\';
result_foler = 'Results\';
mat_tail = '.mat';
%% 

list_name = 'number_hb_13.xlsx';
[hb_13_num,hb_13_list] = xlsread(list_name,'all');
hb_13_num = hb_13_num(2:end,:); hb_13_list = hb_13_list(2:end,:);
T = hb_13_num(:,6)*60;
Embryo0 = find(T>=300 & T<=350);
Channel=hb_13_num(Embryo0,9);
j=0;
figure
foci_inten_all = [];
for list_I = Embryo0(1:end-1)'
j=j+1;
folder_list = [hb_13_list{list_I,1},hb_13_list{list_I,2},'\'];
load([folder_list,result_foler,hb_13_list{list_I,3},'_new',mat_tail],'foci_RNA_profile');

nmax = 100;
bin = 0:1:nmax;
nd_mRNA = hist(foci_RNA_profile(:,3),bin);
foci_inten_all = [foci_inten_all;foci_RNA_profile(foci_RNA_profile(:,1)>0.2&foci_RNA_profile(:,1)<0.4,3)];
% subplot(5,5,j)
% bar(bin,nd_mRNA/sum(nd_mRNA));
% % title([hb_13_list{list_I,3},',T=',num2str(T(list_I))]);
% title(['T=',num2str(T(list_I))]);
end
nd_mRNA = hist(foci_inten_all,bin);
bar(bin,nd_mRNA/sum(nd_mRNA));