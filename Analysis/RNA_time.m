%% mean nascent mRNA for hb, Kr
clear 
bhh_path='Z:\kr-enhancer\';
[num,txt,raw]=xlsread('Z:\kr-enhancer\number_Kr_13new.xlsx','Bcd1x_de','B:C');
%txt=readcell('Z:\bhh-fish\oreR_wjy\number_hb_13new2.xlsx','sheet','bhh2','Range','A:B');
RNAchannel=1;
for j=1:length(txt)
 
 load(['Z:\kr-enhancer\',sprintf('%d',num(j,1)),'\Results2_new\',char(txt(j,1)),'_new.mat']);
 %EL=cytoplasmic_signal2_profile(:,1);
 if RNAchannel==1
     RNAout=nucleus_RNA_profile;
 else
     RNAout=nucleus_signal2_profile;
 end
 
 EL=RNAout(:,1);
%  if mean(EL)>0.5
%     EL=1-EL;
% end
 %% hb
%         M1(j) = mean(RNAout((EL >= 0.2 & EL <= 0.4),4));%nucleus 4
%         %S1(j) = sum(nucleus_RNA_profile((EL >= 0.2 & EL <= 0.4),2));
%         M2(j) = mean(RNAout((EL >= 0.6 & EL <= 0.8),4));
%         %S2(j) = sum(nucleus_RNA_profile((EL >= 0.6 & EL <= 0.8),2));
%         Mhb(j)=max(M1(j),M2(j));
  %% Kr     
        
        Mkr(j) = mean(RNAout((EL >= 0.35 & EL <= 0.55),4));%nucleus 4 oreR:EL >= 0.4 & EL <= 0.6
        nucleus_distance=EL;
   


end
