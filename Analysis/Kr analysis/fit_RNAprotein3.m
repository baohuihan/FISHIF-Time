%% Spatial profile for each embryo with Kr RNA, Bcd protein, Hb protein (bin)
%% also get peak intenisty of protein and mean intensity of Kr
 
clear all
close all
% bhh_path='Z:\kr-enhancer\20210501\'; 
% [num,txt,raw]=xlsread([bhh_path,'stacks\matchlist.xls']);
% num1=txt(:,1);
% num2=split(num1,"/");
% num3=num2(:,1);
bhh_path='Z:\kr-enhancer\';
[num,txt,raw]=xlsread('Z:\kr-enhancer\number_Kr_13new.xlsx','Hb-oreR-add');

figure_tail = '.fig';
protein_add1 = 'Bcd protein';
RNA_add = 'Kr RNA';


%subtp = 'background subtraction'; 
Lmin = 0;
Lnbin = 0.05;
Lmax = 1;
Lnwindow = 0.05;
nu_L = Lmin:Lnbin:Lmax;
for j=2%2:length(txt(:,1))
    if char(txt(j,5))=='H'
       protein_add2 = 'Hb protein';
    else
        protein_add2 = 'Gt protein';
    end
    
    mkdir([bhh_path,sprintf('%d',num(j,1)),'\fit']);
    %% Bcd protein
     load([bhh_path,sprintf('%d',num(j,1)),'\Results2_new\',char(txt(j,3)),'_new.mat']);
     %% emmask imerode
    em_mask1=imerode(em_mask,strel('disk',200));
    load([bhh_path,sprintf('%d',num(j,1)),'\masks\',char(txt(j,3)),'_new\mask.mat']);
    for I_layer = 1:size(mask_stack,3)
        em_mask1(:,:,I_layer)=imerode(em_mask,strel('disk',200));
    end
    mask_stack1=mask_stack.*uint16(em_mask1);
    delete_maskstack=mask_stack-mask_stack1;
    delete_mask=double(unique(delete_maskstack(:,:,:)));
    for i=2:length(delete_mask)
        nucleus_protein_profile_ab(delete_mask(i),:)=0;
    end
     nucleus_protein_profile_ab(find(nucleus_protein_profile_ab(:,1)==0),:)=[];
%%   
for Lcenter = 1:length(nu_L)
    nu_map = (nucleus_protein_profile_ab(:,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_protein_profile_ab(:,1) <= nu_L(Lcenter)+Lnwindow);
    nu_mean(Lcenter) = mean(nucleus_protein_profile_ab(nu_map,5)*1e9);
    nu_std(Lcenter) = std0(nucleus_protein_profile_ab(nu_map,5)*1e9);
end
figure(j);
%  if isempty(subtp)
%         sub_nu = 0;
%     else
%         sub_nu = min(nucleus_protein_profile(:,2));
%  end
nu_out2=[nu_L',nu_mean',nu_std'];
save([bhh_path,sprintf('%d',num(j,1)),'\fit\',char(txt(j,3)),'_',protein_add1,'_intensity_new3.mat'],'nu_out2');
 subplot(3,1,1)
 plot(nucleus_protein_profile_ab(:,1),nucleus_protein_profile_ab(:,5)*1e9,'ro','DisplayName',[protein_add1,'Nucleus concentration']);
 hold on
 errorbar(nu_L,nu_mean,nu_std,'k','DisplayName',[protein_add1,'Averaged nuclear profile']);
 hold on;
 
 [MaxIbcd(j),z]=max(nu_mean);
 [MinIbcd(j),~]=min(nu_mean);
 Ibcd2(j) = nu_mean(2); 
 Ibcd3(j) = nu_mean(3); 
%  if z < 3
%      Type(j)=1;
%  else
%      Type(j)=0;
%  end
 %% Hb protein
  load([bhh_path,sprintf('%d',num(j,1)),'\Results_new\',char(txt(j,3)),'_new.mat']);
  %% imerode
  for i=2:length(delete_mask)
        nucleus_protein_profile_ab(delete_mask(i),:)=0;
  end
     nucleus_protein_profile_ab(find(nucleus_protein_profile_ab(:,1)==0),:)=[];
     %%
for Lcenter = 1:length(nu_L)
    nu_map = (nucleus_protein_profile_ab(:,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_protein_profile_ab(:,1) <= nu_L(Lcenter)+Lnwindow);
    nu_mean(Lcenter) = mean(nucleus_protein_profile_ab(nu_map,5)*1e9);
    nu_std(Lcenter) = std0(nucleus_protein_profile_ab(nu_map,5)*1e9);
end
%  if isempty(subtp)
%         sub_nu = 0;
%     else
%         sub_nu = min(nucleus_protein_profile(:,2));
%  end
nu_out=[nu_L',nu_mean',nu_std'];
save([bhh_path,sprintf('%d',num(j,1)),'\fit\',char(txt(j,3)),'_',protein_add2,'_intensity_new3.mat'],'nu_out');

 subplot(3,1,2)
 plot(nucleus_protein_profile_ab(:,1),nucleus_protein_profile_ab(:,5)*1e9,'go','DisplayName',[protein_add2,'Nucleus concentration']);
 hold on
 errorbar(nu_L,nu_mean,nu_std,'k','DisplayName',[protein_add2,'Averaged nuclear profile']);
 hold on;
 
 [MaxIhb(j),~]=max(nu_mean);
 [MinIhb(j),~]=min(nu_mean);
 Ihb2(j) = nu_mean(2); 
 Ihb3(j) = nu_mean(3); 
 
%  title(['concentration profile: ',char(num3(j))]);
%  xlabel('AP axis (normalized)');
 ylabel('Concentration (A.U.)');
 %saveas(1,[bhh_path,'fit\',char(num3(j)),protein_add ,figure_tail]);
 %% foci RNA intensity
 %figure(2);
 load([bhh_path,sprintf('%d',num(j,1)),'\Results2_new\',char(txt(j,3)),'_new.mat']);
 %% imerode
  for i=2:length(delete_mask)
        nucleus_RNA_profile(delete_mask(i),:)=0;
  end
     nucleus_RNA_profile(find(nucleus_RNA_profile(:,1)==0),:)=[];
     %%
     rate=1;
%      if j>=30&&j<66
%          rate=2;
%      else
%          rate=1;
%      end
     
     Mkr(j) = rate.*mean(nucleus_RNA_profile((nucleus_RNA_profile(:,1) >= 0.4 & nucleus_RNA_profile(:,1) <= 0.6),4));
 for Lcenter = 1:length(nu_L)
    nu_map = (nucleus_RNA_profile(:,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_RNA_profile(:,1) <= nu_L(Lcenter)+Lnwindow);
    nu_meanKr(Lcenter) = mean(rate.*nucleus_RNA_profile(nu_map,4));
    nu_stdKr(Lcenter) = std0(rate.*nucleus_RNA_profile(nu_map,4));
  end  
    fiout=[nu_L',nu_meanKr',nu_stdKr'];
     
         
    save([bhh_path,sprintf('%d',num(j,1)),'\fit\',char(txt(j,3)),'_',RNA_add,'_intensity_new3.mat'],'fiout');
    
    [MaxIKr(j),~]=max(nu_meanKr);
    subplot(3,1,3)
 plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,4),'bo','DisplayName',[RNA_add,'Nucleus foci intensity']);
 hold on
 errorbar(nu_L,nu_meanKr,nu_stdKr,'k','DisplayName',[RNA_add,'Mean nucleus foci intensity']);
 %hold on;
 %title(['concentration profile: ',char(num3(j))]);
 xlabel('AP axis (normalized)');
 %ylabel('Concentration (A.U.)');

 saveas(j,[bhh_path,sprintf('%d',num(j,1)),'\fit\',char(txt(j,3)),'_profile_new3' ,figure_tail]);
end


 function y = std0(x)
        y = std(x)/sqrt(length(x));
 end

 
 