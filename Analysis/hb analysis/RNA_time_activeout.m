%% RNA intensity and active nuclei ratio along EL (bin with EL for each embryo)
clear 
%bhh_path='Z:\bhh-fish\oreR_hb\';
[num,txt,raw]=xlsread('Z:\bhh-fish\oreR_wjy\number_hb_12new.xlsx','all');
%txt=readcell('Z:\bhh-fish\oreR_wjy\number_hb_13new2.xlsx','sheet','bhh2','Range','A:B');
RNAchannel=1;
for j=24:length(txt)
    RNAchannel=num(j,9);
    bhh_path=char(txt(j,1));
    if ismissing(txt(j,11))~= 1
        mkdir([bhh_path,char(txt(j,2)),'\RNAapply']);
        load([char(txt(j,11)),char(txt(j,2)),'\Results\',char(txt(j,3)),'_new.mat']);
    elseif ismissing(txt(j,2))== 1
 mkdir([bhh_path,sprintf('%d',num(j,1)),'\RNAapply']);
 load([bhh_path,sprintf('%d',num(j,1)),'\Results_new\',char(txt(j,3)),'_new.mat']);
    else
       mkdir([bhh_path,char(txt(j,2)),'\RNAapply']);
       load([bhh_path,char(txt(j,2)),'\Results\',char(txt(j,3)),'_new.mat']);
     end
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
        M1(j) = mean(RNAout((EL >= 0.2 & EL <= 0.4),4));%nucleus 4
        %S1(j) = sum(nucleus_RNA_profile((EL >= 0.2 & EL <= 0.4),2));
        M2(j) = mean(RNAout((EL >= 0.6 & EL <= 0.8),4));
        %S2(j) = sum(nucleus_RNA_profile((EL >= 0.6 & EL <= 0.8),2));
        Mhb(j)=max(M1(j),M2(j));
  %% Kr     
        
        %Mkr(j) = mean(RNAout((EL >= 0.4 & EL <= 0.6),4));%nucleus 4
   
    %end
    if M1(j)<M2(j)
        nucleus_distance=1-EL;
    else
        nucleus_distance=EL;
    end
    %% intensity 
    nucleus_bin = 0:0.01:1;%0.025:0.05:0.975;
    average_radius = 0.03;
    bin_max = min(nucleus_bin+average_radius,1);
    bin_min = max(nucleus_bin-average_radius,0);
    fi0 = zeros(size(nucleus_bin));
    fi1 = zeros(size(nucleus_bin));
    for I_bin = 1:length(nucleus_bin)
        fi0(I_bin) = mean(RNAout((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)),4));
        fi1(I_bin) = std(RNAout((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)),4));
    end
    fiout=[nucleus_bin',fi0',fi1'];
    %save([bhh_path,char(txt(j,1)),'\RNAapply\',char(txt(j,2)),'_intensity.mat'],'fiout');
    %% active nuclei
nucleus_bin = 0:0.01:1;%0.025:0.05:0.975;
average_radius = 0.03;
Nfoci_nucleus=RNAout(:,3);
    bin_max = min(nucleus_bin+average_radius,1);
    bin_min = max(nucleus_bin-average_radius,0);
    fn0 = zeros(size(nucleus_bin));
    fn1 = zeros(size(nucleus_bin));
    fn2 = zeros(size(nucleus_bin));
    fn3 = zeros(size(nucleus_bin));
    for I_bin = 1:length(nucleus_bin)
        fn0(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 0) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 0) <= bin_max(I_bin)));
        fn1(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 1) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 1) <= bin_max(I_bin)));
        fn2(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 2) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 2) <= bin_max(I_bin)));
        fn3(I_bin) = sum((nucleus_distance(Nfoci_nucleus  > 2) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus  > 2) <= bin_max(I_bin)));
    end
    %fn0 = hist(nucleus_distance(Nfoci_nucleus == 0), nucleus_bin);
    %fn1 = hist(nucleus_distance(Nfoci_nucleus == 1), nucleus_bin);
    %fn2 = hist(nucleus_distance(Nfoci_nucleus == 2), nucleus_bin);
    %fn3 = hist(nucleus_distance(Nfoci_nucleus > 2), nucleus_bin);
    fn_all = fn0+fn1+fn2+fn3+(fn0+fn1+fn2+fn3 == 0);
    fnout=[nucleus_bin',((1-fn0./fn_all)*100)',(fn1./fn_all*100)',(fn2./fn_all*100)',(fn3./fn_all*100)'];
    %save([bhh_path,char(txt(j,1)),'\RNAapply\',char(txt(j,2)),'_active.mat'],'fnout');
    %save([bhh_path,sprintf('%d',num(j,1)),'\RNAapply\',char(txt(j,1)),'_active.mat'],'fnout');
    %% save
     if ismissing(txt(j,2))== 1
 save([bhh_path,sprintf('%d',num(j,1)),'\RNAapply\',char(txt(j,3)),'_intensity.mat'],'fiout');
 save([bhh_path,sprintf('%d',num(j,1)),'\RNAapply\',char(txt(j,3)),'_active.mat'],'fnout');
    else
      save([bhh_path,char(txt(j,2)),'\RNAapply\',char(txt(j,3)),'_intensity.mat'],'fiout');
 save([bhh_path,char(txt(j,2)),'\RNAapply\',char(txt(j,3)),'_active.mat'],'fnout');
     end

end
% Mout(:,1)=Mhb';
% Mout(:,2)=Mkr';
% plot(Time,M,'o')
% xlabel('T');
% ylabel('Mean Concentration (A.U.)');


