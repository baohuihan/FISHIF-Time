%% time dynamic of peak protein, mean RNA intensity
%close all
clear
[num,txt,raw]=xlsread('Z:\kr-enhancer\number_Kr_11new.xlsx','oreR-all4');
T=num(:,22);
Thb=T;
R=num(:,6);
Phb=num(:,8);
Pbcd=num(:,7);
addt=0;
for j=2:length(txt(:,1))
    if char(txt(j,5))=='G'
       Phb(j)=0;
       Thb(j)=0;
    end
    if char(txt(j,25))=='F'
       T(j)=[];
       R(j)=[];
       Pbcd(j)=[];
       Phb(j)=[];
       Thb(j)=[];
    end
end
Phb(Phb==0)=[];
Thb(Thb==0)=[];
%% number_bin
%  fixR(:,1)=T;
%  fixR(:,2)=Pbcd;
%  fixR(:,3)=R;
%  fixR=rmmissing(fixR);
%  fixR0=sortrows(fixR);
%  fixH(:,1)=Thb;
%  fixH(:,2)=Phb;
%  fixH=rmmissing(fixH);
%  fixH0=sortrows(fixH);
%  for i=1:(length(fixR0(:,1))-3)
%    fiT(i)=mean(fixR0(i:i+3,1));
%    fibcd(i)=mean(fixR0(i:i+3,2));
%    stdbcd(i)=std0(fixR0(i:i+3,2));
%    fikr(i)=mean(fixR0(i:i+3,3));
%    stdkr(i)=std0(fixR0(i:i+3,3));
%  end
%  for i=1:(length(fixH0(:,1))-3)
%    fiThb(i)=mean(fixH0(i:i+3,1));
%    fihb(i)=mean(fixH0(i:i+3,2));
%    stdhb(i)=std0(fixH0(i:i+3,2));
%  end
%  figure;
% subplot(3,1,1);
%  hold on
%  errorbar(fiT,fikr,stdkr,'b')
%  hold on
%  scatter(T,R,'filled','MarkerFaceColor','b','MarkerFaceAlpha',.4); 
%   subplot(3,1,2);
%  errorbar(fiT,fibcd,stdbcd,'r')
%  hold on
%  scatter(T,Pbcd,'filled','MarkerFaceColor','r','MarkerFaceAlpha',.4); 
%  subplot(3,1,3);
%  hold on
%  errorbar(fiThb,fihb,stdhb,'g')
%  hold on
%  scatter(Thb,Phb,'filled','MarkerFaceColor','g','MarkerFaceAlpha',.4); 
%  hold on
%% nucleus _bin
nucleus_bin = 1.5:0.5:7.5;%11:1.5:0.5:7.5 12:0.8:0.5:8.8 13:0.5:0.5:13
average_radius = 1;
 bin_max = min(nucleus_bin+average_radius,7.5);%11:7.5 12:9 13:13.5
 bin_min = max(nucleus_bin-average_radius,0.1);
 dq=jet(14);
%figure;
 for I_bin = 1:length(nucleus_bin)
     %% number
     number(I_bin)=length(T((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
     %% RNA
        RR(I_bin)=mean(R((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
        stdR(I_bin)=std0(R((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
        TT(I_bin)=mean(T((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
%         Rmax(I_bin)=max(R((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
%         Rmin(I_bin)=min(R((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
      %% Bcd
        PPbcd(I_bin)=mean(Pbcd((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
        stdbcd(I_bin)=std0(Pbcd((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
%         Pmaxbcd(I_bin)=max(Pbcd((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
%         Pminbcd(I_bin)=min(Pbcd((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
        %% Hb
        PPhb(I_bin)=mean(Phb((Thb >= bin_min(I_bin))&(Thb <= bin_max(I_bin))));
        stdhb(I_bin)=std0(Phb((Thb >= bin_min(I_bin))&(Thb <= bin_max(I_bin))));
        TThb(I_bin)=mean(Thb((Thb >= bin_min(I_bin))&(Thb <= bin_max(I_bin))));
%         Pmaxhb(I_bin)=max(Phb((Thb >= bin_min(I_bin))&(Thb <= bin_max(I_bin))));
%         Pminhb(I_bin)=min(Phb((Thb >= bin_min(I_bin))&(Thb <= bin_max(I_bin))));
 end
TT=rmmissing(TT);
RR=rmmissing(RR);
stdR=rmmissing(stdR);
PPbcd=rmmissing(PPbcd);
stdbcd=rmmissing(stdbcd);
PPhb=rmmissing(PPhb);
stdhb=rmmissing(stdhb);
TThb=rmmissing(TThb);
TT=TT+addt;
T=T+addt;
Thb=Thb+addt;
TThb=TThb+addt;
%%
% figure
%  subplot(3,1,1);
%  hold on
%  errorbar(TT,RR,stdR,'b')
%  hold on
%  scatter(T,R,'filled','MarkerFaceColor','b','MarkerFaceAlpha',.4); 
% %  ylim([0 20])
% %  xlim([1 9])
% %  %xlabel('T');
% %  ylabel('Mean Intensity(au)');
% %  set(gca,'xtick',1:5:9);
% %  set(gca,'ytick',0:10:20);
% %  title('cycle13')
%  subplot(3,1,2);
%  errorbar(TT,PPbcd,stdbcd,'r')
%  hold on
%  scatter(T,Pbcd,'filled','MarkerFaceColor','r','MarkerFaceAlpha',.4); 
%  subplot(3,1,3);
%  hold on
%  errorbar(TT,PPhb,stdhb,'g')
%  hold on
%  scatter(Thb,Phb,'filled','MarkerFaceColor','g','MarkerFaceAlpha',.4); 
%  hold on
% %  ylim([0 180])
% %  xlim([1 9])
% %  xlabel('T/min');
% %  ylabel('Peak Intensity(nM)');
% %  set(gca,'xtick',1:5:9);
% %  set(gca,'ytick',0:90:180);
 %% hb
 %figure
 subplot(3,1,2)
 color_all = [0.329411764705882,0.650980392156863,0.619607843137255;...
    0.745098039215686,0.870588235294118,0.874509803921569;...
    0.894117647058824,0.905882352941177,0.894117647058824;...
    0.0156862745098039,0.282352941176471,0.419607843137255;...
    0.305882352941177,0.541176470588235,0.317647058823529];
color_between=color_all(3,:);
color_mean=color_all(5,:);
hold on
 Pmaxhb=PPhb+stdhb;
 Pminhb=PPhb-stdhb;
 for i = 1:length(TThb)-1
    x = [TThb(i),TThb(i+1),TThb(i+1),TThb(i)];
    y = [Pmaxhb(i),Pmaxhb(i+1),Pminhb(i+1),Pminhb(i)];
    h1 = fill(x,y,'m');
    set(h1,'Facecolor',color_between,'FaceAlpha',0.5,'EdgeColor','none');
 end
 hold on
 %e=errorbar(TThb,PPhb,stdhb,'g');
 e=plot(TThb,PPhb,'g');
 e.Color = color_mean;
 hold on
 scatter(Thb,Phb,'filled','MarkerFaceColor',color_mean,'MarkerFaceAlpha',.4); 
 hold on
xlabel('T/min');
ylabel('Peak intensity(nM)')
title('oreR-Hb protein')
hold on
% line([6 6],[0 40],'Color','k','LineStyle','--')
% hold on
% line([8 8],[0 40],'Color','k','LineStyle','--')
 %% bcd
% figure
subplot(3,1,1)
 color_all = [0.329411764705882,0.650980392156863,0.619607843137255;...
    0.745098039215686,0.870588235294118,0.874509803921569;...
    0.905882352941177,0.894117647058824,0.894117647058824;...
    0.0156862745098039,0.282352941176471,0.419607843137255;...
    0.541176470588235,0.305882352941177,0.317647058823529];
color_between=color_all(3,:);
color_mean=color_all(5,:);
hold on
 Pmaxbcd=PPbcd+stdbcd;
 Pminbcd=PPbcd-stdbcd;
 for i = 1:length(TThb)-1
    x = [TT(i),TT(i+1),TT(i+1),TT(i)];
    y = [Pmaxbcd(i),Pmaxbcd(i+1),Pminbcd(i+1),Pminbcd(i)];
    h1 = fill(x,y,'m');
    set(h1,'Facecolor',color_between,'FaceAlpha',0.5,'EdgeColor','none');
 end
 hold on
  %e=errorbar(TT,PPbcd,stdbcd,'g');
 e=plot(TT,PPbcd,'g');
 e.Color = color_mean;
 hold on
 scatter(T,Pbcd,'filled','MarkerFaceColor',color_mean,'MarkerFaceAlpha',.4); 
 hold on
xlabel('T/min');
ylabel('Peak intensity(nM)')
title('oreR-Bcd protein')
% hold on
% line([6 6],[0 140],'Color','k','LineStyle','--')
% hold on
% line([8 8],[0 140],'Color','k','LineStyle','--')
 %% Kr
 subplot(3,1,3)
 color_all = [0.329411764705882,0.650980392156863,0.619607843137255;...
    0.745098039215686,0.870588235294118,0.874509803921569;...
    0.894117647058824,0.894117647058824,0.905882352941177;...
    0.0156862745098039,0.282352941176471,0.419607843137255;...
    0.541176470588235,0.305882352941177,0.317647058823529];
color_between=color_all(3,:);
color_mean=color_all(1,:);
 hold on
 Rmax=RR+stdR;
 Rmin=RR-stdR;
 for i = 1:length(TThb)-1
    x = [TT(i),TT(i+1),TT(i+1),TT(i)];
    y = [Rmax(i),Rmax(i+1),Rmin(i+1),Rmin(i)];
    h1 = fill(x,y,'m');
    set(h1,'Facecolor',color_between,'FaceAlpha',0.5,'EdgeColor','none');
 end
 hold on
 %e=errorbar(TT,RR,stdR,'g');
 e=plot(TT,RR,'g');
 e.Color = color_mean;
 hold on
 scatter(T,R,'filled','MarkerFaceColor',color_mean,'MarkerFaceAlpha',.4); 
 hold on
xlabel('T/min');
ylabel('mean intensity(au)')
title('oreR-Kr RNA')
hold on
% line([9 9],[0 25],'Color','k','LineStyle','--')
% hold on
% line([11 11],[0 25],'Color','k','LineStyle','--')


 
 function y = std0(x)
        y = std(x)/sqrt(length(x));
 end