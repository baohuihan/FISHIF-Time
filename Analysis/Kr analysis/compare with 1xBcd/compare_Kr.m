%% compare Kr in oreR and 1xBcd
clear 
%close all
[num,txt,raw]=xlsread('Z:\kr-enhancer\number_Kr_13new.xlsx','Bcd1x_de');
protein_add1 = 'Bcd protein';
RNA_add = 'Kr RNA';
for j=2:length(txt(:,1))
    if char(txt(j,5))=='H'
       protein_add2 = 'Hb protein';
    else
        protein_add2 = 'Gt protein';
    end
bhh_path=char(txt(j,1));
    if ismissing(txt(j,2))== 1
 load([bhh_path,sprintf('%d',num(j,1)),'\fit\',char(txt(j,3)),'_',RNA_add,'_intensity_new3.mat']);
 load([bhh_path,sprintf('%d',num(j,1)),'\fit\',char(txt(j,3)),'_',protein_add1,'_intensity_new3.mat']);
 load([bhh_path,sprintf('%d',num(j,1)),'\fit\',char(txt(j,3)),'_',protein_add2,'_intensity_new3.mat']);
    else
       load([bhh_path,char(txt(j,2)),'\fit\',char(txt(j,3)),'_',RNA_add,'_intensity.mat']);
       load([bhh_path,char(txt(j,2)),'\fit\',char(txt(j,3)),'_',protein_add1,'_intensity.mat']);
       load([bhh_path,char(txt(j,2)),'\fit\',char(txt(j,3)),'_',protein_add2,'_intensity.mat']);
    end
    %out3D=[nu_out2(3:18,2)',nu_out(3:18)',fiout(3:18,2)'];
    peakr(j-1)=max(fiout(:,2));
    T(j-1)=num(j,22);
    R(:,j-1)=fiout(:,2);    
    
    P1(:,j-1)=nu_out(:,2);%Hb Gt
    P2(:,j-1)=nu_out2(:,2);%Bcd
    Tr(:,j-1)=num(j,22)*ones(21,1);%cycle13 13
    Tp(:,j-1)=num(j,22)*ones(21,1);%cycle13 13
end
PHb=P1;
THb=Tp;
outP(:,1)=num(:,22);
outP(:,2)=num(:,8);
for j=2:length(txt(:,1))
    if num(j,6)<6
        Tr(:,j-1)=0;
         %THb(:,j-1)=0;
         %Tp(:,j-1)=0;
    end
    if char(txt(j,5))=='H'%&&num(j,6)>2
       continue
    else
        PHb(:,j-1)=0;
       THb(:,j-1)=0;
       outP(j,:)=0;
    end
end

 color_all = [0.329411764705882,0.650980392156863,0.619607843137255;...
    0.745098039215686,0.870588235294118,0.874509803921569;...
    0.894117647058824,0.905882352941177,0.894117647058824;...
    0.0156862745098039,0.282352941176471,0.419607843137255;...
    0.305882352941177,0.541176470588235,0.317647058823529];
color_1=color_all(3,:);
color_mean1=color_all(4,:);
color_mean2=color_all(5,:); 
hold on
%figure;
EL=0:0.05:1; 
RR=R.*double((round(Tr,1) >= 9.5)&(round(Tr,1) < 10.5));
RR(find(isnan(RR)==1))=0;
RR(:,all(RR==0,1))=[];
fi0 = mean(RR,2);  
fi1=std(RR,0,2);

%figure
for i=1:length(RR(1,:))
    %plot(EL,RR(:,i),'Color',color_1,'o');
    scatter(EL,RR(:,i),'filled','MarkerFaceColor',color_mean1,'MarkerFaceAlpha',.4);
    hold on
end
%plot(EL,fi0,'Color',color_mean1)
%scatter(EL,fi0,'filled','MarkerFaceColor',color_mean2)

fiout(:,1)=EL;
fiout(:,2)=fi0;
fit_initial=[max(fiout(:,2)),0.3,0.1,0.6,0.1,0];
 fit_higher=[max(fiout(:,2))*1.1,0.5,2,1,2,1];
 fit_lower=[0,0,0,0.5,0,0];
 ft = fittype(@(a3,b3,c3,d3,e3,f3,x) a3*(exp((x-b3)/c3)./(exp((x-b3)/c3)+1)).*(exp(-(x-d3)/e3)./(exp(-(x-d3)/e3)+1))+f3 );
 filter=rmmissing(fiout,1);
 fitresult= fit( filter(:,1), filter(:,2), ft, 'StartPoint', fit_initial, 'Upper', fit_higher, 'Lower', fit_lower,'Exclude',filter(:,2)==0);
 fita(1)=fitresult.b3;
 fitp(1)=fitresult.d3;
 xi=0:0.01:1;
 yi=fitresult(xi);
 plot(xi,fitresult(xi),'Color',color_mean1,'DisplayName','fit nuclear profile')

