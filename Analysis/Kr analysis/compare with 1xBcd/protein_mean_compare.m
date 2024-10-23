%% compare protein between oreR and 1xBcd
clear
%close all
[num,txt,raw]=xlsread('Z:\kr-enhancer\number_Kr_13new.xlsx','oreR-all4');
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
    peakp1(j-1)=max(nu_out(:,2));
    T(j-1)=num(j,22);
    R(:,j-1)=fiout(:,2);
    P1(:,j-1)=nu_out(:,2);%Hb Gt
    P2(:,j-1)=nu_out2(:,2);%Bcd
    Tr(:,j-1)=num(j,22)*ones(21,1);%cycle13 13
    Tp(:,j-1)=num(j,22)*ones(21,1);%cycle13 13
end
PHb=P1;
PGt=P1;
THb=Tp;
TGt=Tp;
outP(:,1)=num(:,22);
outP(:,2)=num(:,8);
for j=2:length(txt(:,1))
    if char(txt(j,5))=='H'
       PGt(:,j-1)=0;
       TGt(:,j-1)=0;
    else
        PHb(:,j-1)=0;
       THb(:,j-1)=0;
       outP(j,:)=0;
    end
end
load('Z:\kr-enhancer\fit2\co_13.mat')
outP(outP(:,2)==0,:)=[];
 PPHb=PHb.*double((Tp >= co(18))&(Tp <=co(21)));
PPHb(find(isnan(PPHb)==1))=0;
PPHb(:,all(PPHb==0,1))=[];
nuHb = mean(PPHb,2);

 PP2=P2.*double((Tp >= co(18))&(Tp <= co(21)));
PP2(find(isnan(PP2)==1))=0;
PP2(:,all(PP2==0,1))=[];
nu2 = mean(PP2,2);

EL=0:0.05:1;
figure
plot(EL,nuHb,'g');
hold on
plot(EL,nu2,'r')
xlim([0.2 0.8])
save('Z:\kr-enhancer\Results2\meanBcd_oreR.mat','nu2');
save('Z:\kr-enhancer\Results2\meanHb_oreR.mat','nuHb');
%%
load('Z:\kr-enhancer\Results2\meanBcd_oreR.mat')
load('Z:\kr-enhancer\Results2\meanHb_oreR.mat')
Bcd_oreR=nu2(5:15);
Hb_oreR=nuHb(5:15);
load('Z:\kr-enhancer\Results2\meanBcd_1x.mat')
load('Z:\kr-enhancer\Results2\meanHb_1x.mat')
Bcd_1x=nu2(5:15);
Hb_1x=nuHb(5:15);
figure;h=plot(Bcd_oreR,Bcd_1x,'o');
hold on;g=plot(Hb_oreR,Hb_1x,'o');
x=0:1:18;
y1=x;
y2=0.5*x;
plot(x,y1,'k--');
hold on;
plot(x,y2,'k--')
xlim([0 18])
ylim([0 18])
%%
nu_outhb(:,1)=EL;
nu_outhb(:,2)=nuHb;
fit_initial=[max(nu_outhb(:,2)),0.4,0.1,0];
fit_higher=[max(nu_outhb(:,2))*1.1,1,0.1,1];
fit_lower=[0,0.2,0,0];
ft = fittype(@(a2,b2,c2,d2,x) a2*exp(-(x-b2)/c2)./(exp(-(x-b2)/c2)+1)+d2 );
fitresult= fit( nu_outhb(:,1), nu_outhb(:,2), ft, 'StartPoint', fit_initial, 'Upper', fit_higher, 'Lower', fit_lower,'Exclude',nu_out(:,1)<0.15|nu_out(:,1)>0.85);
cohb(1)=fitresult.a2;
cohb(2)=fitresult.b2;
cohb(3)=fitresult.c2;
cohb(4)=fitresult.d2;

nu_outbcd(:,1)=EL;
nu_outbcd(:,2)=nu2;
fit_initial=[max(nu_outbcd(:,2)),0.2,0];
fit_higher=[max(nu_outbcd(:,2))*1.1,1,10];
fit_lower=[0,0,0];
ft = fittype(@(a1,b1,c1,x) a1*exp(-x/b1)+c1 );
fitresult= fit( nu_out2(:,1), nu_out2(:,2), ft, 'StartPoint', fit_initial, 'Upper', fit_higher, 'Lower', fit_lower,'Exclude',nu_out2(:,1)<0.1|nu_out2(:,1)>0.75);
cobcd(1) = fitresult.a1;
cobcd(2) = fitresult.b1;
cobcd(3) = fitresult.c1;
