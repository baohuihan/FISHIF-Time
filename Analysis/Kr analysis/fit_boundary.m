%% fit bundary for Kr, Bcd, Hb
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
 load([bhh_path,sprintf('%d',num(j,1)),'\fit\',char(txt(j,3)),'_',protein_add1,'_intensity_new3.mat']);%oreR new3
 load([bhh_path,sprintf('%d',num(j,1)),'\fit\',char(txt(j,3)),'_',protein_add2,'_intensity_new3.mat']);%oreR new3
    else
       load([bhh_path,char(txt(j,2)),'\fit\',char(txt(j,3)),'_',RNA_add,'_intensity.mat']);
       load([bhh_path,char(txt(j,2)),'\fit\',char(txt(j,3)),'_',protein_add1,'_intensity.mat']);
       load([bhh_path,char(txt(j,2)),'\fit\',char(txt(j,3)),'_',protein_add2,'_intensity.mat']);
    end
    
    
    figure
    %% Bcd fit
 fit_initial=[max(nu_out2(:,2)),0.2,0];
 fit_higher=[max(nu_out2(:,2))*1.1,1,10];
 fit_lower=[0,0,0];
 ft = fittype(@(a1,b1,c1,x) a1*exp(-x/b1)+c1 );
 fitresult= fit( nu_out2(:,1), nu_out2(:,2), ft, 'StartPoint', fit_initial, 'Upper', fit_higher, 'Lower', fit_lower,'Exclude',nu_out2(:,1)<0.1|nu_out2(:,1)>0.75);
 fitout_Bcd(j,1) = fitresult.a1;
 fitout_Bcd(j,2) = fitresult.b1;
 fitout_Bcd(j,3) = fitresult.c1;
 subplot(3,1,1)
 plot(nu_out2(:,1),nu_out2(:,2),'ro');
 hold on
 xi=0:0.01:1;
 plot(xi,fitresult(xi),'r','DisplayName','fit nuclear profile')
 hold off
 nuBcd(j)=nu_out2(12,2);
  tb(j)=num(j,22);
 %% Hb fit
 if char(txt(j,5))=='H'
 fit_initial=[max(nu_out(:,2)),0.4,0.1,0];
 fit_higher=[max(nu_out(:,2))*1.1,1,0.1,1];
 fit_lower=[0,0.2,0,0];
 ft = fittype(@(a2,b2,c2,d2,x) a2*exp(-(x-b2)/c2)./(exp(-(x-b2)/c2)+1)+d2 );
 fitresult= fit( nu_out(:,1), nu_out(:,2), ft, 'StartPoint', fit_initial, 'Upper', fit_higher, 'Lower', fit_lower,'Exclude',nu_out(:,1)<0.15|nu_out(:,1)>0.85);
 fitout_Hb(j,1) = fitresult.a2;
 fitout_Hb(j,2) = fitresult.b2;
 fitout_Hb(j,3) = fitresult.c2;
 fitout_Hb(j,4) = fitresult.d2;
 subplot(3,1,2)
 plot(nu_out(:,1),nu_out(:,2),'go');
 hold on
 xi=0:0.01:1;
 plot(xi,fitresult(xi),'g','DisplayName','fit nuclear profile')
 hold off
%  kra=roundn(num(j,28),-2);
nuhb(j)=fitresult(0.36);
th(j)=num(j,22);
 else
     fitout_Hb(j,1) = 0;
 fitout_Hb(j,2) = 0;
 fitout_Hb(j,3) = 0;
 fitout_Hb(j,4) = 0;
 
 end
 %% Kr fit
 clear filter
 if num(j,6)>3
 fit_initial=[max(fiout(:,2)),0.3,0.1,0.6,0.1,0];
 fit_higher=[max(fiout(:,2))*1.1,0.5,2,1,2,1];
 fit_lower=[0,0,0,0.5,0,0];
 ft = fittype(@(a3,b3,c3,d3,e3,f3,x) a3*(exp((x-b3)/c3)./(exp((x-b3)/c3)+1)).*(exp(-(x-d3)/e3)./(exp(-(x-d3)/e3)+1))+f3 );
 filter=rmmissing(fiout,1);
 fitresult= fit( filter(:,1), filter(:,2), ft, 'StartPoint', fit_initial, 'Upper', fit_higher, 'Lower', fit_lower,'Exclude',filter(:,2)==0);
 fitout_kr(j,1) = fitresult.a3;
 fitout_kr(j,2) = fitresult.b3;
 fitout_kr(j,3) = fitresult.c3;
 fitout_kr(j,4) = fitresult.d3;
 fitout_kr(j,5) = fitresult.e3;
 fitout_kr(j,6) = fitresult.f3;
 subplot(3,1,3)
 plot(fiout(:,1), fiout(:,2),'bo');
 hold on
 xi=0:0.01:1;
 yi=fitresult(xi);
 plot(xi,fitresult(xi),'b','DisplayName','fit nuclear profile')
 hold off
 else
 fitout_kr(j,1) = 0;
 fitout_kr(j,2) = 0;
 fitout_kr(j,3) = 0;
 fitout_kr(j,4) = 0;
 fitout_kr(j,5) = 0;
 fitout_kr(j,6) = 0;
 end
end
%% bound out
% % load('Z:\kr-enhancer\fit2\co_13.mat')
% % Tstart=co(18);
% % Tend=co(21);
boundkr_a=fitout_kr(:,2);
boundkr_p=fitout_kr(:,4);
boundHb=fitout_Hb(:,2);
boundBcd=fitout_Bcd(:,2);
% 
T=num(:,22);
Thb=T;
Tr=T;

boundHb(boundHb==0)=0;
Thb(boundHb==0)=0;
%% 
boundBcd(boundBcd==0)=0;
T(boundBcd==0)=0;
%%
Tr(boundkr_a==0)=0;
%%
boundkr_a(boundkr_a==0)=[];
boundkr_p(boundkr_p==0)=[];
Tr(Tr==0)=[];
T(T==0)=[];
boundBcd(boundBcd==0)=[];
Thb(Thb==0)=[];
boundHb(boundHb==0)=[];



boundkr_11=[Tr,boundkr_a,boundkr_p];
boundBcd_11=[T,boundBcd];
boundHb_11=[Thb,boundHb];
 save('Z:\kr-enhancer\Results2\boundkr_11.mat','boundkr_11');
 save('Z:\kr-enhancer\Results2\boundBcd_11.mat','boundBcd_11');
 save('Z:\kr-enhancer\Results2\boundHb_11.mat','boundHb_11');
%% for 1xBcd  get protein intensity in boundary
close all
load('Z:\kr-enhancer\fit2\co_13.mat')
nuBcd(tb<co(18)|tb>co(21))=0;
nuhb(th<co(18)|th>co(21))=0;

nuhb(nuhb==0)=[];
nuBcd(nuBcd==0)=[];
figure;
bar([mean(nuBcd) mean(nuhb)],'Facecolor','y')
hold on
errorbar(1,mean(nuBcd),std0(nuBcd),'k*')
hold on

errorbar(2,mean(nuhb),std0(nuhb),'k*')
hold on
plot([1 2],[2.57 7.93],'ro')
%plot([0 3],[2.8 2.8],'r')
%hold on
%plot([0 3],[8.4 8.4],'g')

 
 function y = std0(x)
        y = std(x)/sqrt(length(x));
 end

