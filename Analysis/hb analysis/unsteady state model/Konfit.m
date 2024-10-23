%% analysis for kon Pactive
clear
%close all
load('U:\TimeFit\OutputModel\Ton\CyclePlot\Inherit1-AlpOn-Turn1-TimeBin0.5-Combin3-PonLock1-TimeWindowMingle-PonBoole0-0530-bcd\Results13.mat')
%KGA(:,12)=KGA(:,11);
%KGA(:,13)=KGA(:,11);
%YEL(12)=0.6;
%YEL(13)=0.65;

for i=1:45
    if max(KGA(i,:))<0.01
        continue
    else
    figure
    plot(YEL,KGA(i,:),'go')
    hold on
    %%
co0=[max(KGA(i,:))*0.8,0.4,0.01,0];
lb=[0,0.2,0,0];
ub=[max(KGA(i,:))*1.4,1,0.1,1];
[co,resnorm,residual] = lsqcurvefit(@logistfit,co0,YEL(5:11)',KGA(i,5:11)',lb,ub);
xi=0:0.01:0.6;
 plot(xi,logistfit(co,xi),'g','DisplayName','fit nuclear profile')
 hold off
 A(:,i)=co;
    %% fit
%  fit_initial=[max(KGA(i,:))*0.8,0.4,0.01,0];
%  fit_higher=[max(KGA(i,:))*1.4,1,0.1,1];
%  fit_lower=[0,0.2,0,0];
%  ft = fittype(@(a2,b2,c2,d2,x) a2*exp(-(x-b2)/c2)./(exp(-(x-b2)/c2)+1)+d2 );
%  fitresult= fit( YEL', KGA(i,:)', ft, 'StartPoint', fit_initial, 'Upper', fit_higher, 'Lower', fit_lower, 'Exclude',YEL<=0.35);
%  fitout_Hb(i,1) = fitresult.a2;
%  fitout_Hb(i,2) = fitresult.b2;
%  fitout_Hb(i,3) = fitresult.c2;
%  fitout_Hb(i,4) = fitresult.d2;
%  
%  xi=0:0.01:0.55;
%  plot(xi,fitresult(xi),'g','DisplayName','fit nuclear profile')
%  hold off
    end
end
 close all
 figure;
plot(XTime,KGA(:,9))
figure
plot(XTime,A(1,:))
figure
plot(XTime(5:41),A(2,5:41))
%% timefit
figure;plot(XTime(11:17),KGA(11:17,9))%9:19
co0=[0,0.5,0.1];
lb=[0,0,0];
ub=[1,1.5,1];
[co,resnorm,residual] = lsqcurvefit(@expfit,co0,(XTime(11:17)-XTime(11)),KGA(11:17,9),lb,ub);
%[co,resnorm,residual] = lsqcurvefit(@hillfit,co0,nu2(7:12,12),KGA(27,6:11)',lb,ub);
hold on
xi=XTime(10):5:270;
 plot(xi,expfit(co,xi-XTime(11)),'g','DisplayName','fit nuclear profile')
%%
Bcd_out=mean(nu2(5:12,8:11),2);%10:11 8:11
kout=mean(KGA(18:25,4:11),1);%21:25 18:25
figure;plot(Bcd_out,kout,'o');
%figure;plot(nu2(6:12,12),KGA(27,5:11));
 %%
co0=[0,max(kout),6,5];
lb=[0,0,0,0];
ub=[0,max(kout)*4,10,15];
[co,resnorm,residual] = lsqcurvefit(@hillfit,co0,Bcd_out(2:8),kout(2:8)',lb,ub);
%[co,resnorm,residual] = lsqcurvefit(@hillfit,co0,nu2(7:12,12),KGA(27,6:11)',lb,ub);
hold on
xi=0:1:40;
 plot(xi,hillfit(co,xi),'g','DisplayName','fit nuclear profile')
 %k=co(4).^co(3);
 
color1 = [46,114,188]/255;
color2 = [206,85,255]/255;
color3=[0.1 0.1 0.5];
 figure
 for i=5:10
  pn1=polyfit(time(2:3),nu2(i,2:3),1);
   nu_L1=pn1(1).*XTime(5:9)+pn1(2);
  %  pn2=polyfit(time(15:18),nu2(i,15:18),1);
 %  nu_L2=pn2(1).*XTime(32:37)+pn2(2);
   plot(nu_L1,KGA(5:9,i-1))
   hold on
   pn(i)=pn1(1);
%    plot(nu_L2,KGA(32:37,i-1),'--','Color','g')
%    hold on
 end
 figure;plot(YEL(4:9),pn(5:10)*60)
  figure
 for i=5:10
  pn2=polyfit(time(16:18),nu2(i,16:18),1);
   nu_L2=pn2(1).*XTime(35:37)+pn2(2);
   plot(nu_L2,KGA(35:37,i-1))
   hold on
   pn(i)=pn2(1);
 end
 figure;plot(YEL(4:9),pn(5:10)*60)
 figure;
 for i=6:1:10
 pn1=polyfit(time(2:3),nu2(i,2:3),1);
 pn2=polyfit(XTime(5:10),KGA(5:10,i-1),1);
 
 %p2(i)=KGA(8,i-1)-KGA(5,i-1);
 p1(i)=pn1(1);
 p2(i)=pn2(1)*60;
 scatter(p1(i),p2(i))
 hold on
 end
 figure;
 for i=6:1:10
 pn1=polyfit(time(16:18),nu2(i,16:18),1);
 pn2=polyfit(XTime(33:38),KGA(33:38,i-1),1);
 
 %p2(i)=KGA(8,i-1)-KGA(5,i-1);
 p1(i)=pn1(1);
 p2(i)=pn2(1)*60;
 scatter(p1(i),p2(i))
 hold on
 end
 
 %% time
 ELout=0.2:0.05:0.45;
 kupt=[192.4 192.4 167.4 192.4 167.4 192.4];
 kdownt=[525.1 525.1 525.1 525.1 525.1 525.1];
 Bupt=[141 141 141 141 141 141];
 Bdownt=[533 533 533 533 533 533];
 figure;plot(ELout,kupt,'r','--o')
 hold on;plot(ELout,kdownt,'b','-o')
 hold on;plot(ELout,Bupt,'-o')
 hold on;plot(ELout,Bdownt,'-o')
 
 %% Pactive time
 Tons=[66.8 30.39 -17.5];
 Toffs=[110.7 103.2 61.7];
 Activet=[35.98 99.66 125];
 Inactivet=[27.65 77.46 89.73];
 Windowt=[178.9 224.3 521.1];
 cycle=[11 12 13];
 Maxactive=[0.7515 0.7527 0.7663];
 figure;yyaxis left
 h1=plot(cycle,Tons,'r-o');
 hold on;h2=plot(cycle,Toffs,'r--o');
 hold on;h3=plot(cycle,Activet,'g-o');
 hold on;h4=plot(cycle,Inactivet,'g--o');
 hold on;h5=plot(cycle,Windowt,'b-o');
 ylabel('Time')
 yyaxis right
 hold on;h6=plot(cycle,Maxactive,'k-*');
 ylabel('Probability')
 lgd1=legend([h1,h2],'Ton','Toff','orientation','horizontal','location','north');
 lgd2=legend([h3,h4],'Tup','Tdown','orientation','horizontal','location','north');
 lgd3=legend([h5,h6],'Twindow','Maxactive','orientation','horizontal','location','north');
 
 figure
 x=[11 12 13];
 y=[7 9 13];
 figure;bar(x,y)
 y2=[-1 -1 -1];
 hold on;bar(x,y2)
 yup=CycletimeExtractG(:,1)-60;
 ydown=CycletimeExtractG(:,2)-60;
 ystart=CycletimeExtractG(:,3)-60;
 yend=CycletimeExtractG(:,4)-60;
 hold on; plot(x,yup/60,'ro');
 hold on; plot(x,ydown/60,'bo');
 hold on; plot(x,ystart/60,'yo');
 hold on; plot(x,yend/60,'ko');