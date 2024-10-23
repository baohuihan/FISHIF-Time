%% spatiotemporal profile of hb 
clear 
[num,txt,raw]=xlsread('Z:\bhh-fish\oreR_wjy\number_hb_13new2.xlsx','all');
for j=2:length(txt)
    
bhh_path=char(txt(j,1));
    if ismissing(txt(j,2))== 1
 load([bhh_path,sprintf('%d',num(j,1)),'\RNAapply\',char(txt(j,3)),'_intensity.mat']);
    else
       load([bhh_path,char(txt(j,2)),'\RNAapply\',char(txt(j,3)),'_intensity.mat']); 
    end
    R(:,j-1)=fiout(:,2);
    T(:,j-1)=num(j,8)*ones(101,1);%cycle13 13
    P1(:,j-1)=max(fiout(1:50,2));
    P2(:,j-1)=max(fiout(80:101,2));
    t(j-1)=num(j,8);
    if num(j,5)>5
     %% Hb fit
 fit_initial=[max(fiout(:,2)),0.4,0.1,0];
 fit_higher=[max(fiout(:,2))*1.1,1,0.1,1];
 fit_lower=[0,0.2,0,0];
 ft = fittype(@(a2,b2,c2,d2,x) a2*exp(-(x-b2)/c2)./(exp(-(x-b2)/c2)+1)+d2 );
 fitresult= fit( fiout(:,1), fiout(:,2), ft, 'StartPoint', fit_initial, 'Upper', fit_higher, 'Lower', fit_lower,'Exclude',fiout(:,1)<0.15|fiout(:,1)>0.8);
 fitout_Hb(j-1,1) = fitresult.a2;
 fitout_Hb(j-1,2) = fitresult.b2;
 fitout_Hb(j-1,3) = fitresult.c2;
 fitout_Hb(j-1,4) = fitresult.d2;
 %figure
 plot(fiout(:,1),fiout(:,2),'go');
 hold on
 xi=0:0.01:1;
 plot(xi,fitresult(xi),'g','DisplayName','fit nuclear profile')
 hold off
  else
 fitout_Hb(j-1,1) = 0;
 fitout_Hb(j-1,2) = 0;
 fitout_Hb(j-1,3) = 0;
 fitout_Hb(j-1,4) = 0;
    end
end
boundHb=fitout_Hb(:,2);

%%
    nucleus_bin = 1:1:14;
    average_radius = 1;
 bin_max = min(nucleus_bin+average_radius,13.5);
 bin_min = max(nucleus_bin-average_radius,0.1);
 dq=jet(length(nucleus_bin));
figure;
 for I_bin = 1:length(nucleus_bin)
        %fi0(:,I_bin) = mean(R.*double((T >= bin_min(I_bin))&(T <= bin_max(I_bin))),2);
        RR=R.*double((T >= bin_min(I_bin))&(T <= bin_max(I_bin)));
        RR(find(isnan(RR)==1))=0;
        RR(:,all(RR==0,1))=[];
        fi0(:,I_bin) = mean(RR,2);
        time(I_bin)=mean(T((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
        max_active1(I_bin)=max(fi0(1:50,I_bin));
        max_active2(I_bin)=max(fi0(51:101,I_bin));
        
%         position1=find(diff(fi0(30:90,I_bin))==min(diff(fi0(30:90,I_bin))));
%         position2=find(diff(fi0(30:90,I_bin))==max(diff(fi0(30:90,I_bin))));
%         
%         bound1(I_bin)=fnout(30+position1,1);
%         bound2(I_bin)=fnout(30+position2,1);
        maxP1(I_bin)=mean(P1((t >= bin_min(I_bin))&(t <= bin_max(I_bin))));
        stdP1(I_bin)=std0(P1((t >= bin_min(I_bin))&(t <= bin_max(I_bin))));
        maxP2(I_bin)=mean(P2((t >= bin_min(I_bin))&(t <= bin_max(I_bin))));
        stdP2(I_bin)=std0(P2((t >= bin_min(I_bin))&(t <= bin_max(I_bin))));
        
        bounda(I_bin)=mean(boundHb((t >= bin_min(I_bin))&(t <= bin_max(I_bin))));
        stda(I_bin)=std0(boundHb((t >= bin_min(I_bin))&(t <= bin_max(I_bin))));
        
  %% subplot out      
        %errorbar(fnout(:,1),fi0(:,I_bin),fi1(:,I_bin));
        %subplot(2,7,I_bin)
        plot(fiout(:,1),fi0(:,I_bin),'color',dq(I_bin,:));
        hold on
        
        %hold off
        legend_str{I_bin} = ['T=', num2str(time(I_bin))];
        %legend(['T=', num2str(time(I_bin))])
        %ylim([0 10])
        xlabel('AP axis (normalized)')
        ylabel('Intensity')
        %text(0.05,55,['T=', num2str(time(I_bin))])
        
%% Frame
%         F=getframe(gcf);
%         I=frame2im(F);
%         [I,map]=rgb2ind(I,256);
%         if I_bin == 1
%             imwrite(I,map,'Z:\bhh-fish\RNAapply\cycle13.gif','gif','Loopcount',inf,'DelayTime',1);
%         else
%             imwrite(I,map,'Z:\bhh-fish\RNAapply\cycle13.gif','gif','WriteMode','append','DelayTime',1);
%         end
 end
 %legend(legend_str);
 figure
 errorbar(time,bounda,stda,'r')
 hold on
 scatter(t,boundHb,'filled','MarkerFaceColor','r','MarkerFaceAlpha',.4); 
 hold on
 figure
% plot(time,max_active1,'r')
%  hold on
%  plot(time,max_active2,'b')
%% teo peak of nc13
errorbar(time,maxP1,stdP1,'r')
 hold on
 scatter(t,P1,'filled','MarkerFaceColor','r','MarkerFaceAlpha',.4); 
 hold on
 errorbar(time,maxP2,stdP2,'b')
 hold on
 scatter(t,P2,'filled','MarkerFaceColor','b','MarkerFaceAlpha',.4); 
 hold on
 xlabel('T/min')
 ylabel('intensity(au)')
 
 
  function y = std0(x)
        y = std(x)/sqrt(length(x));
 end