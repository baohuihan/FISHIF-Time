%% dynamic of boundary along time axis
clear
load('Z:\kr-enhancer\Results2\boundBcd_11.mat')
load('Z:\kr-enhancer\Results2\boundHb_11.mat')
load('Z:\kr-enhancer\Results2\boundkr_11.mat')
load('Z:\kr-enhancer\fit2\co_11.mat')
bounda_11=boundkr_11(:,2);
boundp_11=boundkr_11(:,3);
boundbcd_11=boundBcd_11(:,2);
boundhb_11=boundHb_11(:,2);
boundbcd_11_de=boundbcd_11;
boundbcd_11_de(boundBcd_11(:,1)<co(8)|boundBcd_11(:,1)>co(9))=[];
boundhb_11_de=boundhb_11;
boundhb_11_de(boundHb_11(:,1)<co(8)|boundHb_11(:,1)>co(9))=[];
load('Z:\kr-enhancer\Results2\boundBcd_12.mat')
load('Z:\kr-enhancer\Results2\boundHb_12.mat')
load('Z:\kr-enhancer\Results2\boundkr_12.mat')
load('Z:\kr-enhancer\fit2\co_12.mat')
bounda_12=boundkr_12(:,2);
boundp_12=boundkr_12(:,3);
boundbcd_12=boundBcd_12(:,2);
boundhb_12=boundHb_12(:,2);
boundbcd_12_de=boundbcd_12;
boundbcd_12_de(boundBcd_12(:,1)<co(8)|boundBcd_12(:,1)>co(9))=[];
boundhb_12_de=boundhb_12;
boundhb_12_de(boundHb_12(:,1)<co(8)|boundHb_12(:,1)>co(9))=[];
load('Z:\kr-enhancer\Results2\boundBcd_13.mat')
load('Z:\kr-enhancer\Results2\boundHb_13.mat')
load('Z:\kr-enhancer\Results2\boundkr_13.mat')
load('Z:\kr-enhancer\fit2\co_13.mat')
bounda_13=boundkr_13(:,2);
boundp_13=boundkr_13(:,3);
boundbcd_13=boundBcd_13(:,2);
boundhb_13=boundHb_13(:,2);
boundbcd_13_de=boundbcd_13;
boundbcd_13_de(boundBcd_13(:,1)<co(18)|boundBcd_13(:,1)>co(21))=[];
boundhb_13_de=boundhb_13;
boundhb_13_de(boundHb_13(:,1)<co(18)|boundHb_13(:,1)>co(21))=[];
%%
%% 
Tadd=18;
load('Z:\kr-enhancer\fit2\co_13.mat')
Tstart=co(18)+Tadd;
Tend=co(21)+Tadd;

Tr=boundkr_13(:,1);
Thb=boundHb_13(:,1);
T=boundBcd_13(:,1);
boundkr_a=bounda_13;
boundkr_p=boundp_13;
boundBcd=boundbcd_13;
boundHb=boundhb_13;
nucleus_bin = 1:1:13;%11:1:0.5:7 12:4.5:0.5:8
average_radius = 1;
 bin_max = min(nucleus_bin+average_radius,13);
 bin_min = max(nucleus_bin-average_radius,0.1);
 dq=jet(14);
 for I_bin = 1:length(nucleus_bin)
     %% RNA
        bounda(I_bin)=mean(boundkr_a((Tr >= bin_min(I_bin))&(Tr <= bin_max(I_bin))));
        stda(I_bin)=std0(boundkr_a((Tr >= bin_min(I_bin))&(Tr <= bin_max(I_bin))));
        boundp(I_bin)=mean(boundkr_p((Tr >= bin_min(I_bin))&(Tr <= bin_max(I_bin))));
        stdp(I_bin)=std0(boundkr_p((Tr >= bin_min(I_bin))&(Tr <= bin_max(I_bin))));
        TTr(I_bin)=mean(Tr((Tr >= bin_min(I_bin))&(Tr <= bin_max(I_bin))));
      %% Bcd
        PPbcd(I_bin)=mean(boundBcd((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
        stdbcd(I_bin)=std0(boundBcd((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
        TT(I_bin)=mean(T((T >= bin_min(I_bin))&(T <= bin_max(I_bin))));
%         Nubcd(I_bin)=mean(nuBcd((Tr >= bin_min(I_bin))&(Tr <= bin_max(I_bin))));
%         stdnubcd(I_bin)=std0(nuBcd((Tr >= bin_min(I_bin))&(Tr <= bin_max(I_bin))));
        %% Hb
        PPhb(I_bin)=mean(boundHb((Thb >= bin_min(I_bin))&(Thb <= bin_max(I_bin))));
        stdhb(I_bin)=std0(boundHb((Thb >= bin_min(I_bin))&(Thb <= bin_max(I_bin))));
        TThb(I_bin)=mean(Thb((Thb >= bin_min(I_bin))&(Thb <= bin_max(I_bin))));
%         Nuhb(I_bin)=mean(nuhb((Thb >= bin_min(I_bin))&(Thb <= bin_max(I_bin))));
%         stdnuhb(I_bin)=std0(nuhb((Thb >= bin_min(I_bin))&(Thb <= bin_max(I_bin))));
 end
 %figure
 TT=TT+Tadd;
 T=T+Tadd;
 TTr=TTr+Tadd;
 Tr=Tr+Tadd;
 TThb=TThb+Tadd;
 Thb=Thb+Tadd;
 Outa=[TTr',bounda',stda'];
 Outa(Outa(:,3)==0,:)=[];
 Outp=[TTr',boundp',stdp'];
 Outp(Outp(:,3)==0,:)=[];
 Outb=[TT',PPbcd',stdbcd'];
 Outb(Outb(:,3)==0,:)=[];
 Outh=[TThb',PPhb',stdhb'];
 Outh(Outh(:,3)==0,:)=[];
 h1=errorbar(Outa(:,1),Outa(:,2),Outa(:,3),'color',[0 0.4470 0.7410]);
 hold on
 scatter(Tr,boundkr_a,'filled','MarkerFaceColor',[0 0.4470 0.7410],'MarkerFaceAlpha',.2); 
 hold on
 h2= errorbar(Outp(:,1),Outp(:,2),Outp(:,3),'color',[0.3010 0.7450 0.9330]);
 hold on
 scatter(Tr,boundkr_p,'filled','MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerFaceAlpha',.2); 
 hold on
 h3=errorbar(Outb(:,1),Outb(:,2),Outb(:,3),'color',[0.6350 0.0780 0.1840]);
 hold on
 scatter(T,boundBcd,'filled','MarkerFaceColor',[0.6350 0.0780 0.1840],'MarkerFaceAlpha',.2); 
 hold on
 h4=errorbar(Outh(:,1),Outh(:,2),Outh(:,3),'color',[0.4660 0.6740 0.1880]);
 hold on
 scatter(Thb,boundHb,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerFaceAlpha',.2); 
 hold on
line([Tstart Tstart],[0 1]);
hold on
line([Tend Tend],[0 1]);
%% all cycle

figure
errorbar([1:3],[mean(bounda_11),mean(bounda_12),mean(bounda_13)],[std(bounda_11),std(bounda_12),std(bounda_13)],'color',[0 0.4470 0.7410]);
hold on
errorbar([1:3],[mean(boundp_11),mean(boundp_12),mean(boundp_13)],[std(boundp_11),std(boundp_12),std(boundp_13)],'color',[0.3010 0.7450 0.9330]);
hold on
errorbar([1:3],[mean(boundbcd_11),mean(boundbcd_12),mean(boundbcd_13)],[std(boundbcd_11),std(boundbcd_12),std(boundbcd_13)],'color',[0.6350 0.0780 0.1840])
hold on
errorbar([1:3],[mean(boundhb_11),mean(boundhb_12),mean(boundhb_13)],[std(boundhb_11),std(boundhb_12),std(boundhb_13)],'color',[0.4660 0.6740 0.1880])
hold on
xlim([0 4])
ylim([0 1])
xticks([1 2 3])
yticks([0 0.5 1])
xlabel('cycle')
ylabel('AP axis(normalized)')

figure
h1=errorbar([1:3],[mean(bounda_11),mean(bounda_12),mean(bounda_13)],[std(bounda_11),std(bounda_12),std(bounda_13)],'color',[0 0.4470 0.7410]);
hold on
h2=errorbar([1:3],[mean(boundp_11),mean(boundp_12),mean(boundp_13)],[std(boundp_11),std(boundp_12),std(boundp_13)],'color',[0.3010 0.7450 0.9330]);
hold on
h3=errorbar([1:3],[mean(boundbcd_11_de),mean(boundbcd_12_de),mean(boundbcd_13_de)],[std(boundbcd_11_de),std(boundbcd_12_de),std(boundbcd_13_de)],'color',[0.6350 0.0780 0.1840]);
hold on
h4=errorbar([1:3],[mean(boundhb_11_de),mean(boundhb_12_de),mean(boundhb_13_de)],[std(boundhb_11_de),std(boundhb_12_de),std(boundhb_13_de)],'color',[0.4660 0.6740 0.1880]);
hold on
xlim([0 4])
ylim([0 1])
xticks([1 2 3])
yticks([0 0.5 1])
xlabel('cycle')
ylabel('AP axis(normalized)')
% figure
% g1 = repmat({'11'},1,length(bounda_11));
% g2 = repmat({'12'},1,length(bounda_12));
% g3 = repmat({'13'},1,length(bounda_13));
% g=[g1,g2,g3];
% x=[bounda_11,bounda_12,bounda_13];
% boxplot(x,g)

 function y = std0(x)
        y = std(x)/sqrt(length(x));
 end