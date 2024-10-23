%% compare Kr boundary between strains
%clear
[num,txt,~]=xlsread('Z:\kr-enhancer\number_Kr_13new.xlsx','oreR-all4');
T=num(:,22);
Tr=T;
R=num(:,6);
Bkra=num(:,28);
Bkrp=num(:,29);
for j=2:length(txt(:,1))
    if char(txt(j,5))=='G'
       Phb(j)=0;
       Thb(j)=0;
       
    end
    if R(j)<3
        Tr(j)=0;
        
    end
end

%% kr new range
load('Z:\kr-enhancer\fit2\co_13.mat')
outa=Bkra((Tr >= co(18))&(Tr <= co(21)));
outp=Bkrp((Tr >= co(18))&(Tr <= co(21)));

[num,txt,~]=xlsread('Z:\kr-enhancer\number_Kr_13new.xlsx','Bcd1x_de');
T=num(:,22);
Tr=T;
R=num(:,6);
Bkra=num(:,28);
Bkrp=num(:,29);
for j=2:length(txt(:,1))
    if char(txt(j,5))=='G'
       Phb(j)=0;
       Thb(j)=0;
       
    end
    if R(j)<3
        Tr(j)=0;
        
    end
end

%% kr new range
load('Z:\kr-enhancer\fit2\co_13.mat')
outa_1x=Bkra((Tr >= co(18))&(Tr <= co(21)));
outp_1x=Bkrp((Tr >= co(18))&(Tr <= co(21)));

      

%% 
[h1,p1] = ttest2(outa,outa_1x);
[h2,p2] = ttest2(outp,outp_1x);
%%
figure
color1 = [46,114,188]/255;
color2 = [206,85,255]/255;
color3=[0.1 0.1 0.5];
g1 = repmat({'1'},1,length(outa'));
g2 = repmat({'2'},1,length(outa_1x'));
g3 = repmat({'3'},1,length(outp'));
g4 = repmat({'4'},1,length(outp_1x'));
g_1=[g1,g3];
g_2=[g2,g4];
pos1=[1,4];
pos2=[2,5];
x_1=[outa',outp'];
x_2=[outa_1x',outp_1x'];
box1=boxplot(x_1,g_1,'positions',pos1,'Widths',0.5,'Colors',color1,'Symbol','o','OutlierSize',5);
set(box1,'LineWidth',1);
hold on
box2=boxplot(x_2,g_2,'positions',pos2,'Widths',0.5,'Colors',color2,'Symbol','o','OutlierSize',5);
set(box2,'LineWidth',1);
hold on
set(gca,'XTick', (pos1+pos2)/2, 'XTickLabel', ["Kr-A","Kr-P"],'Xlim',[0 6],'Ylim',[0 1]);

%%
figure
h1=bar(1,mean(outa),0.8,'y');
hold on
h2=bar(2,mean(outa_1x),0.8,'b');
hold on
bar(4,mean(outp),0.8,'y')
hold on
bar(5,mean(outp_1x),0.8,'b')
hold on
legend([h1 h2],'WT','1xbcd');

errorbar(1,mean(outa),std(outa),'k*')
hold on
errorbar(2,mean(outa_1x),std(outa_1x),'k*')
hold on
errorbar(4,mean(outp),std(outp),'k*')
hold on
errorbar(5,mean(outp_1x),std(outp_1x),'k*')
hold on


scatter(ones(length(outa),1),outa,'filled','MarkerFaceColor','r','MarkerFaceAlpha',.4);
hold on
scatter(2*ones(length(outa_1x),1),outa_1x,'filled','MarkerFaceColor','r','MarkerFaceAlpha',.4);
hold on
scatter(4*ones(length(outp),1),outp,'filled','MarkerFaceColor','r','MarkerFaceAlpha',.4);
hold on
scatter(5*ones(length(outp_1x),1),outp_1x,'filled','MarkerFaceColor','r','MarkerFaceAlpha',.4);

set(gca,'XTick', (pos1+pos2)/2, 'XTickLabel', ["Kr-A","Kr-P"],'Xlim',[0 6],'Ylim',[0 1]);
 
 function y = std0(x)
        y = std(x)/sqrt(length(x));
 end