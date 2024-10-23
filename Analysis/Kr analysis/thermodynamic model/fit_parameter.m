%% parameter analysis for model fit
load('Z:\kr-enhancer\fit2\co_11.mat');
co_11=co;
load('Z:\kr-enhancer\fit2\co_12.mat');
co_12=co;
load('Z:\kr-enhancer\fit2\co_13.mat');
co_13=co;

load('Z:\kr-enhancer\fit\conf_11.mat');
conf_11=conf;
load('Z:\kr-enhancer\fit\conf_12.mat');
conf_12=conf;
load('Z:\kr-enhancer\fit\conf_13.mat');
conf_13=conf;

figure
%subplot(2,1,1)
%yyaxis left
T=[0 co_11(8) co_11(9) 8+co_12(8) 8+co_12(9) 18+co_13(18) 18+co_13(19) 18+co_13(20) 18+co_13(21) 31];
kon=[0 co_11(2) 0 co_12(2) 0 co_13(2) co_13(7) co_13(12) 0 0];
stairs(T,kon)
hold on
xlabel('T/min')
ylabel('k')
ylim([0 34])
yticks([0 17 34])

% yyaxis right
% % Td=[0 7 8 17 18 31];
% % d=[co_11(7) 0 co_12(7) 0 co_13(17) 0];
% % stairs(Td,d)
% line([0,7],[co_11(7),co_11(7)])
% hold on
% line([8,17],[co_12(7),co_12(7)])
% hold on
% line([18,31],[co_13(17),co_13(17)])
% hold on
% xlabel('T/min')
% ylabel('decay')
% xlim([0 31])
% ylim([0 2])
% yticks([0 1 2])
% xticks([0 7 8 17 18 31])


lowco_11=co_11'-conf_11(:,1);
upco_11=conf_11(:,2)-co_11';
lowco_12=co_12'-conf_12(:,1);
upco_12=conf_12(:,2)-co_12';
lowco_13=co_13'-conf_13(:,1);
upco_13=conf_13(:,2)-co_13';


cbcd=[co_11(3),co_12(3),co_13(3),co_13(8),co_13(13)];
low_cbcd=[lowco_11(3),lowco_12(3),lowco_13(3),lowco_13(8),lowco_13(13)];
up_cbcd=[upco_11(3),upco_12(3),upco_13(3),upco_13(8),upco_13(13)];
chb=[co_11(4),co_12(4),co_13(4),co_13(9),co_13(14)];
low_chb=[lowco_11(4),lowco_12(4),lowco_13(4),lowco_13(9),lowco_13(14)];
up_chb=[upco_11(4),upco_12(4),upco_13(4),upco_13(9),upco_13(14)];

nbcd=[co_11(5),co_12(5),co_13(5),co_13(10),co_13(15)];
low_nbcd=[lowco_11(5),lowco_12(5),lowco_13(5),lowco_13(10),lowco_13(15)];
up_nbcd=[upco_11(5),upco_12(5),upco_13(5),upco_13(10),upco_13(15)];
nhb=[co_11(6),co_12(6),co_13(6),co_13(11),co_13(16)];
low_nhb=[lowco_11(6),lowco_12(6),lowco_13(6),lowco_13(11),lowco_13(16)];
up_nhb=[upco_11(6),upco_12(6),upco_13(6),upco_13(11),upco_13(16)];

d=[co_11(7),co_12(7),co_13(17)];
low_d=[lowco_11(7),lowco_12(7),lowco_13(17)];
up_d=[upco_11(7),upco_12(7),upco_13(17)];

figure
subplot(1,2,1)
bar([mean(cbcd) mean(chb)],'Facecolor','y')
hold on
errorbar(1,mean(cbcd),std(cbcd),'k*')
hold on
errorbar(2,mean(chb),std(chb),'k*')
hold on
scatter(1:2,co_11(3:4),'filled','MarkerFaceColor','r','MarkerFaceAlpha',.4);
hold on
scatter(1:2,co_12(3:4),'filled','MarkerFaceColor','g','MarkerFaceAlpha',.4);
hold on
scatter(1:2,co_13(3:4),'filled','MarkerFaceColor','b','MarkerFaceAlpha',.4);
hold on
scatter(1:2,co_13(8:9),'filled','MarkerFaceColor','c','MarkerFaceAlpha',.4);
hold on
scatter(1:2,co_13(13:14),'filled','MarkerFaceColor','y','MarkerFaceAlpha',.4);
hold on
xticks([1 2])
xlim([0 3])
ylim([0 12])
yticks([0 4 8 12])
ylabel('Coefficient Value(nM)')

subplot(1,2,2)
bar([mean(nbcd) mean(nhb)],'Facecolor','y')
hold on
errorbar(1,mean(nbcd),std(nbcd),'k*')
hold on
errorbar(2,mean(nhb),std(nhb),'k*')
hold on
scatter(1:2,co_11(5:6),'filled','MarkerFaceColor','r','MarkerFaceAlpha',.4);
hold on
scatter(1:2,co_12(5:6),'filled','MarkerFaceColor','g','MarkerFaceAlpha',.4);
hold on
scatter(1:2,co_13(5:6),'filled','MarkerFaceColor','b','MarkerFaceAlpha',.4);
hold on
scatter(1:2,co_13(10:11),'filled','MarkerFaceColor','c','MarkerFaceAlpha',.4);
hold on
scatter(1:2,co_13(15:16),'filled','MarkerFaceColor','y','MarkerFaceAlpha',.4);
hold on
xticks([1 2])
xlim([0 3])
ylim([0 12])
yticks([0 4 8 12])
ylabel('Coefficient Value')
% xbcd=0.5:0.5:2.5;
% for i=1:length(cbcd)
% errorbar(xbcd(i),cbcd(i),low_cbcd(i),up_cbcd(i),'*'), grid
% hold on
% end

% errorbar(0.5:0.5:2.5,cbcd,low_cbcd,up_cbcd,'*'),grid
% hold on
% errorbar(5.5:0.5:7.5,chb,low_chb,up_chb,'*')
% hold on
% errorbar(10.5:0.5:12.5,nbcd,low_nbcd,up_nbcd,'*')
% hold on
% errorbar(15.5:0.5:17.5,nhb,low_nhb,up_nhb,'*')
% hold on
% title('Estimated Regression Coefficients with 90% Confidence Intervals')
% ylabel('Coefficient Value')
% xlabel('Estimated Regression Coefficient \beta_j, j = 1,2,3')
% xticks([1.5 6.5 11.5 16.5])
% xlim([.8 3.2])
