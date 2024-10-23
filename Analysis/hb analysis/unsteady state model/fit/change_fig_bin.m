%% change fig bin
uiopen('Y:\TimeFit\OutputModel\Ton\CyclePlotUse\Inherit1-AlpON-Turn1-TimeBin1-EL[0.2 0.4]bin0.2-DataFilter-ShowInherit--Label[0 0 1]-kfitTime[0 1 1]-KL[0.001 0.001 0.01]-KU[0.02 0.02 0.55]-TestBug3\EL 0.2-0.4\FitAnneal.fig',1);
foir = gcf;
hbin2 = figure;

for i =1:12
    figure(foir)
    ax = nexttile(i);
    figData = get(ax, 'Children');
    h = findobj(figData,'Type','Patch'); 
    ydata=get(h,'YData');
    hb1data=([ydata(2,:)]);
    h = findobj(figData,'Type','Line'); 
    ydata=get(h,'YData');
    hb1fit=([ydata{2}]);
    
    hb1data_bin2=bin1TObin2(hb1data);
    hb1fit_bin2=bin1TObin2(hb1fit);
    
    figure(hbin2)
    subplot(2,6,i)
    b=bar(0:2:2*(length(hb1data_bin2)-1),hb1data_bin2,'hist');
    b.FaceColor = [0.3010 0.7450 0.9330];
    hold on
    plot(0:2:2*(length(hb1fit_bin2)-1),hb1fit_bin2,'r','LineWidth',2)
    xlim([0 65])
    ylim([0 0.08])
    xlabel('RNA Number')
    ylabel('Frequence')
    axis square
    title(ax.Title.String)
end