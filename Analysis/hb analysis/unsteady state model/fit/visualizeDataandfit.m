function visualizeDataandfit(Hanneal, Nbin, Nbar, Pdtb, TimeLabel, T_i)
    % Visualize data and fitting results
    figure(Hanneal);
    nexttile;
    h = bar(Nbin, Nbar, 'hist');
    h.FaceColor = [0.3010 0.7450 0.9330];
    hold on;
    plot(0:70, Pdtb(:,1:71), 'LineWidth', 2);
    ylim([0 0.05]);xlim([-0.5 70])
    xlabel('RNA Number');
    ylabel('Frequency');
    axis square;
    title([num2str(TimeLabel(T_i)), 's']);
end
