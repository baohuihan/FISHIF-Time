function plotTimeWindows(timeWindows, features, TimeElong)
    % 创建图形
    % figure('Position',[800 600 1200 300]);
    hold on;
    
    % 设置颜色映射
    Maxf=1;
    Minf=0;
    colormap('cool');
    clim([Minf, Maxf]); % 根据特征值范围设置颜色轴
    SizePatch=[[-0, -0, -0.4, -0.4];[0, 0, -0.4, -0.4]];
    % 绘制每个时间窗口
    for i = 1:size(timeWindows, 1)
        % 计算颜色
        % colorValue = (features(i) - Minf) / (Maxf - Minf);
        colorValue=0.1;
        color = cool(256);
        color = color(ceil(colorValue * size(color, 1)), :);
        
        % 绘制矩形
        patch([timeWindows(i, 1), timeWindows(i, 2), timeWindows(i, 2), timeWindows(i, 1)],...
              SizePatch(mod(i,2)+1,:), color, 'EdgeColor', 'r', 'FaceAlpha', 0.5);
    
    end
    patch([TimeElong(1), TimeElong(2), TimeElong(2), TimeElong(1)],...
              [0, 0, 0.4, 0.4], 'm', 'FaceAlpha', 0.3);

    % 添加颜色条
    % colorbar;
    
    % 设置图形属性
    xlabel('Time');
    ylabel('Window');
    title('Time Windows Visualization');
    xlim([min([timeWindows(:);TimeElong(1)]), max([timeWindows(:);TimeElong(2)])]);
    ylim([-1, 1]);
    set(gca,'ytick',[],'yticklabel',[])
    hold off;
end
