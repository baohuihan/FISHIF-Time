function [timeBoundaries, featureValues] = createTimelines(timeWindows, features)
% method
Method='Mingle';
% Method='Add';

    % 展平时间窗口并排序
    allTimes = unique(timeWindows(:));
    allTimes = sort(allTimes);

    % 初始化特征值数组
    featureValues = zeros(size(allTimes) - [1 1-size(features,2)]);

    % 遍历每个时间段
    for i = 1:length(allTimes) - 1
        % 找到当前时间段覆盖的所有窗口
        inWindow = timeWindows(:,1) < allTimes(i+1) & timeWindows(:,2) > allTimes(i);
        
        % 计算平均特征值
        if ~any(inWindow)
            inWindow = timeWindows(:,1) < allTimes(i+2) & timeWindows(:,2) > allTimes(i-1);
        end
        if strcmp(Method,'Mingle')
            featureValues(i,:) = mean(features(inWindow,:),1);
        elseif strcmp(Method,'Add')
            featuresInWindow=features(inWindow,:);
            featureValues(i,:) = featuresInWindow(1,:);
        end
    end

    % 返回时间边界和特征值
    timeBoundaries = allTimes.';
end
