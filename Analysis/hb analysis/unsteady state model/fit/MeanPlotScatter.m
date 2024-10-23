function y_element=MeanPlotScatter(x,y,model)
if iscell(y)
    if strcmp(model,'nozero')
        y=cellfun(@(x) Zero2Nan(x),y,'UniformOutput',false);
    end
    y_mean=cellfun(@(x) mean(x,1,'omitnan'),y,'UniformOutput',false);
    y_empty=cell2mat(cellfun(@(x) isempty(x),y,'UniformOutput',false));
    y_mean_ture=cell2mat(y_mean(~y_empty));
    y_element=zeros(size(y_mean,1),size(y_mean_ture,2));
    y_element(~y_empty,:)=y_mean_ture;
else
    y_element=mean(y,1);
end
scatter(x,y_element,75,'filled')
end

function y=Zero2Nan(x)
    x(x==0)=nan;
    y=x;
end