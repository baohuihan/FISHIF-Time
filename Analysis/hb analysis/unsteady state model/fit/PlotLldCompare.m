LLD1=LldG;
LLD2=LldG;
LLDc=[LLD1';LLD2'];
LLDcNorm=LLDc./mean(LLDc,1);
figure
boxplot([LLD1,LLD2],'Notch','on','Labels',{'ActiveIn','TwoState'})
figure
plot(TimeLabelMean,[LLD1';LLD2'])
figure
plot(TimeLabelMean,LLDcNorm)
figure
boxplot(LLDcNorm(:,1:end-1)','Notch','on','Labels',{'ActiveIn','TwoState'})
ylabel('Lost')
title('Compare LOST of Models')