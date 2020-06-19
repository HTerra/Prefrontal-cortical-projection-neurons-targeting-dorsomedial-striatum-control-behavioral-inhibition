load DMSprematureSD
load DMSaccuracySD
load DMSomissionsSD

figure
for i = 1:3
    SalIx = find(prematureTestSD(:,2)==1 & prematureTestSD(:,3)==i);
    CNO5 = find(prematureTestSD(:,2)==2 & prematureTestSD(:,3)==i);
    CNO10 = find(prematureTestSD(:,2)==3 & prematureTestSD(:,3)==i);
    subplot(1,3,i)
    hold on
    plot([prematureTestSD(SalIx)'; prematureTestSD(CNO5)'; prematureTestSD(CNO10)'],'k')
    boxplot([prematureTestSD(SalIx) prematureTestSD(CNO5) prematureTestSD(CNO10)],'Labels',{'Sal','5','10'})
    % scatter(ones(size(GemOm(ai,1))).*(1+(rand(size(GemOm(ai,1)))-0.5)/5),GemOm(ai,1),20,'r','filled')
    % scatter(ones(size(GemOm(ai,2))).*(2+(rand(size(GemOm(ai,2)))-0.5)/5),GemOm(ai,2),20,'r','filled')
    % scatter(ones(size(GemOm(ai,3))).*(3+(rand(size(GemOm(ai,3)))-0.5)/5),GemOm(ai,3),20,'r','filled')
    ylabel('Premature')
    xlabel('Dosage');
    ylim([0 50])
end

figure
for i = 1:3
    SalIx = find(accuracyTestSD(:,2)==1 & accuracyTestSD(:,3)==i);
    CNO5 = find(accuracyTestSD(:,2)==2 & accuracyTestSD(:,3)==i);
    CNO10 = find(accuracyTestSD(:,2)==3 & accuracyTestSD(:,3)==i);
    subplot(1,3,i)
    hold on
    plot([accuracyTestSD(SalIx)'; accuracyTestSD(CNO5)'; accuracyTestSD(CNO10)'],'k')
    boxplot([accuracyTestSD(SalIx) accuracyTestSD(CNO5) accuracyTestSD(CNO10)],'Labels',{'Sal','5','10'})
    % scatter(ones(size(GemOm(ai,1))).*(1+(rand(size(GemOm(ai,1)))-0.5)/5),GemOm(ai,1),20,'r','filled')
    % scatter(ones(size(GemOm(ai,2))).*(2+(rand(size(GemOm(ai,2)))-0.5)/5),GemOm(ai,2),20,'r','filled')
    % scatter(ones(size(GemOm(ai,3))).*(3+(rand(size(GemOm(ai,3)))-0.5)/5),GemOm(ai,3),20,'r','filled')
    ylabel('Accuracy')
    xlabel('Dosage');
    ylim([40 100])
end

figure
for i = 1:3
    SalIx = find(omissionsTestSD(:,2)==1 & omissionsTestSD(:,3)==i);
    CNO5 = find(omissionsTestSD(:,2)==2 & omissionsTestSD(:,3)==i);
    CNO10 = find(omissionsTestSD(:,2)==3 & omissionsTestSD(:,3)==i);
    subplot(1,3,i)
    hold on
    plot([omissionsTestSD(SalIx)'; omissionsTestSD(CNO5)'; omissionsTestSD(CNO10)'],'k')
    boxplot([omissionsTestSD(SalIx) omissionsTestSD(CNO5) omissionsTestSD(CNO10)],'Labels',{'Sal','5','10'})
    % scatter(ones(size(GemOm(ai,1))).*(1+(rand(size(GemOm(ai,1)))-0.5)/5),GemOm(ai,1),20,'r','filled')
    % scatter(ones(size(GemOm(ai,2))).*(2+(rand(size(GemOm(ai,2)))-0.5)/5),GemOm(ai,2),20,'r','filled')
    % scatter(ones(size(GemOm(ai,3))).*(3+(rand(size(GemOm(ai,3)))-0.5)/5),GemOm(ai,3),20,'r','filled')
    ylabel('Omissions')
    xlabel('Dosage');
    ylim([0 60])
end