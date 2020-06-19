load DMSprematureITI
load DMSaccuracyITI
load DMSomissionsITI

figure
for i = 1:3
    SalIx = find(prematureTestITI(:,2)==1 & prematureTestITI(:,3)==i);
    CNO5 = find(prematureTestITI(:,2)==2 & prematureTestITI(:,3)==i);
    CNO10 = find(prematureTestITI(:,2)==3 & prematureTestITI(:,3)==i);
    subplot(1,3,i)
    hold on
    plot([prematureTestITI(SalIx)'; prematureTestITI(CNO5)'; prematureTestITI(CNO10)'],'k')
    boxplot([prematureTestITI(SalIx) prematureTestITI(CNO5) prematureTestITI(CNO10)],'Labels',{'Sal','5','10'})
    % scatter(ones(size(GemOm(ai,1))).*(1+(rand(size(GemOm(ai,1)))-0.5)/5),GemOm(ai,1),20,'r','filled')
    % scatter(ones(size(GemOm(ai,2))).*(2+(rand(size(GemOm(ai,2)))-0.5)/5),GemOm(ai,2),20,'r','filled')
    % scatter(ones(size(GemOm(ai,3))).*(3+(rand(size(GemOm(ai,3)))-0.5)/5),GemOm(ai,3),20,'r','filled')
    ylabel('Premature')
    xlabel('Dosage');
    ylim([0 50])
end

figure
for i = 1:3
    SalIx = find(accuracyTestITI(:,2)==1 & accuracyTestITI(:,3)==i);
    CNO5 = find(accuracyTestITI(:,2)==2 & accuracyTestITI(:,3)==i);
    CNO10 = find(accuracyTestITI(:,2)==3 & accuracyTestITI(:,3)==i);
    subplot(1,3,i)
    hold on
    plot([accuracyTestITI(SalIx)'; accuracyTestITI(CNO5)'; accuracyTestITI(CNO10)'],'k')
    boxplot([accuracyTestITI(SalIx) accuracyTestITI(CNO5) accuracyTestITI(CNO10)],'Labels',{'Sal','5','10'})
    % scatter(ones(size(GemOm(ai,1))).*(1+(rand(size(GemOm(ai,1)))-0.5)/5),GemOm(ai,1),20,'r','filled')
    % scatter(ones(size(GemOm(ai,2))).*(2+(rand(size(GemOm(ai,2)))-0.5)/5),GemOm(ai,2),20,'r','filled')
    % scatter(ones(size(GemOm(ai,3))).*(3+(rand(size(GemOm(ai,3)))-0.5)/5),GemOm(ai,3),20,'r','filled')
    ylabel('Accuracy')
    xlabel('Dosage');
    ylim([40 100])
end

figure
for i = 1:3
    SalIx = find(omissionsTestITI(:,2)==1 & omissionsTestITI(:,3)==i);
    CNO5 = find(omissionsTestITI(:,2)==2 & omissionsTestITI(:,3)==i);
    CNO10 = find(omissionsTestITI(:,2)==3 & omissionsTestITI(:,3)==i);
    subplot(1,3,i)
    hold on
    plot([omissionsTestITI(SalIx)'; omissionsTestITI(CNO5)'; omissionsTestITI(CNO10)'],'k')
    boxplot([omissionsTestITI(SalIx) omissionsTestITI(CNO5) omissionsTestITI(CNO10)],'Labels',{'Sal','5','10'})
    % scatter(ones(size(GemOm(ai,1))).*(1+(rand(size(GemOm(ai,1)))-0.5)/5),GemOm(ai,1),20,'r','filled')
    % scatter(ones(size(GemOm(ai,2))).*(2+(rand(size(GemOm(ai,2)))-0.5)/5),GemOm(ai,2),20,'r','filled')
    % scatter(ones(size(GemOm(ai,3))).*(3+(rand(size(GemOm(ai,3)))-0.5)/5),GemOm(ai,3),20,'r','filled')
    ylabel('Omissions')
    xlabel('Dosage');
    ylim([0 60])
end