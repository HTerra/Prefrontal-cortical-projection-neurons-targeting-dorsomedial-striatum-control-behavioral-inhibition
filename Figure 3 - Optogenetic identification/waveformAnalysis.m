function p = waveformAnalysis(All, parameters_OI, neuronIndex, IntIxVAR_ITI, IntIxVAR_SD, PyrIxVAR_ITI, PyrIxVAR_SD, OIPyrIx_VAR_ITI, OIPyrIx_VAR_SD, OIIntIxVAR_ITI, OIIntIxVAR_SD)

%% Plots Figure 3B

figure
%subplot(1,2,1)
for cluster2 = [IntIxVAR_ITI IntIxVAR_SD]
        scatter(All{15,1}(cluster2,5), All{15,1}(cluster2,6),10,'b','filled')
        hold on
%         if All{15,1}(cluster2,5) > 0.4 || All{15,1}(cluster2,5) < 0.1 || All{15,1}(cluster2,6) > 1.2
%             display(cluster2)
%         end
end
for cluster2 = [PyrIxVAR_ITI PyrIxVAR_SD]
    scatter(All{15,1}(cluster2,5), All{15,1}(cluster2,6),10,'r','filled')
    hold on
%         if All{15,1}(cluster2,5) > 0.4 || All{15,1}(cluster2,5) < 0.1 || All{15,1}(cluster2,6) > 1.2
%             display(cluster2)
%         end
end
for cluster2 = [OIPyrIx_VAR_ITI OIPyrIx_VAR_SD]
        scatter(All{15,1}(cluster2,5), All{15,1}(cluster2,6),10,'r','filled')
        hold on
%         if All{15,1}(cluster2,5) > 0.4 || All{15,1}(cluster2,5) < 0.1 || All{15,1}(cluster2,6) > 1.2
%             display(cluster2)
%         end
end
for cluster2 = [OIIntIxVAR_ITI OIIntIxVAR_SD]
        scatter(All{15,1}(cluster2,5), All{15,1}(cluster2,6),10,'b','filled')
        hold on
%         if All{15,1}(cluster2,5) > 0.4 || All{15,1}(cluster2,5) < 0.1 || All{15,1}(cluster2,6) > 1.2
%             display(cluster2)
%         end
end
for cluster2 = [OIIntIxVAR_ITI OIIntIxVAR_SD]
        scatter(All{15,1}(cluster2,5), All{15,1}(cluster2,6),10,'b','filled')
        hold on
%         if All{15,1}(cluster2,5) > 0.4 || All{15,1}(cluster2,5) < 0.1 || All{15,1}(cluster2,6) > 1.2
%             display(cluster2)
%         end
end
for cluster2 = find(All{23, 1}==3)
        scatter(All{15,1}(cluster2,5), All{15,1}(cluster2,6),10,'k','filled')
        hold on
%         if All{15,1}(cluster2,5) > 0.4 || All{15,1}(cluster2,5) < 0.1 || All{15,1}(cluster2,6) > 1.2
%             display(cluster2)
%         end
end
%axis([0 0.5 0 1.2])
xlabel('Half peak width (ms)')
ylabel('Peak-Valley width (ms)')
zlabel('Mean firing Rate (Hz)')


% subplot(1,2,2)
% for cluster2 = [IntIxVAR_ITI IntIxVAR_SD]
%     scatter3(All{1,1}{cluster2,1}(1,4), All{1,1}{cluster2,1}(1,5), All{1,1}{cluster2,1}(1,6),10,'MarkerFaceColor',[0 0 1])
%     hold on
% end
% for cluster2 = [PyrIxVAR_ITI PyrIxVAR_SD]
%     scatter3(All{1,1}{cluster2,1}(1,4), All{1,1}{cluster2,1}(1,5), All{1,1}{cluster2,1}(1,6),10,'MarkerFaceColor',[1 0 0])
%     hold on
% end
n=1;
for cluster2 = [OIPyrIx_VAR_ITI OIPyrIx_VAR_SD]
    lightIxs = find(All{1,1}{cluster2,1}(:,3) < 0.01 & All{1,1}{cluster2,1}(:,15)>= parameters_OI.PWCor & All{1,1}{cluster2,1}(:,2)>= parameters_OI.minLightInt & All{1,1}{cluster2,1}(:,6)>= parameters_OI.ReliabilityTreshold & All{1,1}{cluster2,1}(:,4) < parameters_OI.LatencyTreshold & All{1,1}{cluster2,1}(:,5) < parameters_OI.JitterTreshold & numel(find(All{24,1}{cluster2,1}(:,5) == 1)) > 0);
    %lightIxs = find(All{1,1}{cluster2,1}(:,3) < 0.01 & All{1,1}{cluster2,1}(:,15)>= PWCor & All{1,1}{cluster2,1}(:,2)>= minLightInt & All{1,1}{cluster2,1}(:,6)>= ReliabilityTreshold);
    %stimBlocksPower = find(All{1,1}{cluster2,1}(lightIxs,2) == 100 | All{1,1}{cluster2,1}(lightIxs,2) == 75);
    [~, stimBlocksDur] = max(All{1,1}{cluster2,1}(lightIxs,1));
    [~, stimBlocksInt] = max(All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),2));
%     lightIxs = find(All{1,1}{cluster2,1}(:,3) < 0.01 & All{1,1}{cluster2,1}(:,15)>= PWCor & All{1,1}{cluster2,1}(:,2)>= minLightInt & All{1,1}{cluster2,1}(:,6)>= ReliabilityTreshold);
%     [~, stimBlocks] = max(All{1,1}{cluster2,1}(lightIxs,15));
%     stimBlock = lightIxs(stimBlocks);
    latencyTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),4);
    JitterTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),5);
    ReliabilityTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),6);
    WvCorTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),15);
    SALTPTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),3);
%     scatter3(All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),4), All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),5), All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),6),60,'d','MarkerFaceColor',[1 0 0])
%     hold on
    n=n+1;
end

n=1;
for cluster2 = [PyrIxVAR_ITI PyrIxVAR_SD]
    %lightIxs = find(All{1,1}{cluster2,1}(:,3) < 0.01 & All{1,1}{cluster2,1}(:,15)>= parameters_OI.PWCor & All{1,1}{cluster2,1}(:,2)>= parameters_OI.minLightInt & All{1,1}{cluster2,1}(:,6)>= parameters_OI.ReliabilityTreshold & All{1,1}{cluster2,1}(:,4) < parameters_OI.LatencyTreshold & All{1,1}{cluster2,1}(:,5) < parameters_OI.JitterTreshold & numel(find(All{24,1}{cluster2,1}(:,5) == 1)) > 0);
    [~, stimBlocksDur] = max(All{1,1}{cluster2,1}(:,1));
    [~, stimBlocksInt] = max(All{1,1}{cluster2,1}(stimBlocksDur,2));
    %latencyTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),4);
    %JitterTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),5);
    %ReliabilityTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),6);
    WvCorTEMP_nonOI(n) = All{1,1}{cluster2,1}(stimBlocksDur(stimBlocksInt),15);
    SALTPTEMP_nonOI(n) = All{1,1}{cluster2,1}(stimBlocksDur(stimBlocksInt),3);
%     scatter3(All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),4), All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),5), All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),6),60,'d','MarkerFaceColor',[1 0 0])
%     hold on
    n=n+1;
end

n = 1;
for cluster2 = [OIIntIxVAR_ITI OIIntIxVAR_SD]
    lightIxs = find(All{1,1}{cluster2,1}(:,3) < 0.01 & All{1,1}{cluster2,1}(:,15)>= parameters_OI.PWCor & All{1,1}{cluster2,1}(:,2)>= parameters_OI.minLightInt & All{1,1}{cluster2,1}(:,6)>= parameters_OI.ReliabilityTreshold & All{1,1}{cluster2,1}(:,4) < parameters_OI.LatencyTreshold & All{1,1}{cluster2,1}(:,5) < parameters_OI.JitterTreshold);
    %stimBlocksPower = find(All{1,1}{cluster2,1}(lightIxs,2) == 100 | All{1,1}{cluster2,1}(lightIxs,2) == 75);
    [~, stimBlocksDur] = max(All{1,1}{cluster2,1}(lightIxs,1));
    [~, stimBlocksInt] = max(All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),2));
%     lightIxs = find(All{1,1}{cluster2,1}(:,3) < 0.01 & All{1,1}{cluster2,1}(:,15)>= 0.85);
%     [~, stimBlocks] = max(All{1,1}{cluster2,1}(lightIxs,15));
%     stimBlock = lightIxs(stimBlocks);
    latencyTEMP_INT(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),4);
    JitterTEMP_INT(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),5);
    ReliabilityTEMP_INT(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),6);
%     scatter3(All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),4), All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),5), All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),6),60,'d','MarkerFaceColor',[0 0 1])
    n = n+1;
%     hold on
end
% cluster2 = find(All{23, 1}==3);
% for n = 1:numel(find(All{23, 1}==3))
%     scatter3(All{1,1}{cluster2(n),1}(1,4), All{1,1}{cluster2(n),1}(1,5), All{1,1}{cluster2(n),1}(1,6),10,'MarkerFaceColor',[0 0 0])
%     hold on
% end
% %axis([0 0.5 0 1.2])
% xlabel('Latency (ms)')
% ylabel('Jitter (ms)')
% zlabel('Reliability')


figure
subplot(1,4,1)
hold on
Lat = latencyTEMP;
Lat_Groups = [ones(1,size(latencyTEMP',1))];
boxplot(Lat,Lat_Groups,'Labels',{'Lat'})
scatter(ones(size(Lat)).*(1+(rand(size(Lat))-0.5)/5),Lat,8,'k','filled')
ylabel('Latency (ms)')
subplot(1,4,2)
hold on
Jit = [JitterTEMP];
Jit_Groups = [ones(1,size(JitterTEMP',1))];
ylabel('Jitter (ms)')
boxplot(Jit,Jit_Groups,'Labels',{'Jit'})
scatter(ones(size(Jit)).*(1+(rand(size(Jit))-0.5)/5),Jit,8,'k','filled')
subplot(1,4,3)
hold on
Rel = [ReliabilityTEMP];
Rel_Groups = [ones(1,size(ReliabilityTEMP',1))];
ylabel('Reliability')
boxplot(Rel,Rel_Groups,'Labels',{'Rel'})
scatter(ones(size(Rel)).*(1+(rand(size(Rel))-0.5)/5),Rel,8,'k','filled')
subplot(1,4,4)
hold on
WvCor = [WvCorTEMP];
WvCor_Groups = [ones(1,size(WvCorTEMP',1))];
ylabel('Waveform correlation')
boxplot(WvCor,WvCor_Groups,'Labels',{'Wv Cor'})
scatter(ones(size(WvCor)).*(1+(rand(size(WvCor))-0.5)/5),WvCor,8,'k','filled')

% n=1;
% for cluster2 = [PyrIxVAR_ITI PyrIxVAR_SD OIPyrIx_VAR_ITI OIPyrIx_VAR_SD]
%     lightIxs = find(All{1,1}{cluster2,1}(:,2)>= parameters_OI.minLightInt);
%     [~, stimBlocksDur] = min(All{1,1}{cluster2,1}(lightIxs,3));
%     [~, stimBlocksInt] = max(All{1,1}{cluster2,1}(lightIxs(stimBlocksDur),2));
%     SALTPTEMP(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),3);
%     WvCorTEMP_INT(n) = All{1,1}{cluster2,1}(lightIxs(stimBlocksDur(stimBlocksInt)),15);
%     %hold on
%     n=n+1;
% end
figure
hold on
scatter(SALTPTEMP_nonOI,WvCorTEMP_nonOI,10)
scatter(SALTPTEMP,WvCorTEMP,'filled')
xlabel('SALT p-value')
ylabel('WV Cor')
%set(gca,'XScale','log')

figure
scatter3(latencyTEMP, JitterTEMP,ReliabilityTEMP,'filled','k')
xlabel('Latency (ms)')
ylabel('Jitter (ms)')
zlabel('Reliability')

figure
subplot(1,3,1)
hold on
Lat_Pyr = histcounts(latencyTEMP,0:0.1:10,'Normalization','cdf');
Lat_Int = histcounts(latencyTEMP_INT,0:0.1:10,'Normalization','cdf');
stairs(0.1:0.1:10,Lat_Pyr,'r')
stairs(0.1:0.1:10,Lat_Int,'b')
% histogram(latencyTEMP,0:0.1:10,'FaceColor','r','Normalization','probability')
% histogram(latencyTEMP_INT,0:0.1:10,'FaceColor','b','Normalization','probability')
%line([median(latencyTEMP) median(latencyTEMP)], [0 0.2],'Color','red')
%line([median(latencyTEMP_INT) median(latencyTEMP_INT)], [0 0.2],'Color','blue')
xlabel('Latency (ms)')
ylabel('Cum. Prob.')
subplot(1,3,2)
hold on
Jit_Pyr = histcounts(JitterTEMP,0:0.1:7,'Normalization','cdf');
Jit_Int = histcounts(JitterTEMP_INT,0:0.1:7,'Normalization','cdf');
stairs(0.1:0.1:7,Jit_Pyr,'r')
stairs(0.1:0.1:7,Jit_Int,'b')
% histogram(JitterTEMP,0:0.1:7,'FaceColor','r','Normalization','probability')
% histogram(JitterTEMP_INT,0:0.1:7,'FaceColor','b','Normalization','probability')
% line([median(JitterTEMP) median(JitterTEMP)], [0 0.2],'Color','red')
% line([median(JitterTEMP_INT) median(JitterTEMP_INT)], [0 0.2],'Color','blue')
xlabel('Jitter (ms)')
ylabel('Cum. Prob.')
subplot(1,3,3)
hold on
Rel_Pyr = histcounts(ReliabilityTEMP,0:0.01:1,'Normalization','cdf');
Rel_Int = histcounts(ReliabilityTEMP_INT,0:0.01:1,'Normalization','cdf');
stairs(0.01:0.01:1,Rel_Pyr,'r')
stairs(0.01:0.01:1,Rel_Int,'b')
%histogram(ReliabilityTEMP,0:0.01:1,'FaceColor','r','Normalization','probability')
%histogram(ReliabilityTEMP_INT,0:0.01:1,'FaceColor','b','Normalization','probability')
%line([median(ReliabilityTEMP) median(ReliabilityTEMP)], [0 0.2],'Color','red')
%line([median(ReliabilityTEMP_INT) median(ReliabilityTEMP_INT)], [0 0.2],'Color','blue')
xlabel('Reliability')
ylabel('Cum. Prob.')
legend({'Pyr','Int'})

[h,p.Lat,p.Lat_ks2stat] = kstest2(latencyTEMP,latencyTEMP_INT);
[h,p.Jit,p.Jit_ks2stat] = kstest2(JitterTEMP,JitterTEMP_INT);
[h,p.Rel,p.Rel_ks2stat] = kstest2(ReliabilityTEMP,ReliabilityTEMP_INT);
