function make_OSM_figure_PhiDistributions()

load('ws_ClineSampleSummaries.mat','mvals','singleLocusQuants','medianVals','meanVals','clineSampleSummary'); % summary of single locus stats
%clineMvals = mvals;

load('ws_NatureEEboxData_202971-203170.mat');

focalQuantile = 0.25;
useCol = find(singleLocusQuants <= focalQuantile, 1, 'last');

primContactRows = (avals <= 0);

primPhiAtGWC = phiAtGWC_singleLocus(primContactRows, useCol+1);

sdmove = mvals(primContactRows);

close all;

figure;
mlevels = unique(sdmove);
ncols = numel(mlevels);
for i = 1:ncols
    subplot(1,ncols,i);
    workingData = primPhiAtGWC(sdmove == mlevels(i));
    workingData = workingData(~isnan(workingData));
    hist(workingData);
    set(gca, 'FontSize', 12, 'FontName', 'Arial');
    xlabel('\phi at GWC');
    ylabel('frequency');
    title(['{\itm} = ' num2str(mlevels(i))]);
    xlim([0,60]);
    ylim([0,13]);
    xl = get(gca,'XLim');
    yl = get(gca,'YLim');
    mm = round(mean(workingData), 2);
    vv = round(var(workingData), 2);
    text(0.3*xl(2), 0.9*yl(2), ['mean = ' num2str(mm)]);
    text(0.3*xl(2), 0.8*yl(2), ['var. = ' num2str(vv)]);
end
