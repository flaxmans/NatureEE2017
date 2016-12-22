function make_SummaryOfClineDistributions()

% function for making a .mat data file summarizing cline with data from the
% sle model (single locus expectations)

mvals = [0.01, 0.02, 0.05, 0.1]; % four different levels of migration

nm = numel(mvals);

singleLocusQuants = 0.01:0.01:0.99; % data summary is about quantiles of the samples
nq = numel(singleLocusQuants);

medianVals = zeros(nm,1);
meanVals = medianVals;

clineSampleSummary = zeros(nm, nq);

for i = 1:nm
    
    filename = ['ClineWidthSamples_m' num2str(mvals(i)) '.txt'];
    
    if exist(filename, 'file')
        
        clineWidthData = importdata(filename);
        
        alleleFreqDiffs = clineWidthData(:,2);
        
        clineSampleSummary(i,:) = quantile(alleleFreqDiffs, singleLocusQuants);
        medianVals(i) = median(alleleFreqDiffs);
        meanVals(i) = mean(alleleFreqDiffs);
        
    else
        error(['File ' filename ' not found!']);
    end
    
end

save('ws_ClineSampleSummaries.mat','mvals','singleLocusQuants','medianVals','meanVals','clineSampleSummary');
    
    

