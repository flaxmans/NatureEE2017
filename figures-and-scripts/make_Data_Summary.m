function make_Data_Summary(indexes, datadir, quantileLimits, rezip, filename, gatherData, pThreshVec)

% Kruuk et al summed coupling coefficient = (L - 1) s / r
% Barton 1983 coupling coefficient = theta = s / r

load('ws_ClineSampleSummaries.mat','mvals','singleLocusQuants','medianVals','meanVals','clineSampleSummary'); % summary of single locus stats
clineMvals = mvals;


if nargin < 4
    rezip = 0;
end 
if nargin < 3
    quantileLimits = 0.01:0.01:0.99; % total percent of extremes included
end
if nargin < 2
    datadir = '/Volumes/4TB_USB_SG2/';
    %datadir = '~/newneutral/';
end
if nargin < 1
    indexes = 202971:203170;
end
if nargin < 5
    filename = ['ws_NatureEEboxData_' num2str(indexes(1)) '-' num2str(indexes(numel(indexes))) '.mat'];
end
if nargin < 6
    gatherData = 1;
end
if nargin < 7
    pThreshVec = 0.01:0.01:0.49;
end

if ~exist(datadir,'dir')
    datadir = './';
    indexes = 203071;
    filename = ['ws_NatureEEboxData_' num2str(indexes(1)) '.mat'];
end

if exist(filename, 'file') && ~gatherData
    load(filename);
else
    orig = pwd;
    cd(datadir);

    disp('Making parameter summary ...');
    oldn = numel(indexes);
    [mvals, svals, avals, mutdists, constantSvals, indexes] = makeParameterSummary(indexes);
    disp('Parameter summary completed.');
    n = numel(indexes);
    disp(['Data present for ' num2str(n) ' of ' num2str(oldn) ' directories']);

    reachedGWC = zeros(n,1);
    GWCtimes = reachedGWC;
    GWCatCntct = reachedGWC;
    
    phiAtGWC_MAF = zeros(n, numel(pThreshVec));
    phiAt2ndCntct_MAF = phiAtGWC_MAF;
    
    phiAtGWC_quant = zeros(n, numel(quantileLimits));
    phiAt2ndCntct_quant = phiAtGWC_quant;
    
    phiAtGWC_singleLocus = zeros(n, numel(singleLocusQuants));
    phiAt2ndCntct_singleLocus = phiAtGWC_singleLocus;
    barton83thetaGWC = phiAtGWC_singleLocus;
    barton83theta2ndCntct = phiAtGWC_singleLocus;
    

    for i = 1:n
        mywd = ['Run' num2str(indexes(i))];
        disp(['Working on ' mywd '...']);
        cd(mywd);

        [reachedGWC(i), GWCtimes(i), GWCatCntct(i), afdata, aftsAFall] = calculateGWCtime(avals(i), mutdists(i), constantSvals(i), mvals(i), rezip);

        [phiAtGWC_MAF(i,:), phiAt2ndCntct_MAF(i,:), phiAtGWC_quant(i,:), phiAt2ndCntct_quant(i,:), phiAtGWC_singleLocus(i,:), phiAt2ndCntct_singleLocus(i,:), barton83thetaGWC(i,:), barton83theta2ndCntct(i,:)] = calculatePhi(afdata, aftsAFall, GWCtimes(i), reachedGWC(i), avals(i), quantileLimits, pThreshVec, singleLocusQuants, clineMvals, clineSampleSummary);
        
        cd ..;
        
        %disp(['GWCtime = ' num2str(GWCtimes(i)) ', SD_MOVE = ' num2str(mvals(i)) ', phi = ' num2str(phiAtGWC_singleLocus(i,50))])
        
    end
    cd(orig);
    gatherData = 0;
    save(filename);
end

close all;

% %np = numel(pThreshVals);
% makeXYplots(mvals, svals, avals, phiAtGWC(:,3), 'Duration of allopatry', '\phi at GWC');
% 
% makeXYplots(mvals, svals, avals, phiAt2ndCntct(:,3), 'Duration of allopatry', '\phi at 2nd cntct', 'bx');
% 
% makeXYplots(mvals, svals, GWCatCntct, phiAt2ndCntct(:,3), 'Stayed diverged at contact?', '\phi at 2nd cntct', 'kv');
% 
% makeXYplotsProps(mvals, svals, avals, GWCatCntct, '2nd contact time', 'proportion congealed', 'kv', 0);

cd(orig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mvals, svals, avals, mutdists, constantSvals, indexes] = makeParameterSummary(indexes)
n = numel(indexes);
% make sure all directories of data are here:
for i = n:-1:1
    wi = ['Run' num2str(indexes(i))];
    if ~exist(wi, 'dir')
        indexes(i) = [];
        disp(['Warning:  ' wi ' does not exist. Removing from list ...']);
    end
end
% now use them:
n = numel(indexes);
mvals = zeros(n,1);
svals = mvals;
avals = mvals;
mutdists = mvals;
constantSvals = mvals;
for i = 1:n
    wi = ['Run' num2str(indexes(i))];
    cd(wi);
    disp(['Gathering parameters from ' wi '...']);
    if exist('parameters.m.bz2', 'file')
        system('bunzip2 parameters.m.bz2');
    elseif ~exist('parameters.m', 'file')
        error(['parameters.m file missing from ' wi]);
    end
    
    run('./parameters.m');
    svals(i) = MEAN_S;
    mvals(i) = SD_MOVE;
    avals(i) = END_PERIOD_ALLOPATRY;
    if avals(i) < 0
        avals(i) = 0;
    end
    mutdists(i) = MUTATION_DISTRIBUTION;
    constantSvals(i) = DEME0_CONSTANT_S;
    
    clear('parameters');
    cd ..;
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [reachedGWC, GWCtime, GWCatCntct, afSelected, aftsAFall] = calculateGWCtime(alloTime, mutdist, constantSval, mval, rezip)

if nargin < 5
    rezip = 0;
end
if exist('FSTtimeSeries.txt.bz2', 'file')
    system('bunzip2 FSTtimeSeries.txt.bz2');
end
fst = importdata('FSTtimeSeries.txt');
if exist('AlleleFreqTimeSeries.txt.bz2', 'file')
    system('bunzip2 AlleleFreqTimeSeries.txt.bz2');
end
afts = importdata('AlleleFreqTimeSeries.txt');
if rezip
    system('bzip2 FSTtimeSeries.txt');
    system('bzip2 AlleleFreqTimeSeries.txt');
end
allFSTtimes = fst(:,1); % time sampling points
allFSTaf = fst(:,4); % allele frequencies (global)
allFSTsvals = fst(:,5); % selection coefficients
isSelected = fst(:,9);
isSelected = (isSelected == 1);

if size(afts,2) ~= 8
    error('Error in calculateGWCtime: afts not size you think!');
end
aftstimes = afts(:,1);
if ~all(aftstimes == allFSTtimes)
    error('afts data and FST data not concordant in times!');
end
%aftsids = afts(:,2);
aftspatches = afts(:,3:4);
aftsloctype = afts(:,6);
aftsAFdiff = afts(:,4) - afts(:,3);
aftsRev = afts(:,5);
for i = 1:numel(aftsRev)
    if aftsRev(i) == 1
        aftsAFdiff(i) = -aftsAFdiff(i);
    end
end
isSelAF = aftsloctype == 1;
if ~all(isSelAF == isSelected)
    error('afts data and FST data not concordant in selected status of sites!');
end


afAll = [allFSTtimes allFSTaf allFSTsvals isSelected];
afSelected = afAll(isSelected, :);

aftsAFall = [aftstimes allFSTsvals aftsAFdiff aftspatches];
aftsAFall = aftsAFall(isSelAF, :);

if size(afSelected,1) ~= sum(isSelected)
    error('Your data on selected sites are not working as planned');
end

afNumSelCounts = calcAfNumVar(afSelected); % times and number variable
masterTimes = afNumSelCounts(:,1); % the time points at which all types of data exist
emReduced = reduceEm(masterTimes, rezip); % eff. mig. filtered down to time points in masterTimes list

% calculate vector of mean s over time
meansvec = calculateMeanS(afSelected, mutdist, constantSval);

% get GWC time
maxEffMig = calcMaxEffMig(afNumSelCounts, meansvec, mval);

[GWCtime, reachedGWC] = calcGWCtime3(emReduced, maxEffMig, alloTime);
if GWCtime == alloTime;
    GWCatCntct = 1;
else
    GWCatCntct = 0;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function afNumVarCounts = calcAfNumVar(af)
times = af(:,1);
utimes = unique(times);
afnv = zeros(size(utimes));
for i = 1:numel(utimes)
    afnv(i) = sum(times == utimes(i));
end
afNumVarCounts = [utimes afnv];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function emReduced = reduceEm(masterTimes, rezip)
if exist('EffectiveMigrationRates.txt.bz2', 'file')
    system('bunzip2 EffectiveMigrationRates.txt.bz2');
end
em = importdata('EffectiveMigrationRates.txt');
if rezip
    system('bzip2 EffectiveMigrationRates.txt');
end
emtimes = em(:,1);

commonTimes = intersect(masterTimes, emtimes);
nct = numel(commonTimes);
if nct ~= numel(masterTimes)
    error('Fewer common times than expected')
end
emReduced = zeros(nct, size(em,2));
for i = 1:nct
    emrows = (emtimes == commonTimes(i));
    emReduced(i,:) = em(emrows,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function meansvec = calculateMeanS(afSelected, mutdist, constantSval)
if any(afSelected(:,4) ~= 1)
    error('Data not working as expected');
end
times = afSelected(:,1);
utimes = unique(times);
nuts = numel(utimes);
meansvec = ones(nuts,1);
if mutdist == 3
    % constant S vals
    meansvec = constantSval .* meansvec;
else
    svals = afSelected(:,3);
    for i = 1:nuts
        wrs = (times == utimes(i));
        wsvals = svals(wrs);
        meansvec(i) = geomean(wsvals);
        if meansvec(i) <= 0
            error('Geomean wsvals included zeros');
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxEffMig = calcMaxEffMig(afNumVarCounts, meansvec, m)
n = size(afNumVarCounts,1); % num time points
if n ~= numel(meansvec)
    error('Bad vector match')
end
maxEffMig = zeros(n,1);
peqvec = zeros(n,1);
for i = 1:n
    peqvec(i) = oneLocusEq(m, meansvec(i));
end
qvec = 1 - peqvec;
for i = 1:n
    peq = peqvec(i);
    q = qvec(i);
    s = meansvec(i);
    L = (afNumVarCounts(i,2));
    
    randomRes = (1 + (peq * s))^L; % full s contribution if p = 1
    randomImm = (1 + (q * s))^L;
    patchAvgFit = (randomImm * m) + ((1 - m) * randomRes);
    maxEffMig(i) = (randomImm * m) / patchAvgFit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gwcTime, gwcReached] = calcGWCtime3(emReduced, maxEffMig, END_PERIOD_ALLOPATRY)
t = emReduced(:,1);
em = emReduced(:,2);
n = numel(t);
if em(n) >= 0.5 * maxEffMig(n) % still close to random at very end
    gwcTime = NaN;
    gwcReached = 0;
else
    gwcIndex = (find(em >= maxEffMig, 1, 'last')) + 1; % last time eff. mig. greater than independent loci threshold
    if size(gwcIndex,1) < 1
        gwcTime = END_PERIOD_ALLOPATRY;
    else
        gwcTime = t(gwcIndex);
    end
    gwcReached = 1;
end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phiAtGWC_MAF, phiAt2ndCntct_MAF, phiAtGWC_quant, phiAt2ndCntct_quant, phiAtGWC_singleLocus, phiAt2ndCntct_singleLocus, barton83thetaGWC, barton83theta2ndCntct] = calculatePhi(afSelected, aftsAFall, GWCtime, reachedGWC, aval, quantileLimits, pThreshVec, singleLocusQuants, mvals, clineSampleSummary)
% afSelected has the following four columns:
% allFSTtimes allFSTaf allFSTsvals isSelected

% USING MAF:
np = numel(pThreshVec);
phiAtGWC_MAF = zeros(1,np);
phiAt2ndCntct_MAF = phiAtGWC_MAF;

nql = numel(quantileLimits);
phiAtGWC_quant = zeros(1,nql);
phiAt2ndCntct_quant = phiAtGWC_quant;

nsql = numel(singleLocusQuants);
phiAtGWC_singleLocus = zeros(1,nsql);
phiAt2ndCntct_singleLocus = phiAtGWC_singleLocus;
barton83theta2ndCntct = phiAtGWC_singleLocus;
barton83thetaGWC = phiAtGWC_singleLocus;

for i = 1:np
    pThresh = pThreshVec(i);
    rows = ( (abs(afSelected(:,2) - 0.5)) <= pThresh );
    L = sum(rows);
    if L < 1
        phiAtGWC_MAF(i) = NaN;
        phiAt2ndCntct_MAF(i) = NaN;
    else
        useData = afSelected(rows,:);
        
        
        run('./parameters.m');
        
        if aval <= 0
            phiAt2ndCntct_MAF(i) = NaN;
        else
            time = find(useData(:,1) >= aval, 1, 'first');
            time = useData(time,1);
            timeRows = (useData(:,1) == time);
            L = sum(timeRows);
            if L < 1
                phiAt2ndCntct_MAF(i) = NaN;
            elseif L == 1
                phiAt2ndCntct_MAF(i) = 0;
            else
                afTimeData = useData(timeRows,:);
                
                allomeanS = geomean(afTimeData(:,3));
                if allomeanS <= 0
                    error('Error in calculating mean s value');
                end
                
                if GAMETE_PRODUCTION_MODE == 1
                    phiAt2ndCntct_MAF(i) = (L - 1) * allomeanS / 0.5; 
                elseif GAMETE_PRODUCTION_MODE == 0
                    allorecomb = calcRecomb(L, 0.02 * TOTAL_MAP_LENGTH, nCHROMOSOMES);
                    phiAt2ndCntct_MAF(i) = (L - 1) * allomeanS / allorecomb;
                else
                    error('GAMETE_PRODUCTION_MODE not 1 or 2');
                end
            end
        end
        if ~reachedGWC
            phiAtGWC_MAF(i) = NaN;
        elseif GWCtime == aval
            phiAtGWC_MAF(i) = phiAt2ndCntct_MAF(i);
        else
            time = find(useData(:,1) >= GWCtime, 1, 'first');
            time = useData(time,1);
            timeRows = (useData(:,1) == time);
            L = sum(timeRows);
            if L < 1
                phiAtGWC_MAF(i) = NaN;
            elseif L == 1
                phiAtGWC_MAF(i) = 0;
            else
                afTimeData = useData(timeRows,:);
                
                GWCmeanS = geomean(afTimeData(:,3));
                if GWCmeanS <= 0
                    error('Error in calculating mean s value');
                end
                
                if GAMETE_PRODUCTION_MODE == 1
                    phiAtGWC_MAF(i) = (L - 1) * GWCmeanS / 0.5;
                elseif GAMETE_PRODUCTION_MODE == 0
                    recomb = calcRecomb(L, 0.02 * TOTAL_MAP_LENGTH, nCHROMOSOMES);
                    phiAtGWC_MAF(i) = (L - 1) * GWCmeanS / recomb;
                else
                    error('GAMETE_PRODUCTION_MODE not 1 or 2');
                end
            end
        end
    end
end

% USING DATA-INFORMED QUANTILES ON CLINE WIDTHS and SINGLE LOCUS:
% columsn of aftsAFall: aftsAFall = [aftstimes allFSTsvals aftsAFdiff aftspatches];
afDiffs = aftsAFall(:,3);

% SINGLE-LOCUS EXPECTATIONS:
nsql = numel(singleLocusQuants);
% allopatry calcs
if aval < 1
    phiAt2ndCntct_quant(:) = NaN;
    phiAt2ndCntct_singleLocus(:) = NaN;
    barton83theta2ndCntct(:) = NaN;
else
    time = find(aftsAFall(:,1) >= aval, 1, 'first');
    time = aftsAFall(time,1);
    useRows = aftsAFall(:,1) == time;
    useData = aftsAFall(useRows, :);
    useAFdiffs = afDiffs(useRows);
    if numel(useRows) < 1
        error('Your data have no rows!');
    elseif numel(useRows) < 2
        phiAt2ndCntct_quant(:) = NaN;
    else
        for i = 1:nql
            qtile = quantileLimits(i);
            afCutoff1 = quantile(useAFdiffs, qtile);
            countedRows = useAFdiffs >= afCutoff1;
            L = sum(countedRows);
            if L < 1
                phiAt2ndCntct_quant(i) = NaN;
            elseif L == 1
                phiAt2ndCntct_quant(i) = 0;
            else
                gmeans = geomean(useData(countedRows, 2));
                if GAMETE_PRODUCTION_MODE > 0
                    phiAt2ndCntct_quant(i) = (L - 1) * gmeans / 0.5;
                else
                    recomb = calcRecomb(L, 0.02 * TOTAL_MAP_LENGTH, nCHROMOSOMES);
                    phiAt2ndCntct_quant(i) = (L - 1) * gmeans / recomb;
                end
            end
        end
        % shifting to single locus:
        slrow = find(mvals == SD_MOVE);
        for i = 1:nsql
            sllim = clineSampleSummary(slrow, i);
            countedRows = useAFdiffs >= sllim; % narrower than cutoff
            L = sum(countedRows);
            if L < 1
                phiAt2ndCntct_singleLocus(i) = NaN;
                barton83theta2ndCntct(i) = NaN;
            elseif L == 1
                phiAt2ndCntct_singleLocus(i) = 0;
                barton83theta2ndCntct(i) = 0;
            else
                gmeans = geomean(useData(countedRows, 2));
                if GAMETE_PRODUCTION_MODE > 0
                    phiAt2ndCntct_singleLocus(i) = (L - 1) * gmeans / 0.5;
                    barton83theta2ndCntct(i) = gmeans / 0.5;
                else
                    recomb = calcRecomb(L, 0.02 * TOTAL_MAP_LENGTH, nCHROMOSOMES);
                    phiAt2ndCntct_singleLocus(i) = (L - 1) * gmeans / recomb;
                    barton83theta2ndCntct(i) = gmeans / recomb;
                end
            end
            
        end
    end
end
% at time of GWC:
if GWCtime == aval
    phiAtGWC_quant = phiAt2ndCntct_quant;
    phiAtGWC_singleLocus = phiAt2ndCntct_singleLocus;
    barton83thetaGWC = barton83theta2ndCntct;
elseif ~reachedGWC
    phiAtGWC_quant(:) = NaN;
    phiAtGWC_singleLocus(:) = NaN;
    barton83thetaGWC(:) = NaN;
else
    time = find(aftsAFall(:,1) >= GWCtime, 1, 'first');
    time = aftsAFall(time,1);
    useRows = aftsAFall(:,1) == time;
    useData = aftsAFall(useRows, :);
    useAFdiffs = afDiffs(useRows);
    if numel(useRows) < 1
        error('Your data have no rows!');
    elseif numel(useRows) < 2
        phiAtGWC_quant(:) = NaN;
    else
        for i = 1:nql
            qtile = quantileLimits(i);
            afCutoff1 = quantile(useAFdiffs, qtile);
            countedRows = useAFdiffs >= afCutoff1;
            L = sum(countedRows);
            if L < 1
                phiAtGWC_quant(i) = NaN;
            elseif L == 1
                phiAtGWC_quant(i) = 0;
            else
                gmeans = geomean(useData(countedRows, 2));
                if GAMETE_PRODUCTION_MODE > 0
                    phiAtGWC_quant(i) = (L - 1) * gmeans / 0.5;
                else
                    recomb = calcRecomb(L, 0.02 * TOTAL_MAP_LENGTH, nCHROMOSOMES);
                    phiAtGWC_quant(i) = (L - 1) * gmeans / recomb;
                end
            end
        end
        % shifting to single locus:
        slrow = find(mvals == SD_MOVE);
        %disp(slrow);
        
        for i = 1:nsql
%             disp(slrow)
%             disp(i)
%             disp(size(clineSampleSummary))
%             
            sllim = clineSampleSummary(slrow, i);
%             if i == 50
%                 disp('hi sllim')
%                 disp(sllim)
%             end
            
            countedRows = useAFdiffs >= sllim; % narrower than cutoff
            %disp(singleLocusQuants(i))
            %disp(useAFdiffs(countedRows))
            L = sum(countedRows);
            if L < 1
                phiAtGWC_singleLocus(i) = NaN;
                barton83thetaGWC(i) = NaN;
            elseif L == 1
                phiAtGWC_singleLocus(i) = 0;
                barton83thetaGWC(i) = 0;
            else
                gmeans = geomean(useData(countedRows, 2));
                if GAMETE_PRODUCTION_MODE > 0
                    phiAtGWC_singleLocus(i) = (L - 1) * gmeans / 0.5;
                    barton83thetaGWC(i) = gmeans / 0.5;
                else
                    recomb = calcRecomb(L, 0.02 * TOTAL_MAP_LENGTH, nCHROMOSOMES);
                    phiAtGWC_singleLocus(i) = (L - 1) * gmeans / recomb;
                    barton83thetaGWC(i) = gmeans / recomb;
%                     if i == 25
%                         disp(sllim)
%                         disp(GWCtime)
%                         disp(L)
%                         disp(gmeans)
%                         disp(recomb)
%                         disp(phiAtGWC_singleLocus(i))
%                     end
                end
            end
            
        end
    end
end

    





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recomb = calcRecomb(L, M, C)
if L <= C
    recomb = 0.5;
else
    nNeighbors = L - 1; % inter-locus distances
    nNeighborsOnSame = L - C; % # distances between loci on same chromosome
    nNeighborsOnDiff = nNeighbors - nNeighborsOnSame;
    sameChromAvgDist = M / L;
    diffChromDist = 0.5;
    recomb = ((nNeighborsOnSame * sameChromAvgDist) + (nNeighborsOnDiff * diffChromDist)) / nNeighbors;
    if recomb > 0.5
        disp(M)
        disp(C)
        disp(L)
        recomb = 0.5;
    elseif recomb < 0
        error('Negative recombination distance')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eq = oneLocusEq(m, s)
eq = (m * (0.5 + 0.75 * s) - 0.125 * s - 0.125 * sqrt(s^2 - 4 * m * s^2 + 4 * m^2 * (2 + s)^2))/((-0.25 + m) * s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeXYplots(mvals, svals, x, y, labx, laby, sym, useLog)
if nargin < 7
    sym = 'ro';
end
if nargin < 8
    useLog = 1;
end

um = unique(mvals);
num = numel(um);
us = unique(svals);
nus = numel(us);
figure;
count = 0;
for i = 1:num;
    wm = um(i);
    wiis = find(mvals == wm);
    for j = 1:nus;
        ws = us(j);
        wjjs = find(svals == ws);
        wrows = intersect(wiis, wjjs);
        wx = x(wrows);
        wy = y(wrows);
        
        count = count + 1;
        subplot(num, nus, count);
        plot(wx, wy, sym);
        set(gca, 'FontSize', 12, 'FontName', 'Arial');
        if useLog
            set(gca,'YScale','log');
        end
        %ylim([1,max(y)]);
        if j == 1
            ylabel(laby);
        end
        if count > (nus * (num-1)) == 1
            xlabel(labx);
        end
        title(['{\itm} = ' num2str(wm) ', {\its} = ' num2str(ws)]);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeXYplotsProps(mvals, svals, x, y, labx, laby, sym, useLog)


um = unique(mvals);
num = numel(um);
us = unique(svals);
nus = numel(us);
figure;
count = 0;
for i = 1:num;
    wm = um(i);
    wiis = find(mvals == wm);
    for j = 1:nus;
        ws = us(j);
        wjjs = find(svals == ws);
        wrows = intersect(wiis, wjjs);
        wx = x(wrows);
        wy = y(wrows);
        ux = unique(x);
        nux = numel(ux);
        prop = ux;
        for k = 1:nux
            wwx = ux(k);
            wxdata = wy(wx == wwx);
            prop(k) = mean(wxdata);
        end
        
        count = count + 1;
        subplot(num, nus, count);
        plot(ux, prop, sym);
        set(gca, 'FontSize', 12, 'FontName', 'Arial');
        if useLog
            set(gca,'YScale','log');
        end
        %ylim([1,max(y)]);
        if j == 1
            ylabel(laby);
        end
        if count > (nus * (num-1)) == 1
            xlabel(labx);
        end
        title(['{\itm} = ' num2str(wm) ', {\its} = ' num2str(ws)]);
        
    end
end



    
