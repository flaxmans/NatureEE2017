function make_OSM_figure_TimeSeriesAndParametric(workingIndexes, gatherDataOverride, rezip, focalQuantile)

% function for making figures of phi, L_e, and m_e over time and for parametric
% plots of phi vs. cline width

load('ws_ClineSampleSummaries.mat','mvals','singleLocusQuants','medianVals','meanVals','clineSampleSummary'); % summary of single locus stats
clineMvals = mvals;

load('ws_NatureEEboxData_202971-203170.mat');
loadedIndexes = indexes;

filename = 'ws_NatureEEboxTimeSeriesPlots.mat';
plotRows = 2;


if nargin < 4
    focalQuantile = 0.25;
end

% graphics directives:
mylw = 1.5; % line width in plots
xticklocs = [10^-7,10^-6,10^-5,0.0001,0.001,0.01,0.1,1,10,100,1000,1E4,1E5];
tickLabs = cell(size(xticklocs));
for i = 1:numel(xticklocs)
    tickLabs{i} = num2str(xticklocs(i));
end
effMigColor = [0 0 0];
deg = sprintf('%c', char(176));
limitEffMig = 1;
lowerLimit = 0.0001;
gwcLineCol = [0.85,0.85,0];
%letters = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U'};
letters = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u'};
clinerefcol = [0.7, 0.7, 0.7];
afdrefcol = [0.9 0.7 0.7];
myms = 7;
simStartFill = [0,1,0];
gwcMark = 'kd';
gwcMarkSize = 8;
gwcFill = [1,1,0];
endAlloSymbol = '>';
endAlloColor = [0 0.8 0];
endAlloSize = 8;
myfs = 12; % font size in plots
lb = 0.8; % bounds on plot limits (graphics)
ub = 1.1;
myfontname = 'Arial';


if nargin < 3
    rezip = 0;
end
if nargin < 2
    gatherDataOverride = 0;
end
if nargin < 1
    % workingIndexes = [202971 203021 203071 203121] + 400; % allopatry 7500
     workingIndexes = [202971 203021 203071 203121]; % primary contact
    %workingIndexes = 203071;
end
plotCols = numel(workingIndexes);

if ~exist(filename, 'file') || gatherDataOverride
    gatherData = 1;
else
    gatherData = 0;
end
n = numel(workingIndexes);

if gatherData
    if exist('/Volumes/4TB_USB_SG2/', 'dir')
        gatherData = 1;
        datadir = '/Volumes/4TB_USB_SG2/';
    elseif exist('/Volumes/SG5TB3/', 'dir')
        datadir = '/Volumes/SG5TB3/';
        gatherData = 1;
    else
        datadir = pwd;
    end
else
    load(filename);
end

orig = pwd;
close all;
figure;
if gatherData == 1
    disp(['Data directory = ' datadir]);
    cd(datadir); 
    for i = 1:n
        % import needed data:
        workingDir = ['Run' num2str(workingIndexes(i))];
        if ~exist(workingDir, 'dir')
            error('You do not have the data directory you need here');
        end
        cd(workingDir);
        if exist('EffectiveMigrationRates.txt.bz2', 'file')
            system('bunzip2 EffectiveMigrationRates.txt.bz2');
        end
        if exist('EffectiveMigrationRates.txt', 'file')
            EffMig = importdata('EffectiveMigrationRates.txt');
        else
            error('Effective migration rate data missing!');
        end
        if exist('FSTtimeSeries.txt.bz2', 'file')
            system('bunzip2 FSTtimeSeries.txt.bz2');
        end
        if exist('FSTtimeSeries.txt', 'file')
            FST = importdata('FSTtimeSeries.txt');
        else
            error('FST data file missing!');
        end
        if exist('AlleleFreqTimeSeries.txt.bz2', 'file')
            system('bunzip2 AlleleFreqTimeSeries.txt.bz2');
        end
        if exist('AlleleFreqTimeSeries.txt', 'file')
            afts = importdata('AlleleFreqTimeSeries.txt');
        else
            error('AlleleFreqTimeSeries.txt missing!');
        end

        if exist('PhiValuesTimeSeries.csv', 'file') && exist('ClineWidthTimeSeries.csv', 'file')  && ~gatherDataOverride
            calcPhi = 0;
            disp(['Reading file ' pwd '/PhiValuesTimeSeries.csv ...']); 
            phiImport = importdata('PhiValuesTimeSeries.csv', ',', 1);
            phiHeaders = phiImport.textdata;
            phiData = phiImport.data;
            aftsImport = importdata('ClineWidthTimeSeries.csv', ',', 1);
            clineWidthHeaders = aftsImport.textdata;
            clineWidthData = aftsImport.data;
            barton83import = importdata('Barton83thetaTimeSeries.csv', ',', 1);
            barton83headers = barton83import.textdata;
            barton83ccData = barton83import.data;
            countedLociImport = importdata('CountedLociTimeSeries.csv', ',', 1);
            Lheaders = countedLociImport.textdata;
            Ldata = countedLociImport.data;
            meanSimport = importdata('MeanStimeSeries.csv', ',', 1);
            meanSheaders = meanSimport.textdata;
            meanSdata = meanSimport.data;
            afdiffimport = importdata('AFdiffTimeSeries.csv', ',', 1);
            afDiffHeaders = afdiffimport.textdata;
            afDiffData = afdiffimport.data;
        else
            calcPhi = 1;
            [phiHeaders, phiData, clineWidthHeaders, clineWidthData, barton83headers, barton83ccData, afDiffHeaders, afDiffData, Lheaders, Ldata, meanSheaders, meanSdata] = calculatePhiValuesAndClineWidths(FST, afts, clineMvals, singleLocusQuants, clineSampleSummary);
        end
        
        % work with data:
        EffMigTimes = EffMig(:,1);
        EffMigRates = EffMig(:,2);
        
        
        % finish and clean up with this directory:
        if calcPhi
            writePhiAndClineWidthTStoFile(phiHeaders, phiData, 'PhiValuesTimeSeries.csv');
            writePhiAndClineWidthTStoFile(clineWidthHeaders, clineWidthData, 'ClineWidthTimeSeries.csv');
            writePhiAndClineWidthTStoFile(barton83headers, barton83ccData, 'Barton83thetaTimeSeries.csv');
            writePhiAndClineWidthTStoFile(afDiffHeaders, afDiffData, 'AFdiffTimeSeries.csv');
            writePhiAndClineWidthTStoFile(Lheaders, Ldata, 'CountedLociTimeSeries.csv');
            writePhiAndClineWidthTStoFile(meanSheaders, meanSdata, 'MeanStimeSeries.csv');
        end
        if rezip
            system('bzip2 EffectiveMigrationRates.txt');
            system('bzip2 FSTtimeSeries.txt');
            system('bzip2 AlleleFreqTimeSeries.txt');
        end
        
        % figures
        subplot(plotRows, plotCols, i);
        
        phiColumn = find(singleLocusQuants <= focalQuantile, 1, 'last');
        if numel(phiColumn) ~= 1
            error('error in finding right quantile');
        end
        
        x1 = EffMigTimes;
        y1 = EffMigRates;
        x2 = phiData(:,1);
        y2 = phiData(:,(phiColumn + 1)); % plus one for first column being times
        x3 = clineWidthData(:,1);
        y3 = clineWidthData(:,(phiColumn + 1));
        x4 = barton83ccData(:,1);
        y4 = barton83ccData(:,(phiColumn + 1));
        if limitEffMig
            [x1,y1,x2,y2,~,y3,~,~] = restrictDataRange(x1,y1,x2,y2,x3,y3,x4,y4,lowerLimit);
        end
        % plot eff mig
        plot(x1,y1, 'Color', effMigColor, 'LineWidth', mylw);
        set(gca,'YScale','log','YTick',xticklocs, 'YTickLabels', tickLabs, 'FontSize', myfs, 'FontName', 'Arial');
        hold on;
        
        % the following is Barton 1983's n_e (my interpretation of it):
        [effNumLoci, ~] = calcEffNumLoci( afDiffData(:,(phiColumn + 1)), meanSdata(:,(phiColumn + 1)) );
        effNumLoci = effNumLoci(1:numel(x2));
        plot(x2, effNumLoci, '--', 'Color', [1 0.5 0], 'LineWidth', mylw);
        
        % y limits:
        maxy = 1.1 * max(max(y2(isfinite(y2))), max(effNumLoci));
        myy = [0.8 * min(y1(y1>0)), maxy];
        ylim(myy);
        
        % plot GWC time reference line:
        gwcIndex = find(loadedIndexes == workingIndexes(i));
        if numel(gwcIndex) ~= 1
            error('run index not properly found in loaded data');
        end
        gwcTime = GWCtimes(gwcIndex);
        plot([gwcTime, gwcTime], myy, '-.','LineWidth',mylw,'Color',gwcLineCol);
        
        xlim([0, max(x1)]);
        run('./parameters.m');
        if i == 1
            h = legend('{\itm}_e', '{\itL}_e', 'GWC time');
            set(h, 'FontSize', myfs, 'FontName', 'Arial');
        end
        h = xlabel('Time (generations)', 'Interpreter', 'LaTeX');
        set(h, 'FontSize', myfs + 1);
        if ( i == 1 )
            h = ylabel('Value of metric', 'Interpreter', 'LaTeX');
            set(h, 'FontSize', myfs + 1);
        end
        [tx,ty] = getTextxy([x1;x2], [y1;y2;effNumLoci]);
        h = text(tx,ty, ['\bf{' letters{i} '}']);    
        set(h, 'FontSize', myfs, 'FontName', myfontname);
        title(['{\itm} = ' num2str(SD_MOVE) ', {\its} = ' num2str(MEAN_S)]);
        
        
        % PARAMETRIC
        % next subplot: parametric cline width vs. phi
        subplot(plotRows, plotCols, i + plotCols);
        maxWidth = getOneLocusEqCline2(meanSdata(:,(phiColumn + 1)));
        minWidth = getFullyCoupledCline(meanSdata(:,(phiColumn + 1)), Ldata(:,(phiColumn + 1)));
        maxWidth = maxWidth(1:numel(y2));
        minWidth = minWidth(1:numel(y2));
        
        
        plot(y2, y3, 'Color', [0.95 0.25 0], 'LineWidth', mylw); % parametric plot with actual cline data
        hold on;
        [xvec, yvec] = sortAFDs(y2, maxWidth);
        plot(xvec, yvec, ':', 'LineWidth', mylw + 0.5, 'Color', clinerefcol);
        [phiSort, minWidthSort] = sortAFDs(y2, minWidth);
        plot(phiSort, minWidthSort, '--', 'LineWidth', mylw-0.5, 'Color', afdrefcol);
        
        startIndex = find((isfinite(y2) + (y2 > 0)) == 2, 1, 'first');
        plot(y2(startIndex), y3(startIndex), 'ko', 'MarkerSize', myms, 'MarkerFaceColor', simStartFill);
        plot(y2(numel(y2)), y3(numel(y2)), 'ks', 'MarkerSize', myms + 1, 'MarkerFaceColor', [1 0 0]);
        set(gca, 'FontSize', myfs, 'FontName', 'Arial', 'YScale', 'Log');
        h = xlabel('$\phi$', 'Interpreter', 'LaTeX');
        set(h, 'FontSize', myfs + 2);
        if i == 1
            h = ylabel('$1/\left(p_2-p_1\right)$', 'Interpreter', 'LaTeX');
            set(h, 'FontSize', myfs + 2);
        end
        limity2 = y2(isfinite(y2));
        limity2 = limity2(limity2 > 0);
        set(gca, 'XScale', 'Log','XTick',xticklocs, ...
            'XLim',[lb*min(limity2), ub*max(y2)], ...
            'XTickLabels', tickLabs, ...
            'YLim', [lb*min(y3), max(10,ub*(max(max(y3(isfinite(y3))),max(yvec(isfinite(yvec))))))], ...
            'YTick', xticklocs, 'YTickLabels', tickLabs);
        gwcIndex = find(phiData(:,1) <= gwcTime, 1, 'last');
        plot(y2(gwcIndex), y3(gwcIndex), gwcMark, 'MarkerSize', gwcMarkSize, 'MarkerFaceColor', gwcFill, 'LineWidth', 1.5);
        if END_PERIOD_ALLOPATRY > 0
            endAlloIndex = find(phiData(:,1) <= END_PERIOD_ALLOPATRY, 1, 'last');
            plot(y2(endAlloIndex), y3(endAlloIndex), endAlloSymbol, 'Color', endAlloColor, 'MarkerSize', endAlloSize, 'LineWidth', 1.5);
        end
        if i == 1
            if RI_REACHED && END_PERIOD_ALLOPATRY > 0
                legend('sim. data','uncoupled','coupled','sim. start','sim. end', 'GWC time', ['2' deg ' cntct']);
            elseif RI_REACHED
                legend('sim. data','uncoupled','coupled','sim. start','sim. end', 'GWC time');
            elseif END_PERIOD_ALLOPATRY > 0
                legend('sim. data','uncoupled','coupled','sim. start','sim. end', ['2' deg ' cntct']);
            else
                legend('sim. data','uncoupled','coupled','sim. start','sim. end');
            end
        end
        
        plot(y2(gwcIndex), y3(gwcIndex), gwcMark, 'MarkerSize', gwcMarkSize, 'MarkerFaceColor', gwcFill, 'LineWidth', 1.5);
        if END_PERIOD_ALLOPATRY > 0
            plot(y2(endAlloIndex), y3(endAlloIndex), endAlloSymbol, 'Color', endAlloColor, 'MarkerSize', endAlloSize, 'LineWidth', 1.5);
        end
        
        [tx,ty] = getTextxy(limity2, [y3;10], 'right');
        h = text(tx,ty, ['\bf{' letters{i+plotCols} '}']);
        set(h, 'FontSize', myfs, 'FontName', myfontname);
        
        clear('parameters');
        
        cd ..;
        
        
    end
    
end    


    
    
cd(orig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phiHeaders, phiData, clineWidthHeaders, clineWidthData, barton83headers, barton83ccData, afDiffHeaders, afDiffData, Lheaders, Ldata, meanSheaders, meanSdata] = calculatePhiValuesAndClineWidths(FST, afts, clineMvals, singleLocusQuants, clineSampleSummary)
if exist('parameters.m.bz2', 'file')
    system('bunzip2 parameters.m.bz2');
end
run('./parameters.m');
useClineRow = find(clineMvals == SD_MOVE);
if numel(useClineRow) ~= 1
    error('clineMvals not working how you expect!');
end
times = unique(FST(:,1));
nt = numel(times);
nslq = numel(singleLocusQuants);
phiData = zeros(nt, (nslq+1));
phiData(:,1) = times;
phiHeaders = cell(1, (nslq+1));
phiHeaders{1} = 'time';
for i = 1:nslq
    phiHeaders{(i+1)} = num2str(singleLocusQuants(i));
end
clineWidthHeaders = phiHeaders;
clineWidthData = phiData;

barton83headers = phiHeaders;
barton83ccData = phiData;

afDiffHeaders = phiHeaders;
afDiffData = phiData;

Lheaders = phiHeaders;
Ldata = phiData;

meanSheaders = phiHeaders;
meanSdata = phiData;

isSelected = FST(:,9);
isSelected = (isSelected == 1);
FST = FST(isSelected,:); % limit to data on selected sites
allFSTtimes = FST(:,1); % time sampling points
allFSTsvals = FST(:,5); % selection coefficients
afts = afts(isSelected,:);
aftsTime = afts(:,1);
if ~all(allFSTtimes == aftsTime)
    error('FST and AF data are not congruent!');
end
aftsDeme0 = afts(:,3);
aftsDeme1 = afts(:,4);
aftsDiff = aftsDeme1 - aftsDeme0;
aftsRev = afts(:,5);
for i = 1:numel(aftsDiff)
    if aftsRev(i) == 1
        aftsDiff(i) = -aftsDiff(i);
    end
end
clineWidthTS = 1 ./ aftsDiff;

for i = 1:nslq % column
    lowerBound = clineSampleSummary(useClineRow, i);
    workingRows = aftsDiff >= lowerBound;
    workingS = allFSTsvals(workingRows);
    workingTimes = allFSTtimes(workingRows);
    workingClineWidth = clineWidthTS(workingRows);
    workingAFdiff = aftsDiff(workingRows);
    for j = 1:nt  % row
        rowsAtTime = (workingTimes == times(j));
        L = sum(rowsAtTime); % number of variable loci at time j
        Ldata(j, (i + 1)) = L;
        if L < 1
            phiData(j, (i+1)) = NaN;
            clineWidthData(j, (i + 1)) = NaN;
            barton83ccData(j, (i + 1)) = NaN;
            afDiffData(j, (i + 1)) = NaN;
            meanSdata(j, (i + 1)) = NaN;
        elseif L == 1
            phiData(j, (i+1)) = 0;
            clineWidthData(j, (i + 1)) = workingClineWidth(rowsAtTime);
            barton83ccData(j, (i + 1)) = 0;
            afDiffData(j, (i + 1)) = workingAFdiff(rowsAtTime);
            meanSdata(j, (i + 1)) = workingS(rowsAtTime);
        else
            if MUTATION_DISTRIBUTION == 3 % constant S values
                meanSval = DEME0_CONSTANT_S;
            else
                meanSval = geomean(workingS(rowsAtTime));
            end
            if GAMETE_PRODUCTION_MODE > 0
                r = 0.5;
            else
                r = calcRecomb(L, 0.02 * TOTAL_MAP_LENGTH, nCHROMOSOMES);
            end
            phiData(j, (i+1)) = meanSval * (L - 1) / r;
            clineWidthData(j, (i + 1)) = mean(workingClineWidth(rowsAtTime));
            afDiffData(j, (i + 1)) = mean(workingAFdiff(rowsAtTime));
            barton83ccData(j, (i + 1)) = meanSval / r;
            meanSdata(j, (i + 1)) = meanSval;
        end
        
    end
    
end
clear('parameters');








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writePhiAndClineWidthTStoFile(phiHeaders, phiData, fname)
fpt = fopen(fname, 'w');
ncols = size(phiHeaders,2);
for i = 1:ncols
    if i == ncols
        fprintf(fpt, '%s\n', phiHeaders{i});
    else
        fprintf(fpt, '%s,', phiHeaders{i});
    end
end
nrows = size(phiData,1);
if size(phiData,2) ~= ncols
    error('Phi data and phi headers do not match!');
end
for i = 1:nrows
    for j = 1:ncols
        if j == ncols
            fprintf(fpt, '%E\n', phiData(i,j));
        else
            fprintf(fpt, '%E,', phiData(i,j));
        end
    end
end
fclose(fpt);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1,y1,x2,y2,x3,y3,x4,y4] = restrictDataRange(x1,y1,x2,y2,x3,y3,x4,y4,lowerLimit)
lastIndex = find((y1 >= lowerLimit), 1, 'last');
lastTime = x1(lastIndex);
x1 = x1(1:lastIndex);
y1 = y1(1:lastIndex);
lastx2 = find((x2 <= lastTime), 1, 'last');
x2 = x2(1:lastx2);
y2 = y2(1:lastx2);
x3 = x3(1:lastx2);
y3 = y3(1:lastx2);
x4 = x4(1:lastx2);
y4 = y4(1:lastx2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y] = getTextxy(xdata, ydata, pos)
if nargin < 3
    pos = 'left';
end
xdata = xdata(:);
xdata = xdata(isfinite(xdata));
xdata = xdata( ~isnan(xdata) );
minx = min(xdata);
maxx = max(xdata);
xl = get(gca, 'XLim');
if isfinite(xl(2))
    maxx = max(maxx, xl(2));
end
mxs = get(gca, 'XScale');
if strcmp(pos, 'right')
    xoffset = 0.92;
else
    xoffset = 0.03;
end
if strcmp(mxs, 'log')
    omrange = log10(maxx) - log10(minx);
    x = 10^( (log10(minx) + xoffset * omrange) );
else
    x = minx + xoffset * (maxx - minx);
end
ydata = ydata(:);
ydata = ydata(isfinite(ydata));
ydata = ydata( ~isnan(ydata) );
mys = get(gca, 'YScale');
if strcmp(mys, 'log')
    ydata = ydata( ydata > 0 );
end
miny = min(ydata);
maxy = max(ydata);
yl = get(gca, 'YLim');
if isfinite(yl(2))
    maxy = min(maxy, yl(2));
end
yoffset = 0.95;
if strcmp(mys, 'log')
    omrange = log10(maxy) - log10(miny);
    y = 10^( (log10(miny) + yoffset * omrange) );
else
    y = miny + yoffset * (maxy - miny);
end
if isnan(y)
    disp('point 2:')
    disp([miny maxy yl])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eq = oneLocusEq(m, s)
% solution from Mathematica notebook:
% "/Users/flaxmans/Documents/Research/ProjectsWithPatrik/Regime_shifts/MigrationSelectionBalance.nb"
eq = (m * (0.5 + 0.75 * s) - 0.125 * s - 0.125 * sqrt(s^2 - 4 * m * s^2 + 4 * m^2 * (2 + s)^2))/((-0.25 + m) * s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ccsort, admaxsort] = sortAFDs(cc, admax)
[ccsort, i] = sort(cc);
admaxsort = admax(i);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxWidth = getOneLocusEqCline2(meanSdata)
maxWidth = zeros(size(meanSdata));
run('./parameters.m');
for i = 1:numel(meanSdata)
    peq = oneLocusEq(SD_MOVE, meanSdata(i));
    maxWidth(i) = abs(1 / (peq - (1 - peq)));
end
clear('parameters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function minWidth = getFullyCoupledCline(meanSdata, Ldata)
minWidth = zeros(size(meanSdata));
run('./parameters.m');
for i = 1:numel(meanSdata)
    maxFit = (1 + meanSdata(i)) ^ (Ldata(i));
    maxS = maxFit - 1;
    peq = oneLocusEq(SD_MOVE, maxS);
    minWidth(i) = abs(1 / (peq - (1 - peq)));
end
clear('parameters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [effNumLoci, effectiveS] = calcEffNumLoci( afDiffData, meanSdata )
effNumLoci = zeros(size(afDiffData));
effectiveS = effNumLoci;
n = numel(afDiffData);
run('./parameters.m');
for i = 1:n
    afd = afDiffData(i); % mean afd at a given time
    if afd <= 0
        effNumLoci(i) = 0;
    else
        avgp = (afd + 1)/2;
        
        effectiveS(i) = invertOneLocusEq(avgp, SD_MOVE);
        
        % note that (1 + meanSdata(i))^effNumLoci(i) == 1 + effectiveS
        
        effNumLoci(i) = log(1 + effectiveS(i)) / log(1 + meanSdata(i));
    end
end
clear('parameters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = invertOneLocusEq(p, m)
% see Mathematica notebook on migration selection balance
s = (m * (p - 0.5)) / (((0.25 - (0.25 * p)) * p) + (m * (0.5 + (p * (p - 1.5)))));
