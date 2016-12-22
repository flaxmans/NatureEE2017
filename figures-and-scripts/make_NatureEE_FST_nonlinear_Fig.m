function make_NatureEE_FST_nonlinear_Fig(indexes)
% script for making panels a, b, and d of Figure 2 in the NatureEE paper.
% These figures use data that was previous published and available at
% DataDryad.org:  
% original publication:
% Flaxman SM, Wacholder AC, Feder JL, Nosil P (2014) Theoretical models of the influence of genomic architecture on the dynamics of speciation. Molecular Ecology 23(16): 4074-4088. http://dx.doi.org/10.1111/mec.12750
% Dryad data package:
% Flaxman SM, Wacholder AC, Feder JL, Nosil P (2014) Data from: Theoretical models of the influence of genomic architecture on the dynamics of speciation. Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.kc596


myfs = 12;
myfontname = 'Arial';
strongDiff = 0.25; % Hartl and Clark's FST definition
myms = 6;
mycolor = [1 0 0];
mylw = 2;

if nargin < 1
    indexes = [25301 26201];
end
filename = ['ws_NatureEE_FST_fig_' num2str(indexes(1)) '_' num2str(indexes(2)) '.mat'];
if exist(filename,'file')
    load(filename);
else
    
    for i = 1:2
        orig = pwd;
        cd(['Flaxman_et_al_MolEcol_Dryad/Figure1_Data_and_Code/Run' num2str(indexes(i))]);
        
        if exist('FSTtimeSeries.txt.bz2', 'file')
            system('bunzip2 FSTtimeSeries.txt.bz2');
        end
        
        if i == 1
            FST1 = importdata('FSTtimeSeries.txt');
            times1 = FST1(:,1);
            FSTvals1 = FST1(:,3);
            lociIDs1 = FST1(:,2);

            utimes1 = unique(times1);
            meanFST1 = zeros(size(utimes1));
            numStrongDiff1 = meanFST1;
            lociCounted1 = [];

            for j = 1:numel(utimes1)
                wrows = times1 == utimes1(j);
                wFSTvals = FSTvals1(wrows);
                meanFST1(j) = mean(wFSTvals);
                wlids = lociIDs1(wrows);
                FSTabove = (wFSTvals >= strongDiff);
                nabove = sum(FSTabove);
                if nabove > 0
                    lociCounted1 = union(lociCounted1, wlids(FSTabove));
                end
                numStrongDiff1(j) = numel(lociCounted1);
                
            end
        else
            FST2 = importdata('FSTtimeSeries.txt');
            times2 = FST2(:,1);
            FSTvals2 = FST2(:,3);
            lociIDs2 = FST2(:,2);

            utimes2 = unique(times2);
            meanFST2 = zeros(size(utimes2));
            numStrongDiff2 = meanFST2;
            lociCounted2 = [];

            for j = 1:numel(utimes2)
                wrows = times2 == utimes2(j);
                wFSTvals = FSTvals2(wrows);
                meanFST2(j) = mean(wFSTvals);
                wlids = lociIDs2(wrows);
                FSTabove = (wFSTvals >= strongDiff);
                nabove = sum(FSTabove);
                if nabove > 0
                    lociCounted2 = union(lociCounted2, wlids(FSTabove));
                end
                numStrongDiff2(j) = numel(lociCounted2);
                
            end
        end
        
        cd(orig);
    end
    
    if ~exist(filename, 'file')
        save(filename);
    end
end

close all;

figure;
threshLocCount = 200;
bot = -5;
tx = 0.03;
ty = 0.9 * threshLocCount;
% first panel
subplot(3,1,1);
l = find(numStrongDiff1 <= threshLocCount, 1, 'last');
plot(utimes1(1:l), numStrongDiff1(1:l), 'Color', mycolor, 'LineWidth', mylw);
hold on;
set(gca, 'FontSize', myfs, 'FontName', myfontname, 'XLim', [0,utimes1(l)], 'YLim', [bot, threshLocCount]);
h = xlabel('Time (generations)');
set(h, 'FontSize', myfs, 'FontName', myfontname);
h = ylabel('# strongly diff. loci');
set(h, 'FontSize', myfs, 'FontName', myfontname);
run(['Flaxman_et_al_MolEcol_Dryad/Figure1_Data_and_Code/Run' num2str(indexes(1)) '/parameters.m']);
h = text(tx * utimes1(l), ty, ['{\itm} = ' num2str(SD_MOVE) ', {\its} = ' num2str(MEAN_S)]);
set(h, 'FontSize', myfs - 1, 'FontName', myfontname);
% second panel
subplot(3,1,2);
l = find(numStrongDiff2 <= threshLocCount, 1, 'last');
plot(utimes2(1:l), numStrongDiff2(1:l), 'Color', mycolor, 'LineWidth', mylw);
hold on;
set(gca, 'FontSize', myfs, 'FontName', myfontname, 'XLim', [0,utimes2(l)], 'YLim', [bot, threshLocCount]);
h = xlabel('Time (generations)');
set(h, 'FontSize', myfs, 'FontName', myfontname);
h = ylabel('# strongly diff. loci');
set(h, 'FontSize', myfs, 'FontName', myfontname);
run(['Flaxman_et_al_MolEcol_Dryad/Figure1_Data_and_Code/Run' num2str(indexes(2)) '/parameters.m']);
h = text(tx * utimes2(l), ty, ['{\itm} = ' num2str(SD_MOVE) ', {\its} = ' num2str(MEAN_S)]);
set(h, 'FontSize', myfs - 1, 'FontName', myfontname);

thirdPanel(myfs, myfontname, myms)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thirdPanel(myfs, myfontname, myms)
grwt = 0.7;
load('Flaxman_et_al_MolEcol_Dryad/Figure4_Data_and_code/BistabilityFigureData.mat');
subplot(3,1,3);
plot(s_hyst,(1-hyst_above),'>', 'Color', [1 0.5 0], 'MarkerSize', myms);
hold on;
plot(s_hyst,(1-hyst_below),'<', 'Color', [0 0.5 1], 'MarkerSize', myms);
% h=legend('{\it s} increasing','{\it s} decreasing');
% set(h,'FontSize',10);
plot(s_hyst,(1-hyst_above),'--','Color',[grwt,grwt,grwt]);
plot(s_hyst,(1-hyst_below),'--','Color',[grwt,grwt,grwt]);
set(gca,'FontSize',myfs,'FontName',myfontname);
h = xlabel('{\it s} (per locus)');
set(h, 'FontSize', myfs, 'FontName', myfontname);
h = ylabel('Allele freq. in favored deme');
set(h, 'FontSize', myfs, 'FontName', myfontname);
ylim([0.5,1.0]);
xlim([min(s_hyst)-0.0002,max(s_hyst)+0.0002]);
h = text(0.035, 0.95, ['{\itm} = 0.1, {\itL} = ' num2str(nLOCI_Hysteresis)]);
set(h, 'FontSize', myfs - 1, 'FontName', myfontname);

