function make_OSM_fig_HistogramsOfAFdiffs()

nbins = 100;
mvals = [0.01, 0.02, 0.05, 0.1];
fs = 14;
n = numel(mvals);

close all;
figure;

currmin = 0;
currmax = 0;

for i = 1:n
    
    wm = mvals(i);
    
    data = importdata(['ClineWidthSamples_m' num2str(wm) '.txt']);
    subplot(2,2,i);
    
    currmin = min(currmin, min(data(:,2)));
    currmax = max(currmax, max(data(:,2)));
    
    
    [n,x] = hist(data(:,2), nbins);
    h = bar(x,n);
    set(h,'FaceColor',[1 0.5 0],'EdgeColor','none')
    
    
    set(gca, 'FontSize', fs, 'FontName', 'Arial');
    set(h, 'FaceColor', [1 0.5 0]);
    h = xlabel('Allele frequency difference');
    set(h, 'FontSize', fs, 'FontName', 'Arial');
    h = ylabel('Frequency');
    set(h, 'FontSize', fs, 'FontName', 'Arial');
    h = title(['{\itm} = ' num2str(wm) ', {\its} = 0.02']);
    set(h, 'FontSize', fs, 'FontName', 'Arial');
    
    xlim([-0.09, 0.47]);
    ylim([0,1.01*max(n)]);
    if i == 1
        tx = -0.07;
    else
        tx = 0.2;
    end
    yl = get(gca,'YLim');
    h = text(tx, 0.9 * yl(2), ['mean = ' num2str(round(mean(data(:,2)),3))]);
    set(h, 'FontSize', fs-2, 'FontName', 'Arial');
    h = text(tx, 0.8 * yl(2), ['s.d. = ' num2str(round(std(data(:,2)), 3))]);
    set(h, 'FontSize', fs-2, 'FontName', 'Arial');
    h = text(tx, 0.7 * yl(2), ['# samples = ' num2str(size(data,1))]);
    set(h, 'FontSize', fs-2, 'FontName', 'Arial');
    
end

disp(currmin);
disp(currmax);