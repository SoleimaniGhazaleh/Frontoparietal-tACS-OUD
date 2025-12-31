%% ============================================================
%  Exploratory EF–Craving Correlation (Group B only)
%  - Uncorrected p < 0.05 screening across 210 ROIs
%  - Reports mean±SD EF in Group B for significant ROIs
%  - Generates color-coded scatterplots (Group B color) with:
%       * fixed window size
%       * y-axis labels that INCLUDE timing (T1..T4 definitions)
%  - Saves results table + PNG figures
%
%  Inputs:
%    EF:  /Volumes/ExtremeSSDD/LIBR_tACS/Final_AB_data/EF.mat   (EF_B, EF_A)
%    VAS: /Volumes/ExtremeSSDD/LIBR_tACS/VAS/Craving_Final_comparison.xlsx
%         columns: ID, Group, T0..T5
%
%  Contrasts:
%    deltaCue_Post  = T4 - T2
%    deltaTask_Pre  = T2 - T1
%    deltaTask_Post = T4 - T3
%    deltaDeltaTask = (T4 - T3) - (T2 - T1)
%% ============================================================

clear; clc; close all;

%% ------------------ Paths ------------------
efFile   = '/Volumes/ExtremeSSDD/LIBR_tACS/Final_AB_data/EF.mat';
vasFile  = '/Volumes/ExtremeSSDD/LIBR_tACS/VAS/Craving_Final_comparison.xlsx';

outDir   = '/Volumes/ExtremeSSDD/LIBR_tACS/VAS/EF_Craving_GroupB_uncorrected/';
figDir   = fullfile(outDir, 'ScatterPlots');
if ~exist(outDir, 'dir'); mkdir(outDir); end
if ~exist(figDir, 'dir'); mkdir(figDir); end

outCSV   = fullfile(outDir, 'EF_Craving_GroupB_uncorrected_p05_withEFmean.csv');

%% ------------------ Colors ------------------
colA = [248 118 109]/255; % #F8766D  Group A = Sham (not used in Group-B-only scatter)
colB = [0 191 196]/255;   % #00BFC4  Group B = Active

%% ------------------ Figure window size ------------------
% Position = [left bottom width height] in pixels
figPos = [100 100 300 300];

%% ------------------ Load EF ------------------
S = load(efFile);  % expects EF_B and EF_A
if isfield(S,'EF_B')
    EF = S.EF_B;
elseif isfield(S,'EF')
    EF = S.EF;
else
    error('EF file does not contain EF_B or EF.');
end

if ~ismatrix(EF)
    error('EF must be a 2D matrix [Nsubj x 210].');
end

nROI = size(EF,2);
if nROI ~= 210
    warning('EF has %d columns (expected 210). Continuing...', nROI);
end

%% ------------------ Load VAS ------------------
T = readtable(vasFile);

reqVars = {'ID','Group','T0','T1','T2','T3','T4','T5'};
missingVars = setdiff(reqVars, T.Properties.VariableNames);
if ~isempty(missingVars)
    error('Missing columns in VAS Excel: %s', strjoin(missingVars, ', '));
end

T.ID    = string(T.ID);
T.Group = categorical(string(T.Group));

cats = categories(T.Group);
if numel(cats) ~= 2
    error('Expected exactly 2 groups in VAS file. Found %d: %s', numel(cats), strjoin(cats, ', '));
end

% Assumption: 1st category = Group A (Sham), 2nd category = Group B (Active)
grpB = cats{2};

% Keep Group B only
T = T(T.Group == grpB, :);

%% ------------------ Build contrasts ------------------
T.deltaCue_Post   = T.T4 - T.T2;
T.deltaTask_Pre   = T.T2 - T.T1;
T.deltaTask_Post  = T.T4 - T.T3;
T.deltaDeltaTask  = T.deltaTask_Post - T.deltaTask_Pre;

%% ------------------ Sanity check: rows match ------------------
if size(EF,1) ~= height(T)
    error(['Mismatch between EF rows (%d) and Group B VAS rows (%d).\n' ...
           'Make sure EF_B rows are in the SAME subject order as VAS Group B rows.\n' ...
           'If not, align by ID before correlating.'], size(EF,1), height(T));
end

%% ------------------ Variables to test + timing-aware labels ------------------
cravingVars = { ...
    'deltaCue_Post', ...
    'deltaTask_Pre', ...
    'deltaTask_Post', ...
    'deltaDeltaTask'};

% Human-readable y-axis labels including timing
cravingLabels = containers.Map( ...
    {'deltaCue_Post', ...
     'deltaTask_Pre', ...
     'deltaTask_Post', ...
     'deltaDeltaTask'}, ...
    { ...
     'Cue-induced craving change (Post-fMRI: T4 − T2)', ...
     'Task-induced craving change (Pre-stimulation: T2 − T1)', ...
     'Task-induced craving change (Post-stimulation: T4 − T3)', ...
     'Task-induced craving ΔΔ (Post − Pre): (T4 − T3) − (T2 − T1)' ...
    });

%% ------------------ Screen correlations (uncorrected p<0.05) ------------------
results = table();

for v = 1:numel(cravingVars)

    varName = cravingVars{v};
    y = T.(varName);

    for r = 1:nROI

        x = EF(:,r);

        ok = ~isnan(x) & ~isnan(y);
        if sum(ok) < 6
            continue
        end

        [rho, p] = corr(x(ok), y(ok), 'Type','Pearson');

        if p < 0.05
            results = [results;
                table( ...
                    r, ...
                    string(varName), ...
                    rho, ...
                    p, ...
                    sum(ok), ...
                    mean(x(ok)), ...
                    std(x(ok)), ...
                    'VariableNames', ...
                    {'ROI','CravingContrast','r','p_uncorrected','N','EF_mean','EF_std'})];
        end
    end
end

%% ------------------ Sort, display, save results ------------------
if isempty(results)
    disp('No uncorrected p < 0.05 correlations found in Group B.');
    writetable(results, outCSV);
else
    results = sortrows(results, 'p_uncorrected');
    disp('=== Significant EF–Craving correlations (Group B, uncorrected) ===');
    disp(results);
    writetable(results, outCSV);
    fprintf('Saved results table: %s\n', outCSV);
end

%% ------------------ Scatter plots (Group B color-coded + timing in y-label) ------------------
if ~isempty(results)

    for i = 1:height(results)

        roiIdx = results.ROI(i);
        varName = char(results.CravingContrast(i)); % e.g., 'deltaCue_Post'

        x = EF(:, roiIdx);
        y = T.(varName);

        ok = ~isnan(x) & ~isnan(y);
        if sum(ok) < 6
            continue
        end

        % Recompute correlation on plotted points (matches scatter)
        [rho, p] = corr(x(ok), y(ok), 'Type','Pearson');
        N = sum(ok);

        % Figure with fixed window size
        f = figure('Color','w', 'Position', figPos);

        scatter(x(ok), y(ok), 60, 'filled', ...
            'MarkerFaceColor', colB, 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.85);
        hold on;

        % Regression line
        pp = polyfit(x(ok), y(ok), 1);
        xx = linspace(min(x(ok)), max(x(ok)), 100);
        yy = polyval(pp, xx);
        plot(xx, yy, 'k-', 'LineWidth', 2);

        % Axis labels (y includes timing)
        xlabel(sprintf('Electric field (Group B) — ROI %d', roiIdx), 'Interpreter','none');

        if isKey(cravingLabels, varName)
            ylabel(cravingLabels(varName), 'Interpreter','none');
        else
            ylabel(varName, 'Interpreter','none');
        end

        % EF mean±SD on plotted points
        efMean = mean(x(ok));
        efStd  = std(x(ok));

        title(sprintf('ROI %d | r=%.3f, p=%.4f | EF=%.3f±%.3f', ...
            roiIdx, rho, p, efMean, efStd), 'Interpreter','none');

        grid on; box on;

        % Save figure
        safeContrast = regexprep(varName, '[^a-zA-Z0-9_]', '_');
        outPng = fullfile(figDir, sprintf('GroupB_ROI%03d_%s_scatter.png', roiIdx, safeContrast));
        exportgraphics(f, outPng, 'Resolution', 300);
        close(f);
    end

    fprintf('Saved scatterplots to: %s\n', figDir);
end

%% ------------------ Done ------------------
disp('DONE.');
