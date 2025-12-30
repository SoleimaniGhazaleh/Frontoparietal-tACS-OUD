%% ============================================================
%  Scatter plots: Craving vs BOLD (Active group only)
%  UPDATED (reviewer-safe, fully automatic):
%   - Automatically selects metricsToPlot from your saved COLUMN_LISTS CSV
%   - Automatically selects roiList PER METRIC from ColumnList_ROI_Indices
%   - Fixed figure size in pixels (600x600)
%   - Active color (#00BFC4)
%   - Dot size + edge color/width configurable
%   - Robust bivariate outlier removal (MAD/median) + optional plotting
%   - Uses the SAME correlation type as in your CSV name (Pearson/Spearman)
%
%  IMPORTANT:
%   - This script only makes scatter plots for ROIs that were nominal
%     (raw p<0.05) in your screening step.
%   - It will stay consistent even if screening results change.
%% ============================================================

clear; clc;
 load('/Volumes/ExtremeSSDD/LIBR_tACS/OrganizedMATLAB/BNA_ROI_parcellation.mat');


%% ------------------ USER SETTINGS ------------------
% Use same method as your screening outputs (match your filename)
method = 'Pearson';  % 'Pearson' or 'Spearman'

% Choose BOLD measure used in screening:
% If your CSV is ...vs_BOLD_Post... then set BOLD = post
% If your CSV is ...vs_BOLD_PostMinusPre... then set BOLD = post-pre, etc.
BOLD_choice = 'Post';   % 'Post' | 'Pre' | 'Base' | 'PostMinusPre' | 'PostMinusBase'

% Outlier removal
removeOutliers = true;        % apply outlier removal before computing correlation
showOutliers   = true;        % plot excluded points in gray
outlierMethod  = 'median';    % isoutlier(...,'median') = MAD/median rule

% Minimum N after NaN + outlier removal
minN = 8;

% Marker / style
colB = hex2rgb('#00BFC4');    % Active face color
dotSize   = 40;               % scatter area in points^2
edgeColor = 'k';              % marker edge
edgeWidth = 0.8;              % edge width

% Figure size (pixels)
figPos = [100 100 250 250];

% Paths
cravingFile = '/Volumes/ExtremeSSDD/LIBR_tACS/VAS/Craving_Final.xlsx';

% This is the key file created by your screening script:
colListCSV  = '/Volumes/ExtremeSSDD/LIBR_tACS/SecondaryAnalyses/Craving_BOLD_Correlations/COLUMN_LISTS_rawP_lt_0.05_Active_vs_BOLD_Post_Pearson.csv';

% Output directory for plots
outDir = '/Volumes/ExtremeSSDD/LIBR_tACS/SecondaryAnalyses/Craving_BOLD_Correlations/ScatterPlots_AutoMetrics_AutoROIs';
if ~exist(outDir,'dir'); mkdir(outDir); end

%% ------------------ Load craving data ------------------
craving_data = readtable(cravingFile);

%% ------------------ Load BOLD data ------------------
% base, pre, post must exist OR load them here:
% load('/path/to/BOLD_Activation_60x246.mat','base','pre','post');

assert(exist('base','var')==1 && exist('pre','var')==1 && exist('post','var')==1, ...
    'Load base, pre, post matrices first.');
assert(size(base,2)==246 && size(pre,2)==246 && size(post,2)==246, ...
    'base/pre/post must be [Nsub x 246].');
assert(height(craving_data)==size(base,1), ...
    'Row mismatch: craving_data height must equal number of rows in base/pre/post.');

% Choose BOLD measure to match screening
switch lower(BOLD_choice)
    case 'post'
        BOLD = post;          BOLD_name = 'BOLD_Post';
    case 'pre'
        BOLD = pre;           BOLD_name = 'BOLD_Pre';
    case 'base'
        BOLD = base;          BOLD_name = 'BOLD_Base';
    case 'postminuspre'
        BOLD = post - pre;    BOLD_name = 'BOLD_PostMinusPre';
    case 'postminusbase'
        BOLD = post - base;   BOLD_name = 'BOLD_PostMinusBase';
    otherwise
        error('Unknown BOLD_choice: %s', BOLD_choice);
end

%% ------------------ Active group mask ------------------
grp = string(craving_data.Group);
isActive = (grp == "B") | contains(grp,"Active","IgnoreCase",true);
fprintf('Active group N (raw): %d\n', sum(isActive));

%% ------------------ Build craving metrics ------------------
C = struct();
C.Baseline     = craving_data.Baseline;
C.Before_Pre   = craving_data.Before_Pre;
C.After_Pre    = craving_data.After_Pre;
C.Before_Post  = craving_data.Before_Post;
C.After_Post   = craving_data.After_Post;
C.Day_After    = craving_data.Day_After;

C.DeltaCue_Pre   = craving_data.After_Pre  - craving_data.Before_Pre;
C.DeltaCue_Post  = craving_data.After_Post - craving_data.Before_Post;

C.Delta_BeforePostMinusBeforePre = craving_data.Before_Post - craving_data.Before_Pre;
C.Delta_AfterPostMinusAfterPre   = craving_data.After_Post  - craving_data.After_Pre;

C.Delta_DayAfterMinusBaseline    = craving_data.Day_After   - craving_data.Baseline;
C.Delta_AfterPostMinusBaseline   = craving_data.After_Post  - craving_data.Baseline;
C.Delta_BeforePostMinusBaseline  = craving_data.Before_Post - craving_data.Baseline;

C.DeltaDelta_CuePostMinusCuePre  = C.DeltaCue_Post - C.DeltaCue_Pre;

%% ------------------ Read screening output and define WHAT TO PLOT ------------------
Tcols = readtable(colListCSV);

% Ensure string columns
Tcols.CravingMetric = string(Tcols.CravingMetric);

% Keep only metrics with at least one ROI
hasAny = Tcols.NumCols_rawP_lt_0p05 > 0;
Tcols = Tcols(hasAny, :);

% OPTIONAL safety check: method/BOLD consistency (won't error, just warn)
if ismember('CorrType', Tcols.Properties.VariableNames)
    if any(string(Tcols.CorrType) ~= string(method))
        warning('COLUMN_LISTS CorrType does not fully match method="%s". Make sure you are consistent.', method);
    end
end
if ismember('BOLDMetric', Tcols.Properties.VariableNames)
    if any(~contains(string(Tcols.BOLDMetric), strrep(BOLD_name,'BOLD_',''), 'IgnoreCase', true))
        warning('COLUMN_LISTS BOLDMetric may not match BOLD_choice="%s". Check consistency.', BOLD_choice);
    end
end

metricsToPlot = cellstr(Tcols.CravingMetric);
fprintf('Metrics to plot (auto): %s\n', strjoin(metricsToPlot, ', '));

%% ------------------ MAIN LOOP (metric-specific ROI lists) ------------------
for i = 1:height(Tcols)

    metric = string(Tcols.CravingMetric(i));

    if ~isfield(C, metric)
        warning('Metric "%s" not found in struct C. Skipping.', metric);
        continue;
    end

    % Parse ROI list string like "[90 112 126]" or "177" or "zeros(1,0)"
    roiList = parse_roi_list(Tcols.ColumnList_ROI_Indices(i));

    if isempty(roiList)
        fprintf('[SKIP] %s: empty ROI list\n', metric);
        continue;
    end

    % Prepare X and Y for active group
    X_all = C.(metric);
    X = X_all(isActive);
    Y = BOLD(isActive, :);

    for rIdx = 1:numel(roiList)
        roi = roiList(rIdx);

        if roi < 1 || roi > size(Y,2)
            warning('ROI %d out of range (1..%d). Skipping.', roi, size(Y,2));
            continue;
        end

        y = Y(:, roi);

        % Remove NaNs
        ok0 = isfinite(X) & isfinite(y);
        X0 = X(ok0);
        y0 = y(ok0);

        if numel(X0) < minN
            fprintf('[SKIP] %s | ROI %d: N=%d < %d (after NaN removal)\n', metric, roi, numel(X0), minN);
            continue;
        end

        % Outlier removal (bivariate)
        outMask = false(size(X0));
        if removeOutliers
            outMask = any(isoutlier([X0 y0], outlierMethod), 2);
        end

        keep = ~outMask;
        Xp = X0(keep);
        yp = y0(keep);

        N_keep = numel(Xp);
        N_out  = sum(outMask);

        if N_keep < minN
            fprintf('[SKIP] %s | ROI %d: N=%d < %d (after outlier removal)\n', metric, roi, N_keep, minN);
            continue;
        end

        % Correlation on cleaned data
        [R,P] = corr(Xp, yp, 'Type', method, 'Rows','complete');

        %% --------- Plot ---------
        f = figure('Color','w','Visible','off');
        set(f, 'Units','pixels', 'Position', figPos);

        hold on;

        % Plot outliers in gray (optional)
        if showOutliers && removeOutliers && N_out > 0
            scatter(X0(outMask), y0(outMask), dotSize, [0.75 0.75 0.75], ...
                'filled', 'MarkerEdgeColor', [0.55 0.55 0.55], 'LineWidth', edgeWidth);
        end

        % Plot kept points (Active color)
        scatter(Xp, yp, dotSize, colB, 'filled', ...
            'MarkerEdgeColor', edgeColor, 'LineWidth', edgeWidth);

        % Fit line on kept points
        if N_keep >= 2
            pfit = polyfit(Xp, yp, 1);
            xx = linspace(min(Xp), max(Xp), 100);
            yy = polyval(pfit, xx);
            plot(xx, yy, 'k-', 'LineWidth', 1.2);
        end

        grid on;

        xlabel(strrep(metric,'_',' '), 'Interpreter','none');
        ylabel(sprintf('%s – ROI %d', BOLD_name, roi), 'Interpreter','none');

        ttl = sprintf('%s | ROI %d', strrep(metric,'_',' '), roi);
        title(ttl, 'Interpreter','none');

        if removeOutliers
            txt = sprintf('%s (outliers removed: %s)\nr = %.2f\np = %.3f\nN kept = %d, out = %d', ...
                method, outlierMethod, R, P, N_keep, N_out);
        else
            txt = sprintf('%s (no outlier removal)\nr = %.2f\np = %.3f\nN = %d', ...
                method, R, P, N_keep);
        end

        text(0.05, 0.95, txt, 'Units','normalized', ...
            'VerticalAlignment','top', 'FontSize', 11);

        if showOutliers && removeOutliers && N_out > 0
            %legend({'Excluded outliers','Included points','Fit (included)'}, 'Location','best');
        else
            %legend({'Included points','Fit (included)'}, 'Location','best');
        end

        % Save figure
        fname = sprintf('Scatter_Active_%s_%s_ROI_%03d_%s.png', ...
            metric, BOLD_name, roi, method);
        saveas(f, fullfile(outDir, fname));
        close(f);

        fprintf('[SAVED] %s | ROI %d | r=%.2f | p=%.3f | N kept=%d | out=%d\n', ...
            metric, roi, R, P, N_keep, N_out);
    end
end

fprintf('\nAll scatter plots saved to:\n%s\n', outDir);

%% ============================================================
% Helper: parse ROI list stored as string in CSV
% Handles: "[90 112 126]" or "177" or "zeros(1,0)"
%% ============================================================
function roiList = parse_roi_list(s)
    s = string(s);
    s = strtrim(s);

    if strlength(s)==0 || contains(lower(s), "zeros")
        roiList = [];
        return;
    end

    % remove brackets if present
    s = erase(s, "[");
    s = erase(s, "]");
    s = strtrim(s);

    if strlength(s)==0
        roiList = [];
        return;
    end

    nums = sscanf(s, '%f');
    roiList = unique(round(nums(:)')) ; % row vector of integers
    roiList = roiList(roiList>=1 & roiList<=246);
end

%% ------------------ Helper: hex to rgb ------------------
function rgb = hex2rgb(hex)
    hex = erase(string(hex),'#');
    rgb = reshape(sscanf(hex,'%2x')/255, 1, 3);
end
