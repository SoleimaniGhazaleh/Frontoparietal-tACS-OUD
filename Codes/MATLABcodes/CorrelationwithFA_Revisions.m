%% ============================================================
%  FULL PIPELINE (ACTIVE GROUP): Craving ↔ BOLD correlations
%  (UPDATED: fixes table() error by storing column lists as STRING)
%
%  What this script does:
%   1) For ACTIVE group only (Group B):
%      - Computes ROI-wise correlations (246 Brainnetome columns) between:
%          X = multiple craving metrics (levels + change scores)
%          Y = a chosen BOLD metric (e.g., post-pre)
%      - Saves per-metric CSV (ROI index = column number), includes:
%          ROI, ROI_Name (optional), r, p_raw, q_FDR, N
%      - Saves MASTER summary CSV across metrics
%      - Saves UNCORRECTED figures (raw p) for each metric
%      - Saves EXACT ROI lists (raw p < 0.05) for each metric
%      - Saves a "COLUMN LISTS" summary with column indices as a string "[..]"
%
%  REQUIRED INPUTS:
%   - Craving table:
%       /Volumes/ExtremeSSDD/LIBR_tACS/VAS/Craving_Final.xlsx
%     columns: ID, Group, Baseline, Before_Pre, After_Pre, Before_Post,
%              After_Post, Day_After
%
%   - BOLD matrices in workspace: base, pre, post  (each 60 x 246)
%     (load them below if needed)
%
%  OPTIONAL:
%   - Brainnetome ROI labels (used only for naming, not required):
%       /Volumes/ExtremeSSDD/LIBR_tACS/OrganizedMATLAB/BNA_ROI_parcellation.mat
%
%  NOTE:
%   - Subject row order MUST match between craving_data and base/pre/post.
%     If not, align by ID before running.
%% ============================================================

clear; clc;
load('/Volumes/ExtremeSSDD/LIBR_tACS/OrganizedMATLAB/BNA_ROI_parcellation.mat');

%% ------------------ User settings ------------------
method    = "Pearson";   % "Pearson" or "Spearman"
minN      = 10;          % minimum N per ROI correlation
alpha_raw = 0.05;        % nominal threshold (uncorrected)
Ntop      = 15;          % top-N ROIs shown in bar plot by smallest raw p

outRoot = '/Volumes/ExtremeSSDD/LIBR_tACS/SecondaryAnalyses/Craving_BOLD_Correlations';

%% ------------------ Load craving table ------------------
filename = '/Volumes/ExtremeSSDD/LIBR_tACS/VAS/Craving_Final.xlsx';
craving_data = readtable(filename);

%% ------------------ Load ROI labels (best effort; optional) ------------------
parcFile = '/Volumes/ExtremeSSDD/LIBR_tACS/OrganizedMATLAB/BNA_ROI_parcellation.mat';
roiNames = load_bna_roi_names(parcFile, 246);  % string(246,1), defaults to "ROI_#"

%% ------------------ Load BOLD activation matrices ------------------
% If needed, load your base/pre/post here:
% load('/path/to/BOLD_Activation_60x246.mat','base','pre','post');

assert(exist('base','var')==1 && exist('pre','var')==1 && exist('post','var')==1, ...
    'Missing base/pre/post variables in workspace. Load them first.');
assert(size(base,2)==246 && size(pre,2)==246 && size(post,2)==246, ...
    'BOLD matrices must be [Nsub x 246].');
assert(height(craving_data)==size(base,1), ...
    'craving_data rows must match BOLD subject count (rows). If not, align by ID.');

%% ------------------ Choose BOLD variable (Y) ------------------
% Pick ONE definition and label:
% Y      = base;        Y_name = "BOLD_Base";
% Y      = pre;         Y_name = "BOLD_Pre";
 Y      = post;        Y_name = "BOLD_Post";
%Y      = post - pre;    Y_name = "BOLD_PostMinusPre";  % recommended
% Y      = post - base; Y_name = "BOLD_PostMinusBase";

%% ------------------ Define ACTIVE group mask ------------------
if iscell(craving_data.Group) || isstring(craving_data.Group)
    grp = string(craving_data.Group);
    isActive = (grp == "B") | contains(grp,"Active","IgnoreCase",true);
else
    % If numeric coding differs, edit this line:
    isActive = (craving_data.Group == 2) | (craving_data.Group == 1);
end
fprintf('Active group N (raw): %d\n', sum(isActive));

%% ------------------ Build craving metrics (levels + changes) ------------------
C = struct();

% Levels
C.Baseline     = craving_data.Baseline;
C.Before_Pre   = craving_data.Before_Pre;
C.After_Pre    = craving_data.After_Pre;
C.Before_Post  = craving_data.Before_Post;
C.After_Post   = craving_data.After_Post;
C.Day_After    = craving_data.Day_After;

% Cue-induced changes
C.DeltaCue_Pre   = craving_data.After_Pre  - craving_data.Before_Pre;
C.DeltaCue_Post  = craving_data.After_Post - craving_data.Before_Post;

% Longitudinal changes
C.Delta_BeforePostMinusBeforePre = craving_data.Before_Post - craving_data.Before_Pre;
C.Delta_AfterPostMinusAfterPre   = craving_data.After_Post  - craving_data.After_Pre;

C.Delta_DayAfterMinusBaseline    = craving_data.Day_After   - craving_data.Baseline;
C.Delta_AfterPostMinusBaseline   = craving_data.After_Post  - craving_data.Baseline;
C.Delta_BeforePostMinusBaseline  = craving_data.Before_Post - craving_data.Baseline;

% Difference-in-differences
C.DeltaDelta_CuePostMinusCuePre  = C.DeltaCue_Post - C.DeltaCue_Pre;

metricNames = fieldnames(C);

%% ------------------ Prepare output folders ------------------
if ~exist(outRoot,'dir'); mkdir(outRoot); end

csvDir = fullfile(outRoot, 'PerMetricCSVs');
figDir = fullfile(outRoot, 'Figures_UNCORRECTED');
repDir = fullfile(outRoot, 'Reports_RawP_ExactROIs');

if ~exist(csvDir,'dir'); mkdir(csvDir); end
if ~exist(figDir,'dir'); mkdir(figDir); end
if ~exist(repDir,'dir'); mkdir(repDir); end

%% ------------------ Initialize outputs ------------------
ROI = (1:246).';
Master = table();
AllNominal = table();        % consolidated list (raw p<0.05) across metrics
ColListSummary = table();    % one row per metric: string "[col1 col2 ...]"

%% ------------------ Run correlations for all metrics ------------------
for m = 1:numel(metricNames)

    X_name = string(metricNames{m});  % ensure STRING
    X_all  = C.(metricNames{m});

    % Subset active group
    X = X_all(isActive);
    Y_act = Y(isActive, :);

    % Remove subjects with missing craving for THIS metric
    valid = isfinite(X);
    X = X(valid);
    Y_act = Y_act(valid, :);

    NactiveUsed = numel(X);
    if NactiveUsed < minN
        fprintf('[SKIP] %s: N=%d < %d\n', X_name, NactiveUsed, minN);
        continue;
    end

    % ROI-wise correlations
    r = nan(246,1);
    p = nan(246,1);
    nUsed = nan(246,1);

    for roi = 1:246
        y = Y_act(:, roi);
        ok = isfinite(X) & isfinite(y);
        nUsed(roi) = sum(ok);

        if nUsed(roi) >= minN
            [R,P] = corr(X(ok), y(ok), 'Type', method, 'Rows', 'complete');
            r(roi) = R;
            p(roi) = P; % raw p
        end
    end

    % FDR across 246 ROIs
    p_for_fdr = p; p_for_fdr(~isfinite(p_for_fdr)) = 1;
    if exist('mafdr','file') == 2
        q = mafdr(p_for_fdr, 'BHFDR', true);
    else
        q = bh_fdr(p_for_fdr);
    end

    % Save per-metric CSV (ROI index = column number)
    ROI_Name = roiNames(ROI);
    T = table(ROI, ROI_Name, r, p, q, nUsed, ...
        'VariableNames', {'ROI','ROI_Name','r','p_raw','q_FDR','N'});
    outCSV = fullfile(csvDir, sprintf('Active_%s_vs_%s_%s.csv', X_name, Y_name, method));
    writetable(sortrows(T,'p_raw','ascend'), outCSV);

    % Nominal hits (raw p < alpha_raw) and COLUMN LIST
    hit = isfinite(p) & (p < alpha_raw);
    colList = ROI(hit);                       % column indices (1..246)
    colListStr = string(mat2str(colList'));   % <-- FIX: store as STRING, not char

    % Add row to Column-list summary (safe table concatenation)
    newRow = table(X_name, string(Y_name), string(method), NactiveUsed, sum(hit), colListStr, ...
        'VariableNames', {'CravingMetric','BOLDMetric','CorrType','N_Active','NumCols_rawP_lt_0p05','ColumnList_ROI_Indices'});
    ColListSummary = [ColListSummary; newRow];

    % Save per-metric exact ROI list (raw p<0.05)
    if any(hit)
        H = table();
        H.CravingMetric = repmat(X_name, sum(hit), 1);
        H.BOLDMetric    = repmat(string(Y_name), sum(hit), 1);
        H.CorrType      = repmat(string(method), sum(hit), 1);

        H.ROI      = ROI(hit);               % <-- column numbers
        H.ROI_Name = roiNames(ROI(hit));     % optional naming
        H.r        = r(hit);
        H.p_raw    = p(hit);
        H.q_FDR    = q(hit);
        H.N        = nUsed(hit);

        H = sortrows(H, 'p_raw', 'ascend');

        outExact = fullfile(repDir, sprintf('ExactROIs_rawP_lt_%.2f_Active_%s_vs_%s_%s.csv', ...
            alpha_raw, X_name, Y_name, method));
        writetable(H, outExact);

        AllNominal = [AllNominal; H]; %#ok<AGROW>
    end

    % MASTER summary row
    nSigFDR  = sum(isfinite(q) & q < 0.05);
    nNominal = sum(hit);

    [~, idxMin] = min(p_for_fdr);
    topROI = idxMin;
    topR   = r(idxMin);
    topP   = p(idxMin);
    topQ   = q(idxMin);

    Master = [Master; table(X_name, string(Y_name), string(method), NactiveUsed, nNominal, nSigFDR, topROI, topR, topP, topQ, ...
        'VariableNames', {'CravingMetric','BOLDMetric','CorrType','N_Active',...
                          'NumNominal_p005','NumSigFDR_q005','TopROI','Top_r','Top_p_raw','Top_q_FDR'})];

    %% ------------------ UNCORRECTED VISUALS (raw p) ------------------
    p_plot = p; p_plot(~isfinite(p_plot)) = 1;
    r_plot = r; r_plot(~isfinite(r_plot)) = NaN;
    sigNominal = (p_plot < alpha_raw) & isfinite(r_plot);

    % 1) Manhattan plot
    f1 = figure('Color','w','Visible','off');
    plot(ROI, -log10(p_plot), 'k.-'); grid on;
    yline(-log10(alpha_raw), '--', sprintf('raw p=%.2f', alpha_raw));
    xlabel('Brainnetome ROI index (column number 1..246)');
    ylabel('-log_{10}(p) [raw]');
    title(sprintf('Active group: %s vs %s (%s)\nUNCORRECTED p-values', X_name, Y_name, method));
    saveas(f1, fullfile(figDir, sprintf('Manhattan_rawP_Active_%s_vs_%s.png', X_name, Y_name)));
    close(f1);

    % 2) r plot + highlights
    f2 = figure('Color','w','Visible','off');
    plot(ROI, r_plot, 'k.-'); hold on; grid on;
    plot(ROI(sigNominal), r_plot(sigNominal), 'ro', 'MarkerFaceColor','r');
    yline(0,'-');
    xlabel('Brainnetome ROI index (column number 1..246)');
    ylabel('Correlation (r)');
    title(sprintf('Active group: %s vs %s (%s)\nRed = raw p<%.2f (UNCORRECTED)', X_name, Y_name, method, alpha_raw));
    legend({'r (all ROIs)','raw p<0.05'}, 'Location','best');
    saveas(f2, fullfile(figDir, sprintf('R_withNominalHighlights_Active_%s_vs_%s.png', X_name, Y_name)));
    close(f2);

    % 3) Top-N by smallest raw p
    f3 = figure('Color','w','Visible','off');
    [~, ord] = sort(p_plot, 'ascend');
    ord = ord(isfinite(r_plot(ord)));
    ordTop = ord(1:min(Ntop, numel(ord)));

    bar(categorical(string(ordTop)), r_plot(ordTop));
    grid on;
    xlabel(sprintf('Top %d ROIs by smallest raw p (ROI indices = columns)', numel(ordTop)));
    ylabel('Correlation (r)');
    title(sprintf('Active group: %s vs %s (%s)\nTOP ROIs (ranked by raw p; UNCORRECTED)', X_name, Y_name, method));
    saveas(f3, fullfile(figDir, sprintf('TopROIs_rawP_Active_%s_vs_%s.png', X_name, Y_name)));
    close(f3);

    fprintf('[DONE] %s | Nominal cols(raw p<0.05)=%d | Cols=%s\n', X_name, nNominal, colListStr);
end

%% ------------------ Save outputs ------------------
masterCSV = fullfile(outRoot, sprintf('MASTER_Active_AllCravingMetrics_vs_%s_%s.csv', Y_name, method));
writetable(Master, masterCSV);
fprintf('\nSaved MASTER summary:\n%s\n', masterCSV);

colListCSV = fullfile(outRoot, sprintf('COLUMN_LISTS_rawP_lt_%.2f_Active_vs_%s_%s.csv', alpha_raw, Y_name, method));
writetable(ColListSummary, colListCSV);
fprintf('Saved COLUMN LISTS summary:\n%s\n', colListCSV);

if ~isempty(AllNominal)
    outAllNom = fullfile(outRoot, sprintf('ALL_Metrics_ExactROIs_rawP_lt_%.2f_Active_vs_%s_%s.csv', ...
        alpha_raw, Y_name, method));
    writetable(AllNominal, outAllNom);
    fprintf('Saved consolidated nominal ROI list:\n%s\n', outAllNom);

    fprintf('\n=== SUMMARY: nominal ROIs (raw p < %.2f) ===\n', alpha_raw);
    [G, metricList] = findgroups(AllNominal.CravingMetric);
    nPer = splitapply(@numel, AllNominal.ROI, G);
    disp(table(metricList, nPer, 'VariableNames', {'CravingMetric','NumROIs_rawP_lt_0p05'}));
else
    fprintf('No ROIs with raw p < %.2f found across metrics.\n', alpha_raw);
end

%% ============================================================
% Helper: BH-FDR fallback (if mafdr is unavailable)
%% ============================================================
function q = bh_fdr(pvals)
    p = pvals(:);
    [ps, idx] = sort(p, 'ascend');
    m = numel(ps);
    q_sorted = ps .* m ./ (1:m)';
    q_sorted = min(1, flipud(cummin(flipud(q_sorted))));
    q = nan(size(p));
    q(idx) = q_sorted;
    q = reshape(q, size(pvals));
end

%% ============================================================
% Helper: robust ROI label loader (best effort)
% - Returns string(nROI,1): uses ROI_### if labels not found
%% ============================================================
function roiNames = load_bna_roi_names(parcFile, nROI)
    roiNames = "ROI_" + string((1:nROI)');

    if ~exist(parcFile,'file')
        warning('Parcellation file not found: %s. Using ROI indices only.', parcFile);
        return;
    end

    S = load(parcFile);
    vars = fieldnames(S);

    % Unwrap common "single struct variable" pattern
    if numel(vars)==1 && isstruct(S.(vars{1}))
        P = S.(vars{1});
    else
        P = S;
    end

    candidateFields = ["ROI_names","roi_names","labels","label","RegionName","RegionNames", ...
                       "name","names","BNA_names","BNA_ROI_names","BNA_label","BNA_labels"];

    for f = candidateFields
        if isfield(P, f)
            tmp = coerce_to_strings(P.(f), nROI);
            if numel(tmp)==nROI
                roiNames = strip(tmp);
                return;
            end
        end
    end

    tmp = try_find_label_vector(P, nROI);
    if ~isempty(tmp)
        roiNames = strip(tmp);
    else
        warning('Could not find ROI names in %s. Using ROI indices only.', parcFile);
    end
end

function names = try_find_label_vector(P, nROI)
    names = [];
    fn = fieldnames(P);
    for k = 1:numel(fn)
        v = P.(fn{k});
        if isstruct(v)
            tmp = try_find_label_vector(v, nROI);
            if ~isempty(tmp), names = tmp; return; end
        else
            if (iscell(v) || isstring(v) || ischar(v)) && numel(v)==nROI
                tmp = coerce_to_strings(v, nROI);
                if numel(tmp)==nROI, names = tmp; return; end
            end
        end
    end
end

function s = coerce_to_strings(x, nROI)
    if ischar(x)
        if size(x,1) == nROI
            s = string(cellstr(x));
        else
            s = strings(nROI,1);
        end
    elseif iscell(x)
        s = string(x(:));
    elseif isstring(x)
        s = x(:);
    else
        s = strings(nROI,1);
    end
end
