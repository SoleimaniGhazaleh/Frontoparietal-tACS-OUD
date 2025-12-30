%% ============================================================
%  FULL PIPELINE:
%   1) Load EF + FA
%   2) Compute ΔFA (post-pre) for Group A and B
%   3) ROI-wise correlation (Group B only) between EF and ΔFA (p-uncorrected)
%   4) Report significant ROIs
%   5) Scatterplots EF vs ΔFA with color coding (Sham vs Active)
%% ============================================================

clear; clc; close all;

%% ------------------ Paths ------------------
efPath = '/Volumes/ExtremeSSDD/LIBR_tACS/Final_AB_data/EF.mat';
faPath = '/Volumes/ExtremeSSDD/LIBR_tACS/Final_AB_data/FA.mat';

outDir   = '/Volumes/ExtremeSSDD/LIBR_tACS/Final_AB_data/EF_FA_Corr_GroupB/';
figDir   = fullfile(outDir, 'ScatterPlots');
if ~exist(outDir,'dir'), mkdir(outDir); end
if ~exist(figDir,'dir'), mkdir(figDir); end

%% ------------------ Settings ------------------
nROI     = 210;
pThresh  = 0.05;         % uncorrected threshold
corrType = 'Pearson';    % 'Pearson' or 'Spearman'
minN     = 6;            % minimum subjects required per ROI to compute correlation

%% ------------------ Colors (your exact codes) ------------------
colA = hex2rgb('#F8766D'); % Sham (Group A)
colB = hex2rgb('#00BFC4'); % Active (Group B)

%% ------------------ Load ------------------
load(efPath);  % expects EF_A/EF_B or EF (or similar)
load(faPath);  % expects pre_A, post_A, pre_B, post_B

%% ------------------ Compute ΔFA ------------------
FA_A = post_A;
FA_B = post_B;

FA_A = FA_A(:,1:nROI);
FA_B = FA_B(:,1:nROI);

%% ------------------ Find EF_A and EF_B ------------------
EF_A_use = pickEFmatrix({'EF_A','ef_A','EFA','efa','EF_sham','ef_sham','EF_A_use'}, size(FA_A,1), nROI);
EF_B_use = pickEFmatrix({'EF_B','ef_B','EFB','efb','EF_active','ef_active','EF_B_use'}, size(FA_B,1), nROI);

% If EF_A/EF_B are not present but EF exists, try to infer:
if isempty(EF_A_use) || isempty(EF_B_use)
    if exist('EF','var') && isnumeric(EF)
        % Common cases:
        %  (1) EF is a struct with fields EF_A / EF_B (handled above)
        %  (2) EF is concatenated [A; B]
        %  (3) EF is only group-specific and you loaded wrong file
        if size(EF,2) >= nROI
            if isempty(EF_A_use) && size(EF,1) == size(FA_A,1)
                EF_A_use = EF(:,1:nROI);
            end
            if isempty(EF_B_use) && size(EF,1) == size(FA_B,1)
                EF_B_use = EF(:,1:nROI);
            end
            if isempty(EF_A_use) && isempty(EF_B_use)
                if size(EF,1) == size(FA_A,1) + size(FA_B,1)
                    EF_A_use = EF(1:size(FA_A,1), 1:nROI);
                    EF_B_use = EF(size(FA_A,1)+1:end, 1:nROI);
                end
            end
        end
    end
end

if isempty(EF_A_use) || isempty(EF_B_use)
    whos('-file', efPath)
    error(['Could not identify EF_A and EF_B matrices automatically. ' ...
           'Open EF.mat and tell me the EF variable names (or edit pickEFmatrix list).']);
end

%% ------------------ ROI labels (optional) ------------------
roiNames = strings(nROI,1);
if exist('roi_labels','var') && numel(roi_labels) >= nROI
    roiNames = string(roi_labels(1:nROI));
elseif exist('labels','var') && numel(labels) >= nROI
    roiNames = string(labels(1:nROI));
else
    for r = 1:nROI
        roiNames(r) = sprintf('ROI_%03d', r);
    end
end

%% ============================================================
%  STEP 1: Correlation in Group B (Active) only
%% ============================================================

rVals = nan(nROI,1);
pVals = nan(nROI,1);
nUsed = nan(nROI,1);

for roi = 1:nROI
    x = EF_B_use(:,roi);   % EF in Group B
    y = FA_B(:,roi);       % ΔFA in Group B

    ok = isfinite(x) & isfinite(y);
    nUsed(roi) = sum(ok);

    if nUsed(roi) >= minN
        [rVals(roi), pVals(roi)] = corr(x(ok), y(ok), 'Type', corrType, 'Rows', 'complete');
    end
end

T_all = table((1:nROI)', roiNames, rVals, pVals, nUsed, ...
    'VariableNames', {'ROI_Index','ROI_Name','r_GroupB','p_uncorrected','N_used'});

sigIdx = isfinite(pVals) & pVals < pThresh;
T_sig  = sortrows(T_all(sigIdx,:), 'p_uncorrected', 'ascend');

%% ------------------ Print + Save results ------------------
fprintf('\n=== Group B (Active): EF vs ΔFA correlations | %s | p < %.3f (uncorrected) ===\n', corrType, pThresh);
if isempty(T_sig)
    fprintf('No ROIs reached p < %.3f (uncorrected).\n', pThresh);
else
    disp(T_sig);
end

writetable(T_all, fullfile(outDir, 'GroupB_EF_FA_Corr_AllROIs.csv'));
writetable(T_sig, fullfile(outDir, 'GroupB_EF_FA_Corr_SignificantROIs_pUncorrected.csv'));

%% ============================================================
%  STEP 2: Scatterplots for significant ROIs (A vs B colored)
%% ============================================================

if isempty(T_sig)
    fprintf('\nNo scatter plots created because no ROIs were significant.\n');
else
    for i = 1:height(T_sig)

        roi = T_sig.ROI_Index(i);
        roiName = T_sig.ROI_Name(i);

        % Data
        xA = EF_A_use(:,roi);  yA = FA_A(:,roi);
        xB = EF_B_use(:,roi);  yB = FA_B(:,roi);

        okA = isfinite(xA) & isfinite(yA);
        okB = isfinite(xB) & isfinite(yB);

        % Figure
        f = figure('Color','w','Position',[300 300 560 460]); hold on;

        scatter(xA(okA), yA(okA), 55, ...
            'MarkerFaceColor', colA, 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.80);

        scatter(xB(okB), yB(okB), 55, ...
            'MarkerFaceColor', colB, 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.80);

        % Regression lines (separate per group)
        addFitLine(xA(okA), yA(okA), colA);
        addFitLine(xB(okB), yB(okB), colB);

        xlabel('Electric Field Strength', 'FontSize', 12);
        ylabel('\DeltaFA (Post - Pre)', 'FontSize', 12);

        title(sprintf('%s | Group B: r = %.2f, p = %.4f', ...
            roiName, T_sig.r_GroupB(i), T_sig.p_uncorrected(i)), ...
            'Interpreter','none', 'FontSize', 12);

        legend({'Sham (A)','Active (B)'}, 'Location','best', 'Box','off');

        grid on; box on;

        % Save
        safeName = regexprep(char(roiName), '[^\w\-]', '_'); % filename-safe
        pngName = sprintf('ROI_%03d_%s_EF_vs_dFA.png', roi, safeName);
        saveas(f, fullfile(figDir, pngName));
        close(f);
    end

    fprintf('\nSaved scatter plots to:\n%s\n', figDir);
end

%% ============================================================
%  Local helper functions
%% ============================================================

function rgb = hex2rgb(hex)
    hex = strrep(hex,'#','');
    if numel(hex) ~= 6
        error('hex2rgb: hex string must be 6 characters (e.g., #F8766D).');
    end
    rgb = reshape(sscanf(hex,'%2x'),1,[])/255;
end

function M = pickEFmatrix(nameList, nSub, nROI)
    % Searches workspace variables by name and returns a numeric [nSub x nROI] matrix if found
    M = [];
    for i = 1:numel(nameList)
        nm = nameList{i};
        if evalin('base', sprintf('exist(''%s'',''var'')', nm))
            X = evalin('base', nm);
            if isnumeric(X) && size(X,1) == nSub && size(X,2) >= nROI
                M = X(:,1:nROI);
                return;
            elseif isstruct(X)
                % If it’s a struct, attempt to find a numeric field matching [nSub x >=nROI]
                fns = fieldnames(X);
                for j = 1:numel(fns)
                    v = X.(fns{j});
                    if isnumeric(v) && size(v,1) == nSub && size(v,2) >= nROI
                        M = v(:,1:nROI);
                        return;
                    end
                end
            end
        end
    end
end

function addFitLine(x, y, col)
    % Adds a simple least-squares line if enough points
    if numel(x) < 4, return; end
    p = polyfit(x, y, 1);
    xx = linspace(min(x), max(x), 50);
    yy = polyval(p, xx);
    plot(xx, yy, '-', 'LineWidth', 2, 'Color', col);
end
