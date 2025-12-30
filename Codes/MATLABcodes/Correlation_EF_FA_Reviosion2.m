%% ============================================================
%  Top-5 ROIs (211–230) by highest E-field (Group B)
%  and report ΔFunctional Activity (Post - Pre) in those ROIs
%% ============================================================

clear; clc; close all;

%% ------------------ INPUTS ------------------
load('/Volumes/ExtremeSSDD/LIBR_tACS/Final_AB_data/EF.mat');   % expects EF_B or EF_A or EF

preFile  = '/Users/ghazaleh/Documents/1.PROJECTS/PUBLICATIONS/MolecularPsychiatry/MolecularPsy/Review_round1/BNA_Activation/Final60_Pre.xlsx';
postFile = '/Users/ghazaleh/Documents/1.PROJECTS/PUBLICATIONS/MolecularPsychiatry/MolecularPsy/Review_round1/BNA_Activation/Final60_Post.xlsx';

roiList = 211:230;
topK    = 5;
GROUP   = 'B';   % 'B' Active, 'A' Sham

%% ------------------ GET EF MATRIX ------------------
if exist('EF_B','var') && strcmp(GROUP,'B')
    EFmat = EF_B;
elseif exist('EF_A','var') && strcmp(GROUP,'A')
    EFmat = EF_A;
elseif exist('EF','var')
    EFmat = EF;
else
    error('Cannot find EF matrix. Expected EF_B / EF_A / EF in EF.mat');
end

EF_roi = EFmat(:, roiList);

%% ------------------ LOAD FUNCTIONAL ACTIVITY ------------------
Tpre  = readtable(preFile);
Tpost = readtable(postFile);

% Build column names as cell arrays: {'Pre_211','Pre_212',...}
preVars  = arrayfun(@(r) sprintf('Pre_%d',  r), roiList, 'UniformOutput', false);
postVars = arrayfun(@(r) sprintf('Post_%d', r), roiList, 'UniformOutput', false);

% Select Group B rows if Group column exists
if any(strcmpi(Tpre.Properties.VariableNames,'Group')) && any(strcmpi(Tpost.Properties.VariableNames,'Group'))
    gpre  = Tpre.Group;
    gpost = Tpost.Group;

    % Works whether Group is numeric, char, string, or categorical
    if isnumeric(gpre)
        idxG = (gpre == 2); % <-- only if your coding is 1/2 (change if needed)
        warning('Group is numeric; assuming B=2. If not true, change this line.');
    else
        idxG = strcmp(string(gpre), GROUP);
    end

    TpreG  = Tpre(idxG, :);
    TpostG = Tpost(idxG, :);
else
    warning('No "Group" column found; assuming the rows already correspond to Group B.');
    TpreG  = Tpre;
    TpostG = Tpost;
end

% Validate required columns exist
missingPre  = preVars(~ismember(preVars,  TpreG.Properties.VariableNames));
missingPost = postVars(~ismember(postVars, TpostG.Properties.VariableNames));
if ~isempty(missingPre)
    error('Missing Pre columns: %s', strjoin(missingPre, ', '));
end
if ~isempty(missingPost)
    error('Missing Post columns: %s', strjoin(missingPost, ', '));
end

Act_pre  = TpreG{:, preVars};
Act_post = TpostG{:, postVars};
dAct_roi = Act_post - Act_pre;

%% ------------------ ALIGN SUBJECT COUNTS ------------------
N = min(size(EF_roi,1), size(dAct_roi,1));
if size(EF_roi,1) ~= size(dAct_roi,1)
    warning('EF N=%d, Activity N=%d. Truncating to N=%d.', size(EF_roi,1), size(dAct_roi,1), N);
end
EF_roi   = EF_roi(1:N, :);
dAct_roi = dAct_roi(1:N, :);

%% ------------------ TOP-5 ROIs BY MEAN EF ------------------
meanEF = mean(EF_roi, 1, 'omitnan');
[~, order] = sort(meanEF, 'descend');

topIdxLocal = order(1:topK);
topROIs     = roiList(topIdxLocal);

%% ------------------ SUMMARIZE ΔACTIVITY IN TOP-5 ------------------
meanDAct = mean(dAct_roi(:, topIdxLocal), 1, 'omitnan');
sdDAct   = std(dAct_roi(:, topIdxLocal), 0, 'omitnan');
semDAct  = sdDAct ./ sqrt(sum(~isnan(dAct_roi(:, topIdxLocal)), 1));

Top5 = table(topROIs(:), meanEF(topIdxLocal)', meanDAct(:), semDAct(:), ...
    'VariableNames', {'ROI','MeanEF','MeanDeltaActivity','SEM_DeltaActivity'});

disp('===== Group B: Top-5 ROIs (211–230) by Mean E-field, with ΔActivity =====');
disp(Top5);

%% ------------------ OPTIONAL: SCATTERS (EF vs ΔActivity) ------------------
for k = 1:topK
    roi = topROIs(k);
    x = EF_roi(:, topIdxLocal(k));
    y = dAct_roi(:, topIdxLocal(k));

    ok = isfinite(x) & isfinite(y);
    r = NaN; p = NaN;
    if nnz(ok) >= 5
        [r,p] = corr(x(ok), y(ok), 'Type','Pearson');
    end

    figure('Color','w'); hold on;
    scatter(x, y, 55, 'filled');
    lsline;
    grid on;
    xlabel('Electric Field Strength');
    ylabel('\Delta Functional Activity (Post - Pre)');
    title(sprintf('ROI_%03d | Group B | r=%.2f, p=%.3g', roi, r, p));
end
