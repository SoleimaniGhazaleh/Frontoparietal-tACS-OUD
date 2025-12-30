%% ============================================================
%  Boarder_and_Craving_Plot.m  (UPDATED FULL SCRIPT)
%
%  Update requested:
%   - Always SHOW blocks 1..8 on the x-axis for BOTH groups,
%     so you can see Neutral/Cue 7–8 even if some data are missing.
%   - Labels: Odd blocks = Neutral 1..4, Even blocks = Cue 1..4
%
%  Outputs:
%   1) RT to yellow border (BOX_RESPONSE, event_code=5) per block:
%        - BLOCK_BoxRT_long_bySubject.csv
%        - Figure: Sham (A) bar+scatter by block (Neutral 1, Cue 1, ...)
%        - Figure: Active (B) bar+scatter by block
%
%   2) Craving ratings (RATING_RESPONSE, event_code=7) per block:
%        - BLOCK_CravingRatings_long_bySubject.csv
%        - Figure: Sham (A) bar+scatter by block
%        - Figure: Active (B) bar+scatter by block
%
%  Notes:
%   - RT block assignment uses RATING_ONSET (event_code=6) as block ends
%   - Ratings use their sequential order (1..#ratings in a session)
% ============================================================

clear; clc;

%% ------------------ Paths ------------------
rootDir  = '/Volumes/ExtremeSSDD/LIBR_tACS/InsideScannerData/Events_Pre/Ekhtiari-2019-tACS';
taskStem = 'task-opioidcuereactivity2_events.tsv';

%% ------------------ Force number of blocks to display ------------------
N_BLOCKS_TO_SHOW = 8;  % Neutral 1..4 and Cue 1..4

%% ------------------ Colors ------------------
colA = hex2rgb('#F8766D'); % Sham
colB = hex2rgb('#00BFC4'); % Active

%% ------------------ Participant IDs ------------------
participant_ids_to_include = { ...
    'sub-AO136', 'sub-AO931', 'sub-AP385', 'sub-AQ090', 'sub-AQ822', 'sub-AS590', 'sub-AU880', 'sub-AY910', ...
    'sub-BH088', 'sub-BI343', 'sub-BI346', 'sub-BI360', 'sub-BI392', 'sub-BI449', 'sub-BI460', 'sub-BI461', ...
    'sub-BI491', 'sub-BI508', 'sub-BI513', 'sub-BI550', 'sub-BI565', 'sub-BI632', 'sub-BI634', 'sub-BI699', ...
    'sub-BI875', 'sub-BJ130', 'sub-BJ131', 'sub-BJ162', 'sub-BJ386', 'sub-BJ469', 'sub-BJ493', 'sub-BJ510', ...
    'sub-BJ542', 'sub-BJ601', 'sub-BJ643', 'sub-BJ660', 'sub-BL307', 'sub-BL501', 'sub-BL749', 'sub-BL779', ...
    'sub-BL884', 'sub-BM282', 'sub-BM419', 'sub-BN573', 'sub-BN643', 'sub-BP433', 'sub-BP446', 'sub-BP880', ...
    'sub-BP917', 'sub-BQ069', 'sub-BQ850', 'sub-BQ976', 'sub-BR275', 'sub-BR276', 'sub-BR496', 'sub-BR527', ...
    'sub-BR665', 'sub-BR952', 'sub-BR953', 'sub-XX361'};

%% ------------------ Group assignments (A=Sham, B=Active) ------------------
group_assignments = { ...
    'sub-AO136', 'B'; 'sub-AO931', 'B'; 'sub-AP385', 'B'; 'sub-AQ090', 'A'; 'sub-AQ822', 'B'; 'sub-AS590', 'B'; 'sub-AU880', 'B'; ...
    'sub-AY910', 'B'; 'sub-BH088', 'B'; 'sub-BI343', 'A'; 'sub-BI346', 'A'; 'sub-BI360', 'A'; 'sub-BI392', 'B'; 'sub-BI449', 'B'; ...
    'sub-BI460', 'A'; 'sub-BI461', 'B'; 'sub-BI491', 'A'; 'sub-BI508', 'B'; 'sub-BI513', 'A'; 'sub-BI550', 'B'; 'sub-BI565', 'A'; ...
    'sub-BI632', 'A'; 'sub-BI634', 'B'; 'sub-BI699', 'A'; 'sub-BI875', 'A'; 'sub-BJ130', 'B'; 'sub-BJ131', 'A'; 'sub-BJ162', 'B'; ...
    'sub-BJ386', 'A'; 'sub-BJ469', 'A'; 'sub-BJ493', 'B'; 'sub-BJ510', 'A'; 'sub-BJ542', 'B'; 'sub-BJ601', 'B'; 'sub-BJ643', 'A'; ...
    'sub-BJ660', 'A'; 'sub-BL307', 'A'; 'sub-BL501', 'A'; 'sub-BL749', 'A'; 'sub-BL779', 'B'; 'sub-BL884', 'A'; 'sub-BM282', 'B'; ...
    'sub-BM419', 'A'; 'sub-BN573', 'B'; 'sub-BN643', 'B'; 'sub-BP433', 'A'; 'sub-BP446', 'A'; 'sub-BP880', 'A'; 'sub-BP917', 'A'; ...
    'sub-BQ069', 'A'; 'sub-BQ850', 'B'; 'sub-BQ976', 'B'; 'sub-BR275', 'A'; 'sub-BR276', 'B'; 'sub-BR496', 'B'; 'sub-BR527', 'B'; ...
    'sub-BR665', 'B'; 'sub-BR952', 'A'; 'sub-BR953', 'A'; 'sub-XX361', 'A'};

groupMap = containers.Map(group_assignments(:,1), group_assignments(:,2)); % char 'A'/'B'

%% ------------------ Output tables ------------------
BlockRT_long  = table();   % SubjectID, SessionID, Group, ConditionLabel, BlockIndex, MeanBoxRT, Ntrials
Craving_long  = table();   % SubjectID, SessionID, Group, ConditionLabel, BlockIndex, Rating, RatingOnset
missingFiles  = {};

%% ============================================================
%  Loop subjects and sessions
%% ============================================================
for i = 1:numel(participant_ids_to_include)
    sid = participant_ids_to_include{i};

    if ~isKey(groupMap, sid)
        warning('No group assignment for %s (skipping)', sid);
        continue
    end

    grpChar = groupMap(sid); % 'A' or 'B'
    if strcmp(grpChar,'A')
        condLabel = 'Sham';
    else
        condLabel = 'Active';
    end

    subDir = fullfile(rootDir, sid);
    sesList = dir(fullfile(subDir, 'ses-*'));
    sesList = sesList([sesList.isdir]);

    if isempty(sesList)
        warning('No sessions found for %s', sid);
        continue
    end

    for s = 1:numel(sesList)
        ses = sesList(s).name;  % e.g., 'ses-v1'
        tsvPath = fullfile(subDir, ses, 'func', sprintf('%s_%s_%s', sid, ses, taskStem));

        if ~isfile(tsvPath)
            missingFiles{end+1,1} = tsvPath; %#ok<AGROW>
            continue
        end

        T = readtable(tsvPath, 'FileType','text', 'Delimiter','\t');
        T = harmonizeTypes(T);

        % --------------------------
        % Define block ends by RATING_ONSET (event_code=6)
        % --------------------------
        ratingOn = T(T.event_code == 6, :);
        ratingTimes = ratingOn.onset;
        ratingTimes = ratingTimes(isfinite(ratingTimes));
        ratingTimes = sort(ratingTimes);

        % We'll use as many blocks as exist, but cap/plot to N_BLOCKS_TO_SHOW
        nBlocks = min(numel(ratingTimes), N_BLOCKS_TO_SHOW);

        % --------------------------
        % (1) RT to yellow border: BOX_RESPONSE (event_code=5)
        % --------------------------
        Box = T(T.event_code == 5, :);

        % NOTE: We do NOT require result==1; this avoids dropping valid RTs
        boxOn = Box.onset;
        boxRT = toNumericNA(Box.response_time);

        good = isfinite(boxOn) & isfinite(boxRT) & boxRT >= 0;
        boxOn = boxOn(good);
        boxRT = boxRT(good);

        if ~isempty(boxRT) && nBlocks > 0
            edges = [-Inf; ratingTimes(1:nBlocks)];
            blockIdx = discretize(boxOn, edges); % 1..nBlocks (others NaN)

            for b = 1:nBlocks
                inB = (blockIdx == b);
                if ~any(inB), continue; end

                tmp = table();
                tmp.SubjectID      = {sid};
                tmp.SessionID      = {ses};
                tmp.Group          = {grpChar};
                tmp.ConditionLabel = {condLabel};
                tmp.BlockIndex     = b;
                tmp.MeanBoxRT      = mean(boxRT(inB), 'omitnan');
                tmp.Ntrials        = sum(inB);

                BlockRT_long = [BlockRT_long; tmp]; %#ok<AGROW>
            end
        end

        % --------------------------
        % (2) Craving rating: RATING_RESPONSE (event_code=7)
        % --------------------------
        R = T(T.event_code == 7, :);

        rVals = toNumericNA(R.response);
        rOn   = R.onset;

        validR = isfinite(rVals) & ismember(rVals, [1 2 3 4]) & isfinite(rOn);
        rVals = rVals(validR);
        rOn   = rOn(validR);

        % Keep only first N_BLOCKS_TO_SHOW ratings (if extra exist)
        nR = min(numel(rVals), N_BLOCKS_TO_SHOW);

        for b = 1:nR
            tmp = table();
            tmp.SubjectID      = {sid};
            tmp.SessionID      = {ses};
            tmp.Group          = {grpChar};
            tmp.ConditionLabel = {condLabel};
            tmp.BlockIndex     = b;
            tmp.Rating         = rVals(b);
            tmp.RatingOnset    = rOn(b);

            Craving_long = [Craving_long; tmp]; %#ok<AGROW>
        end
    end
end

%% ------------------ Save tables ------------------
outRT   = fullfile(rootDir, 'BLOCK_BoxRT_long_bySubject.csv');
outRate = fullfile(rootDir, 'BLOCK_CravingRatings_long_bySubject.csv');
writetable(BlockRT_long, outRT);
writetable(Craving_long, outRate);

fprintf('\nSaved:\n%s\n%s\n', outRT, outRate);

if ~isempty(missingFiles)
    fprintf('\nMissing %d TSV files (showing up to 10):\n', numel(missingFiles));
    disp(missingFiles(1:min(10,end)));
end

%% ------------------ Plot RT (bar+scatter) ------------------
plot_group_boxrt(BlockRT_long, 'A', 'Sham',   colA, N_BLOCKS_TO_SHOW);
plot_group_boxrt(BlockRT_long, 'B', 'Active', colB, N_BLOCKS_TO_SHOW);

%% ------------------ Plot Craving ratings (bar+scatter) ------------------
plot_group_rating(Craving_long, 'A', 'Sham',   colA, N_BLOCKS_TO_SHOW);
plot_group_rating(Craving_long, 'B', 'Active', colB, N_BLOCKS_TO_SHOW);

%% ================== Local helper functions ==================

function plot_group_boxrt(BlockRT_long, grpLetter, grpName, col, NblocksShow)

if isempty(BlockRT_long)
    warning('No RT data to plot.');
    return
end

isG = strcmp(BlockRT_long.Group, grpLetter);
D = BlockRT_long(isG, :);
if isempty(D)
    warning('No RT data for group %s', grpName);
    return
end

% Average within Subject × Block (handles multiple sessions)
D2 = groupsummary(D, {'SubjectID','BlockIndex'}, 'mean', 'MeanBoxRT');
valVar = 'mean_MeanBoxRT';

blocks = (1:NblocksShow)';

[ybar, sem, labels] = block_mean_sem_labels(D2, blocks, valVar);

figure('Color','w'); hold on;
bar(blocks, ybar, 'FaceColor', col, 'FaceAlpha', 0.60, 'EdgeColor', 'none');
errorbar(blocks, ybar, sem, 'k.', 'LineWidth', 1);

for iB = 1:numel(blocks)
    b = blocks(iB);
    if any(D2.BlockIndex == b)
        vals = D2.(valVar)(D2.BlockIndex == b);
        vals = vals(isfinite(vals));
        x = b + (rand(size(vals)) - 0.5) * 0.18;
        scatter(x, vals, 30, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.90);
    end
end

xlabel('Block');
ylabel('RT to Yellow Border (s)');
title(sprintf('%s (Group %s): Box RT by Block', grpName, grpLetter));
xticks(blocks);
xticklabels(labels);
xtickangle(30);
ylim([0 5])
grid on; box on;

end

function plot_group_rating(Craving_long, grpLetter, grpName, col, NblocksShow)

if isempty(Craving_long)
    warning('No rating data to plot.');
    return
end

isG = strcmp(Craving_long.Group, grpLetter);
D = Craving_long(isG, :);
if isempty(D)
    warning('No rating data for group %s', grpName);
    return
end

% Average within Subject × Block (handles multiple sessions)
D2 = groupsummary(D, {'SubjectID','BlockIndex'}, 'mean', 'Rating');
valVar = 'mean_Rating';

blocks = (1:NblocksShow)';

[ybar, sem, labels] = block_mean_sem_labels(D2, blocks, valVar);

figure('Color','w'); hold on;
bar(blocks, ybar, 'FaceColor', col, 'FaceAlpha', 0.60, 'EdgeColor', 'none');
errorbar(blocks, ybar, sem, 'k.', 'LineWidth', 1);

for iB = 1:numel(blocks)
    b = blocks(iB);
    if any(D2.BlockIndex == b)
        vals = D2.(valVar)(D2.BlockIndex == b);
        vals = vals(isfinite(vals));
        x = b + (rand(size(vals)) - 0.5) * 0.18;
        scatter(x, vals, 30, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.90);
    end
end

xlabel('Block');
ylabel('Craving Rating (1–4)');
title(sprintf('%s (Group %s): Craving Rating by Block', grpName, grpLetter));
xticks(blocks);
xticklabels(labels);
xtickangle(30);
ylim([1 4]);
grid on; box on;

end

function [ybar, sem, labels] = block_mean_sem_labels(D2, blocks, valVar)

ybar = nan(size(blocks));
sem  = nan(size(blocks));

for iB = 1:numel(blocks)
    b = blocks(iB);
    vals = D2.(valVar)(D2.BlockIndex == b);
    vals = vals(isfinite(vals));
    if ~isempty(vals)
        ybar(iB) = mean(vals);
        sem(iB)  = std(vals) / sqrt(numel(vals));
    else
        ybar(iB) = NaN; % will show as missing bar
        sem(iB)  = NaN;
    end
end

labels = strings(size(blocks));
nNeutral = 0; nCue = 0;
for iB = 1:numel(blocks)
    if mod(blocks(iB),2)==1
        nNeutral = nNeutral + 1;
        labels(iB) = sprintf('Neutral %d', nNeutral);
    else
        nCue = nCue + 1;
        labels(iB) = sprintf('Cue %d', nCue);
    end
end

end

function T = harmonizeTypes(T)

req = {'onset','duration','trial_type','trial_number','event_code','response_time','response','result'};
for k = 1:numel(req)
    if ~ismember(req{k}, T.Properties.VariableNames)
        error('Missing column "%s" in a TSV file.', req{k});
    end
end

T.onset        = toNumericNA(T.onset);
T.duration     = toNumericNA(T.duration);
T.trial_number = toNumericNA(T.trial_number);
T.event_code   = toNumericNA(T.event_code);

T.response_time = toCellStr(T.response_time);
T.response      = toCellStr(T.response);
T.result        = toCellStr(T.result);
T.trial_type    = toCellStr(T.trial_type);

end

function c = toCellStr(v)

if iscell(v)
    c = cell(size(v));
    for i = 1:numel(v)
        if isempty(v{i})
            c{i} = '';
        elseif isnumeric(v{i})
            c{i} = num2str(v{i});
        else
            c{i} = char(v{i});
        end
    end
    return
end

if isnumeric(v)
    c = arrayfun(@(x) num2str(x), v, 'UniformOutput', false);
    return
end

c = cellstr(v);
end

function x = toNumericNA(v)

if isempty(v)
    x = NaN;
    return
end

if isnumeric(v)
    x = double(v);
    return
end

c = toCellStr(v);
x = NaN(size(c));
for i = 1:numel(c)
    s = strtrim(c{i});
    if isempty(s) || strcmpi(s,'NA')
        x(i) = NaN;
    else
        x(i) = str2double(s);
    end
end
end

function rgb = hex2rgb(hex)

if hex(1) == '#'
    hex = hex(2:end);
end
if numel(hex) ~= 6
    error('hex2rgb expects a 6-digit hex color, e.g., #F8766D');
end

r = hex2dec(hex(1:2));
g = hex2dec(hex(3:4));
b = hex2dec(hex(5:6));

rgb = [r g b] / 255;
end
