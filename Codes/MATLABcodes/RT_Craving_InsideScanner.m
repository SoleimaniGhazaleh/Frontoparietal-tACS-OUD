%% ============================================================
%  CueTask_RT_Craving_Time_Group_LME.m
%
%  Extract:
%    (1) RT to yellow border: BOX_RESPONSE (event_code=5) -> response_time (sec)
%    (2) Craving rating:      RATING_RESPONSE (event_code=7) -> response (1-4)
%
%  Adds factors:
%    - Condition: Neutral vs Cue (odd vs even blocks)
%    - Time: Pre vs Post
%    - Group: Sham (A) vs Active (B)
%
%  Outputs:
%    - Long RT table:      RT_long_Time.csv
%    - Long rating table:  Craving_long_Time.csv
%    - LME results printed to command window
%    - Figures for RT and Craving (each group separately)
% ============================================================

clear; clc;

%% ------------------ USER SETTINGS ------------------
taskStem = 'task-opioidcuereactivity2_events.tsv';
N_BLOCKS_TO_SHOW = 8; % Neutral1/Cue1 ... Neutral4/Cue4

% >>> IMPORTANT: Update these if needed <<<
rootPre  = '/Volumes/ExtremeSSDD/LIBR_tACS/InsideScannerData/Events_Pre/Ekhtiari-2019-tACS';
rootPost = '/Volumes/ExtremeSSDD/LIBR_tACS/InsideScannerData/Events_Post/Ekhtiari-2019-tACS'; % <-- change if your folder name differs

% Colors (ggplot2 default)
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

groupMap = containers.Map(group_assignments(:,1), group_assignments(:,2)); % 'A'/'B'

%% ------------------ Extract both times ------------------
RT_long      = table();   % SubjectID, Group, Time, BlockIndex, Condition, MeanBoxRT
Craving_long = table();   % SubjectID, Group, Time, BlockIndex, Condition, Rating

missingFiles = {};

times = {'Pre','Post'};
roots = {rootPre, rootPost};

for t = 1:2
    timeLabel = times{t};
    rootDir   = roots{t};

    for i = 1:numel(participant_ids_to_include)
        sid = participant_ids_to_include{i};

        if ~isKey(groupMap, sid)
            continue
        end
        grpChar = groupMap(sid); % 'A' or 'B'

        subDir = fullfile(rootDir, sid);
        sesList = dir(fullfile(subDir, 'ses-*'));
        sesList = sesList([sesList.isdir]);

        if isempty(sesList)
            % some datasets might have only ses-v1
            continue
        end

        for s = 1:numel(sesList)
            ses = sesList(s).name;
            tsvPath = fullfile(subDir, ses, 'func', sprintf('%s_%s_%s', sid, ses, taskStem));

            if ~isfile(tsvPath)
                missingFiles{end+1,1} = tsvPath; %#ok<AGROW>
                continue
            end

            T = readtable(tsvPath, 'FileType','text', 'Delimiter','\t');
            T = harmonizeTypes(T);

            % ------------------ Block ends (RATING_ONSET = 6) ------------------
            ratingOn = T(T.event_code==6,:);
            ratingTimes = ratingOn.onset;
            ratingTimes = ratingTimes(isfinite(ratingTimes));
            ratingTimes = sort(ratingTimes);

            nBlocks_here = min(numel(ratingTimes), N_BLOCKS_TO_SHOW);

            % ------------------ RT (BOX_RESPONSE = 5) ------------------
            Box = T(T.event_code==5,:);
            boxOn = Box.onset;
            boxRT = toNumericNA(Box.response_time);

            good = isfinite(boxOn) & isfinite(boxRT) & boxRT>=0;
            boxOn = boxOn(good);
            boxRT = boxRT(good);

            if ~isempty(boxRT) && nBlocks_here>0
                edges = [-Inf; ratingTimes(1:nBlocks_here)];
                blockIdx = discretize(boxOn, edges); % 1..nBlocks_here

                for b = 1:nBlocks_here
                    inB = (blockIdx==b);
                    if ~any(inB), continue; end

                    tmp = table();
                    tmp.SubjectID = {sid};
                    tmp.Group     = {grpChar};
                    tmp.Time      = {timeLabel};
                    tmp.BlockIndex= b;
                    tmp.Condition = {block_condition_label(b)}; % Neutral/Cue
                    tmp.MeanBoxRT = mean(boxRT(inB), 'omitnan');

                    RT_long = [RT_long; tmp]; %#ok<AGROW>
                end
            end

            % ------------------ Craving rating (RATING_RESPONSE = 7) ------------------
            R = T(T.event_code==7,:);
            rVals = toNumericNA(R.response);
            valid = isfinite(rVals) & ismember(rVals,[1 2 3 4]);

            rVals = rVals(valid);
            nR = min(numel(rVals), N_BLOCKS_TO_SHOW);

            for b = 1:nR
                tmp = table();
                tmp.SubjectID = {sid};
                tmp.Group     = {grpChar};
                tmp.Time      = {timeLabel};
                tmp.BlockIndex= b;
                tmp.Condition = {block_condition_label(b)}; % Neutral/Cue
                tmp.Rating    = rVals(b);

                Craving_long = [Craving_long; tmp]; %#ok<AGROW>
            end
        end
    end
end

%% ------------------ Save extracted tables ------------------
outDir = rootPre; % save next to Pre directory (adjust if you want)
writetable(RT_long,      fullfile(outDir, 'RT_long_Time.csv'));
writetable(Craving_long, fullfile(outDir, 'Craving_long_Time.csv'));

fprintf('\nSaved:\n  %s\n  %s\n', ...
    fullfile(outDir,'RT_long_Time.csv'), fullfile(outDir,'Craving_long_Time.csv'));

if ~isempty(missingFiles)
    fprintf('\nMissing %d TSV files (showing up to 10):\n', numel(missingFiles));
    disp(missingFiles(1:min(10,end)));
end

%% ============================================================
%  Prepare subject-level summaries for LME:
%   Use Cue vs Neutral within each Time, so Condition is a factor.
%% ============================================================

% Average within Subject × Time × Condition (collapse blocks)
RT_sub = groupsummary(RT_long, {'SubjectID','Group','Time','Condition'}, 'mean', 'MeanBoxRT');
RT_sub.Properties.VariableNames{end} = 'RT'; % last col is mean_MeanBoxRT (rename safely)

Cr_sub = groupsummary(Craving_long, {'SubjectID','Group','Time','Condition'}, 'mean', 'Rating');
Cr_sub.Properties.VariableNames{end} = 'Craving';

% Convert to categorical for fitlme
RT_sub.SubjectID = categorical(RT_sub.SubjectID);
RT_sub.Group     = categorical(RT_sub.Group);   % 'A'/'B'
RT_sub.Time      = categorical(RT_sub.Time);    % Pre/Post
RT_sub.Condition = categorical(RT_sub.Condition); % Neutral/Cue

Cr_sub.SubjectID = categorical(Cr_sub.SubjectID);
Cr_sub.Group     = categorical(Cr_sub.Group);
Cr_sub.Time      = categorical(Cr_sub.Time);
Cr_sub.Condition = categorical(Cr_sub.Condition);

%% ------------------ LME: RT (seconds) ------------------
% Optional: log-transform if skewed (comment/uncomment)
% RT_sub.logRT = log(RT_sub.RT);
% lmeRT = fitlme(RT_sub, 'logRT ~ Group*Time*Condition + (1|SubjectID)');
lmeRT = fitlme(RT_sub, 'RT ~ Group*Time*Condition + (1|SubjectID)');

disp('================ LME: RT (seconds) ================');
disp(anova(lmeRT));
disp(lmeRT);

%% ------------------ LME: Craving (Likert 1–4) ------------------
% Treat as continuous for primary analysis (common); ordinal model would require other tools.
lmeCr = fitlme(Cr_sub, 'Craving ~ Group*Time*Condition + (1|SubjectID)');

disp('================ LME: Craving (Likert 1–4) ================');
disp(anova(lmeCr));
disp(lmeCr);

%% ============================================================
%  Visualization (bar + scatter) for each group:
%   - RT and Craving across (Time × Condition)
%% ============================================================

plot_time_condition(RT_sub, 'RT',      'Reaction time (s)',        colA, colB);
plot_time_condition(Cr_sub, 'Craving', 'Craving rating (1–4)',     colA, colB, [1 4]);

%% ============================ FUNCTIONS ============================

function lbl = block_condition_label(b)
% Odd blocks = Neutral, Even blocks = Cue
if mod(b,2)==1
    lbl = 'Neutral';
else
    lbl = 'Cue';
end
end

function plot_time_condition(Tsub, yvar, ylab, colA, colB, ylims)
% Creates two figures (Sham/Active) showing Pre/Post × Neutral/Cue bars + scatter.
if nargin < 6, ylims = [0 5]; end

groups = {'A','B'};
groupNames = {'Sham','Active'};
cols = {colA, colB};

for gi = 1:2
    g = groups{gi};
    D = Tsub(Tsub.Group==g, :);
    if isempty(D), continue; end

    % Ensure order: Pre then Post; Neutral then Cue
    times = categories(D.Time); %#ok<NASGU>
    conds = categories(D.Condition); %#ok<NASGU>

    % Force desired plotting order:
    timeOrder = {'Pre','Post'};
    condOrder = {'Neutral','Cue'};

    % Compute mean/SEM per cell using subject-level values
    means = nan(2,2);
    sems  = nan(2,2);

    for ti = 1:2
        for ci = 1:2
            sel = (D.Time==timeOrder{ti}) & (D.Condition==condOrder{ci});
            vals = D.(yvar)(sel);
            vals = vals(isfinite(vals));
            if ~isempty(vals)
                means(ti,ci) = mean(vals);
                sems(ti,ci)  = std(vals)/sqrt(numel(vals));
            end
        end
    end

    figure('Color','w'); hold on;

    % x positions: [Pre-Neutral, Pre-Cue, Post-Neutral, Post-Cue]
    x = [1 2 4 5];
    y = [means(1,1) means(1,2) means(2,1) means(2,2)];
    e = [sems(1,1)  sems(1,2)  sems(2,1)  sems(2,2)];

    bar(x, y, 'FaceColor', cols{gi}, 'FaceAlpha', 0.60, 'EdgeColor','none');
    errorbar(x, y, e, 'k.', 'LineWidth', 1);

    % Scatter: one dot per subject per cell
    % (We already have one row per Subject×Time×Condition in Tsub)
    for xi = 1:numel(x)
        if xi==1, sel = (D.Time=='Pre')  & (D.Condition=='Neutral'); end
        if xi==2, sel = (D.Time=='Pre')  & (D.Condition=='Cue');     end
        if xi==3, sel = (D.Time=='Post') & (D.Condition=='Neutral'); end
        if xi==4, sel = (D.Time=='Post') & (D.Condition=='Cue');     end
        vals = D.(yvar)(sel);
        vals = vals(isfinite(vals));
        xx = x(xi) + (rand(size(vals))-0.5)*0.18;
        scatter(xx, vals, 30, 'MarkerFaceColor', cols{gi}, 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.9);
    end

    xticks(x);
    xticklabels({'Pre Neutral','Pre Cue','Post Neutral','Post Cue'});
    xtickangle(20);
    ylabel(ylab);
    title(sprintf('%s: %s by Time × Condition', groupNames{gi}, yvar));
    grid on; box on;

    if ~isempty(ylims)
        ylim(ylims);
    end
end
end

function T = harmonizeTypes(T)
req = {'onset','duration','trial_type','trial_number','event_code','response_time','response','result'};
for k = 1:numel(req)
    if ~ismember(req{k}, T.Properties.VariableNames)
        error('Missing column "%s" in TSV.', req{k});
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
if isempty(v), x = NaN; return; end
if isnumeric(v), x = double(v); return; end
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
if hex(1)=='#', hex = hex(2:end); end
r = hex2dec(hex(1:2)); g = hex2dec(hex(3:4)); b = hex2dec(hex(5:6));
rgb = [r g b]/255;
end
