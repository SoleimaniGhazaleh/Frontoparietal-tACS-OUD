%% ============================================================
%  FULL SCRIPT: PLUS Questionnaire Visualization (Abstinence / Last Use)
%  - Load participant_last_use_summary_plus.tsv
%  - Keep only participants in participant_ids_to_include
%  - Keep only session == "before_pre_fmri_arm_1"
%  - Add Group A (Sham) vs Group B (Active) from group_assignments
%  - Numeric vars: bar (mean±SEM) + scatter (colored by group)
%  - Binary vars (plus_1a/2a/4a/6a): percent "Yes" + scatter (0/100)
%  - Text vars (plus_5, plus_7, *_otherdesc): frequency (top N + Other) by group
%  - Adds human-readable topic labels based on PLUS JSON descriptions
%  - Saves PNG + PDF to <phenotype>/PLUS_Plots/
%
%  Ghazaleh - Updated Dec 2025
%% ============================================================

clear; clc; close all;

%% ------------------ Paths ------------------
tsvPath = "/Volumes/ExtremeSSDD/LIBR_tACS/Demographic/phenotype/participant_last_use_summary_plus.tsv";
outDir  = '/Users/ghazaleh/Documents/1.PROJECTS/PUBLICATIONS/MolecularPsychiatry/MolecularPsy/Review_round1/Abstinence/PLUS_Plots';
if ~exist(outDir, "dir"); mkdir(outDir); end

%% ------------------ Session to use ------------------
targetSession = "before_pre_fmri_arm_1";

%% ------------------ Colors ------------------
colA = hex2rgb('#F8766D'); % Sham (Group A)
colB = hex2rgb('#00BFC4'); % Active (Group B)

%% ------------------ Participants to include ------------------
participant_ids_to_include = { ...
    'sub-AO136', 'sub-AO931', 'sub-AP385', 'sub-AQ090', 'sub-AQ822', 'sub-AS590', 'sub-AU880', 'sub-AY910', ...
    'sub-BH088', 'sub-BI343', 'sub-BI346', 'sub-BI360', 'sub-BI392', 'sub-BI449', 'sub-BI460', 'sub-BI461', ...
    'sub-BI491', 'sub-BI508', 'sub-BI513', 'sub-BI550', 'sub-BI565', 'sub-BI632', 'sub-BI634', 'sub-BI699', ...
    'sub-BI875', 'sub-BJ130', 'sub-BJ131', 'sub-BJ162', 'sub-BJ386', 'sub-BJ469', 'sub-BJ493', 'sub-BJ510', ...
    'sub-BJ542', 'sub-BJ601', 'sub-BJ643', 'sub-BJ660', 'sub-BL307', 'sub-BL501', 'sub-BL749', 'sub-BL779', ...
    'sub-BL884', 'sub-BM282', 'sub-BM419', 'sub-BN573', 'sub-BN643', 'sub-BP433', 'sub-BP446', 'sub-BP880', ...
    'sub-BP917', 'sub-BQ069', 'sub-BQ850', 'sub-BQ976', 'sub-BR275', 'sub-BR276', 'sub-BR496', 'sub-BR527', ...
    'sub-BR665', 'sub-BR952', 'sub-BR953', 'sub-XX361'};

%% ------------------ Group assignments ------------------
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

Gmap = containers.Map(group_assignments(:,1), group_assignments(:,2));

%% ------------------ Question labels (topics) from PLUS JSON ------------------
Qlabel = containers.Map;
Qlabel('plus_1')  = 'Days since last marijuana use';
Qlabel('plus_1a') = 'Ever used marijuana (Yes/No)';
Qlabel('plus_2')  = 'Days since last alcohol use';
Qlabel('plus_2a') = 'Ever drank alcohol (Yes/No)';
Qlabel('plus_3')  = 'Days since last alcohol binge';
Qlabel('plus_4')  = 'Days since last illicit substance use';
Qlabel('plus_4a') = 'Ever used illicit substance (Yes/No)';
Qlabel('plus_5')  = 'Last used illicit substance (name)';
Qlabel('plus_5otherdesc') = 'Other illicit substance (description)';
Qlabel('plus_6')  = 'Days since last medication use';
Qlabel('plus_6a') = 'Medication in past 2 weeks (Yes/No)';
Qlabel('plus_7')  = 'Medications used in past 2 weeks (free text)';
Qlabel('plus_7otherdesc') = 'Other medication (description)';
Qlabel('plus_8')  = 'Caffeine consumed today (mg)';
Qlabel('plus_9')  = 'Hours since last caffeine beverage';
Qlabel('plus_10') = 'Times tobacco used today';
Qlabel('plus_11') = 'Hours since last tobacco use';
Qlabel('plus_assessed_at') = 'Assessed at (date/time)';

%% ------------------ Load TSV ------------------
opts = detectImportOptions(tsvPath, "FileType","text", "Delimiter","\t");
% Make sure key columns are strings if possible
try
    opts = setvartype(opts, 'participant_id', 'string');
    opts = setvartype(opts, 'session', 'string');
catch
    % older MATLAB may not have setvartype; that's OK
end

T = readtable(tsvPath, opts);

% Ensure participant_id and session exist
assert(ismember('participant_id', T.Properties.VariableNames), "Missing participant_id column.");
assert(ismember('session', T.Properties.VariableNames), "Missing session column.");

% Convert to string robustly
T.participant_id = stringifyColumn(T.participant_id);
T.session        = stringifyColumn(T.session);

% Keep only target session
T = T(T.session == targetSession, :);

% Keep only included participants
incl = ismember(cellstr(T.participant_id), participant_ids_to_include);
T = T(incl,:);

% Add Group
grp = strings(height(T),1);
for i = 1:height(T)
    pid = char(T.participant_id(i));
    if isKey(Gmap, pid)
        grp(i) = string(Gmap(pid));
    else
        grp(i) = "NA";
    end
end
T.Group = grp;

% Keep only A/B rows
T = T(T.Group=="A" | T.Group=="B", :);

fprintf("Loaded %d rows for session %s (Sham A=%d, Active B=%d)\n", ...
    height(T), targetSession, sum(T.Group=="A"), sum(T.Group=="B"));

%% ------------------ Identify variables ------------------
metaCols = ["participant_id","session","Group"];
vars  = string(T.Properties.VariableNames);
qvars = setdiff(vars, metaCols, 'stable');

% Text variables (explicit)
textVars = intersect(["plus_5","plus_5otherdesc","plus_7","plus_7otherdesc","plus_assessed_at"], qvars, 'stable');

% Binary yes/no variables (explicit)
binaryVars = intersect(["plus_1a","plus_2a","plus_4a","plus_6a"], qvars, 'stable');

% Numeric variables = everything else (including binaries, but we treat them separately)
numVars = setdiff(qvars, textVars, 'stable');

% Coerce numeric vars to double (robust)
for v = numVars
    if ~isnumeric(T.(v))
        T.(v) = coerceToDouble(T.(v));
    end
end

%% ------------------ Plot numeric variables ------------------
for v = numVars

    topic = char(v);
    niceTitle = topic;
    if isKey(Qlabel, topic), niceTitle = Qlabel(topic); else, niceTitle = strrep(topic,'_',' '); end

    xA = T.(v)(T.Group=="A");
    xB = T.(v)(T.Group=="B");

    xA = xA(isfinite(xA));
    xB = xB(isfinite(xB));

    if isempty(xA) && isempty(xB)
        continue;
    end

if any(v == binaryVars)
    plotYesNoCounts(T, v, niceTitle, targetSession, outDir, colA, colB);
    continue;
end

    % ---- Continuous numeric handling ----
    mA = mean(xA,'omitnan');  sA = std(xA,'omitnan')/sqrt(max(1,numel(xA)));
    mB = mean(xB,'omitnan');  sB = std(xB,'omitnan')/sqrt(max(1,numel(xB)));

    fig = figure('Color','w','Position',[100 100 200 200]); hold on;

    % Bars with fixed colors (robust)
    bar(1, mA, 0.6, 'FaceColor', colA, 'EdgeColor','none');
    bar(2, mB, 0.6, 'FaceColor', colB, 'EdgeColor','none');

    % Error bars
    errorbar([1 2], [mA mB], [sA sB], 'k', 'LineStyle','none', 'LineWidth',1.5);

    % Scatter with jitter colored by group
    jitterA = (rand(size(xA))-0.5)*0.18;
    jitterB = (rand(size(xB))-0.5)*0.18;
    scatter(1 + jitterA, xA, 42, 'filled', ...
        'MarkerFaceColor', colA, 'MarkerFaceAlpha',0.85, ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);
    scatter(2 + jitterB, xB, 42, 'filled', ...
        'MarkerFaceColor', colB, 'MarkerFaceAlpha',0.85, ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);

    set(gca,'XTick',[1 2],'XTickLabel',{'Sham (A)','Active (B)'});
    ylabel(niceTitle);
    title(niceTitle);
    grid on; box off;

    % Reasonable y-lims
    allx = [xA; xB];
    if ~isempty(allx)
        pad = 0.05 * range(allx);
        if pad==0, pad = 1; end
        ylim([min(allx)-pad, max(allx)+pad]);
    end

    fbase = fullfile(outDir, "NUM_" + v + "_" + targetSession);
    saveas(fig, fbase + ".png");
    %saveas(fig, fbase + ".pdf");
    close(fig);
end

%% ------------------ Plot TEXT variables (count bar plots only) ------------------
topN = 10;   % number of categories to show per group

for v = textVars

    topic = char(v);
    if isKey(Qlabel, topic)
        niceTitle = Qlabel(topic);
    else
        niceTitle = strrep(topic,'_',' ');
    end

    % Convert column to string
    s = stringifyColumn(T.(v));

    % Split by group
    sA = cleanTextEntries(s(T.Group=="A"));
    sB = cleanTextEntries(s(T.Group=="B"));

    if all(sA=="") && all(sB=="")
        continue;
    end

    % Split comma/semicolon-separated lists if needed
    splitMode = (v=="plus_5" | v=="plus_5otherdesc" | ...
                 v=="plus_7" | v=="plus_7otherdesc");
    if splitMode
        sA = splitTokens(sA);
        sB = splitTokens(sB);
    end

    % Count categories
    [catsA, cntA] = topCounts(sA, topN);
    [catsB, cntB] = topCounts(sB, topN);

    fig = figure('Color','w','Position',[100 100 1000 450]);

    % -------- Sham (A) --------
    subplot(1,2,1); hold on;
    bhA = barh(cntA, 0.75);
    bhA.FaceColor = colA;
    bhA.EdgeColor = 'none';
    set(gca,'YDir','reverse');
    yticks(1:numel(catsA));
    yticklabels(catsA);
    xlabel('Count');
    title(niceTitle + " — Sham (A)");
    grid on; box off;

    % -------- Active (B) --------
    subplot(1,2,2); hold on;
    bhB = barh(cntB, 0.75);
    bhB.FaceColor = colB;
    bhB.EdgeColor = 'none';
    set(gca,'YDir','reverse');
    yticks(1:numel(catsB));
    yticklabels(catsB);
    xlabel('Count');
    title(niceTitle + " — Active (B)");
    grid on; box off;

    sgtitle(niceTitle + " (" + targetSession + ")");

    % Save
    fbase = fullfile(outDir, "TXT_" + v + "_" + targetSession);
    saveas(fig, fbase + ".png");
    %saveas(fig, fbase + ".pdf");
    close(fig);
end

fprintf("Done. Figures saved to:\n%s\n", outDir);



%% ============================================================
% Helper functions
%% ============================================================

function rgb = hex2rgb(hex)
hex = char(hex);
hex = strrep(hex,'#','');
rgb = [hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))] / 255;
end

function x = coerceToDouble(col)
% Convert table column to double robustly
if isnumeric(col)
    x = col;
    return;
end

% Convert to string first
s = stringifyColumn(col);
% Common non-numeric tokens to empty
s = strtrim(s);
bad = (s=="" | lower(s)=="na" | lower(s)=="n/a" | lower(s)=="nan" | lower(s)=="none");
s(bad) = "";

x = nan(size(s));
for i = 1:numel(s)
    if s(i)==""; continue; end
    x(i) = str2double(s(i));
end
end

function s = stringifyColumn(col)
% Convert almost anything to string array
if isstring(col)
    s = col;
elseif iscell(col)
    s = string(col);
elseif isnumeric(col)
    s = string(col);
elseif ischar(col)
    s = string(cellstr(col));
else
    s = string(col);
end
s = reshape(s, [], 1);
end

function s = cleanTextEntries(s)
% Normalize empties and common NA strings; remove basic punctuation
s = strtrim(lower(s));
s = regexprep(s, '[\.\(\)\[\]"]', '');
bad = (s=="" | s=="na" | s=="n/a" | s=="nan" | s=="none" | s=="null");
s(bad) = "";
end

function tokens = splitTokens(s)
% MATLAB-version-safe tokenization using regexprep + regexp split
tokens = strings(0,1);

for i = 1:numel(s)
    if s(i) == ""; continue; end
    tmp = char(s(i));

    % Turn common separators into commas
    tmp = regexprep(tmp, '[;/\|]+', ',');
    % Normalize spaces
    tmp = regexprep(tmp, '\s+', ' ');
    % Split on commas
    parts = strtrim(regexp(tmp, ',', 'split')); % cell array of char
    parts = string(parts);
    parts = strtrim(parts);
    parts = parts(parts ~= "");

    tokens = [tokens; parts(:)]; %#ok<AGROW>
end

% Final cleanup: remove trivial tokens
tokens = cleanTextEntries(tokens);
tokens = tokens(tokens~="");
end

function [cats, cnt] = topCounts(s, topN)
% Return topN categories and counts, with "Other" pooled
s = s(s~="");
if isempty(s)
    cats = "No data";
    cnt  = 0;
    return;
end

[catsAll, ~, idx] = unique(s);
cntAll = accumarray(idx, 1);

[sortedCnt, order] = sort(cntAll, 'descend');
sortedCats = catsAll(order);

k = min(topN, numel(sortedCats));
cats = sortedCats(1:k);
cnt  = sortedCnt(1:k);

if numel(sortedCats) > k
    otherCnt = sum(sortedCnt(k+1:end));
    cats = [cats; "Other"];
    cnt  = [cnt; otherCnt];
end

% Make ytick labels readable
cats = replace(cats, "_", " ");
end

function plotYesNoCounts(T, v, niceTitle, targetSession, outDir, colA, colB)
% Plot counts of Yes/No per group for a binary variable v

xA = T.(v)(T.Group=="A");
xB = T.(v)(T.Group=="B");

% Keep finite
xA = xA(isfinite(xA));
xB = xB(isfinite(xB));

% Interpret: nonzero = Yes, zero = No
yesA = sum(xA ~= 0);
noA  = sum(xA == 0);

yesB = sum(xB ~= 0);
noB  = sum(xB == 0);

fig = figure('Color','w','Position',[100 100 200 200]); hold on;

% Stacked bars: [No Yes]
Y = [noA yesA; noB yesB];

hb = bar(Y, 'stacked', 'BarWidth', 0.65);
% Color stacks: No=light gray, Yes=group color
hb(1).FaceColor = [0.85 0.85 0.85];  % No
hb(2).FaceColor = [0.25 0.25 0.25];  % Yes (will recolor by group below)

% To color "Yes" separately for each group, we overlay a second bar:
cla; hold on;
% Base: No (gray) stacked with Yes (group-colored) using two stacked bars per group
bar(1, noA, 0.65, 'FaceColor',[0.85 0.85 0.85], 'EdgeColor','none');
bar(1, yesA,0.65, 'FaceColor',colA,             'EdgeColor','none', 'BaseValue', noA);

bar(2, noB, 0.65, 'FaceColor',[0.85 0.85 0.85], 'EdgeColor','none');
bar(2, yesB,0.65, 'FaceColor',colB,             'EdgeColor','none', 'BaseValue', noB);

set(gca,'XTick',[1 2],'XTickLabel',{'Sham (A)','Active (B)'});
ylabel('Count');
title(niceTitle + " (" + targetSession + ")");
legend({'No','Yes'}, 'Location','northeast');
grid on; box off;

% Add labels on top of stacks
text(1, noA/2, sprintf('No=%d', noA), 'HorizontalAlignment','center');
text(1, noA + yesA/2, sprintf('Yes=%d', yesA), 'HorizontalAlignment','center', 'Color','w');

text(2, noB/2, sprintf('No=%d', noB), 'HorizontalAlignment','center');
text(2, noB + yesB/2, sprintf('Yes=%d', yesB), 'HorizontalAlignment','center', 'Color','w');

% Save
fbase = fullfile(outDir, "BIN_COUNTS_" + string(v) + "_" + targetSession);
saveas(fig, fbase + ".png");
%saveas(fig, fbase + ".pdf");
close(fig);
end
