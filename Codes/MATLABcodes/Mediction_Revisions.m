%% ============================================================
%  Medication Use History (before_pre_fmri_arm_1) - Group Plots
%  - Reads medical_history.tsv
%  - Filters to participant list + session
%  - Adds Group (A=Sham, B=Active)
%  - Visualizations (ALL color-consistent):
%       Fig1: mean±SEM + scatter (count)
%       Fig2: histograms (count)
%       Fig3: 0/1/2/3+ burden (stacked; group-colored shades)
%       Fig4: top medication names (barh)
%% ============================================================

clear; clc; close all;

%% ------------------ Path ------------------
tsvPath = '/Volumes/ExtremeSSDD/LIBR_tACS/Demographic/phenotype/medical_history.tsv';

%% ------------------ Colors ------------------
colA = hex2rgb('#F8766D'); % Sham (Group A)
colB = hex2rgb('#00BFC4'); % Active (Group B)

%% ------------------ Participant IDs ------------------
participant_ids_to_include = { ...
    'sub-AO136','sub-AO931','sub-AP385','sub-AQ090','sub-AQ822','sub-AS590','sub-AU880','sub-AY910', ...
    'sub-BH088','sub-BI343','sub-BI346','sub-BI360','sub-BI392','sub-BI449','sub-BI460','sub-BI461', ...
    'sub-BI491','sub-BI508','sub-BI513','sub-BI550','sub-BI565','sub-BI632','sub-BI634','sub-BI699', ...
    'sub-BI875','sub-BJ130','sub-BJ131','sub-BJ162','sub-BJ386','sub-BJ469','sub-BJ493','sub-BJ510', ...
    'sub-BJ542','sub-BJ601','sub-BJ643','sub-BJ660','sub-BL307','sub-BL501','sub-BL749','sub-BL779', ...
    'sub-BL884','sub-BM282','sub-BM419','sub-BN573','sub-BN643','sub-BP433','sub-BP446','sub-BP880', ...
    'sub-BP917','sub-BQ069','sub-BQ850','sub-BQ976','sub-BR275','sub-BR276','sub-BR496','sub-BR527', ...
    'sub-BR665','sub-BR952','sub-BR953','sub-XX361'};

%% ------------------ Group assignments (A=Sham, B=Active) ------------------
group_assignments = { ...
    'sub-AO136','B'; 'sub-AO931','B'; 'sub-AP385','B'; 'sub-AQ090','A'; 'sub-AQ822','B'; ...
    'sub-AS590','B'; 'sub-AU880','B'; 'sub-AY910','B'; 'sub-BH088','B'; 'sub-BI343','A'; ...
    'sub-BI346','A'; 'sub-BI360','A'; 'sub-BI392','B'; 'sub-BI449','B'; 'sub-BI460','A'; ...
    'sub-BI461','B'; 'sub-BI491','A'; 'sub-BI508','B'; 'sub-BI513','A'; 'sub-BI550','B'; ...
    'sub-BI565','A'; 'sub-BI632','A'; 'sub-BI634','B'; 'sub-BI699','A'; 'sub-BI875','A'; ...
    'sub-BJ130','B'; 'sub-BJ131','A'; 'sub-BJ162','B'; 'sub-BJ386','A'; 'sub-BJ469','A'; ...
    'sub-BJ493','B'; 'sub-BJ510','A'; 'sub-BJ542','B'; 'sub-BJ601','B'; 'sub-BJ643','A'; ...
    'sub-BJ660','A'; 'sub-BL307','A'; 'sub-BL501','A'; 'sub-BL749','A'; 'sub-BL779','B'; ...
    'sub-BL884','A'; 'sub-BM282','B'; 'sub-BM419','A'; 'sub-BN573','B'; 'sub-BN643','B'; ...
    'sub-BP433','A'; 'sub-BP446','A'; 'sub-BP880','A'; 'sub-BP917','A'; 'sub-BQ069','A'; ...
    'sub-BQ850','B'; 'sub-BQ976','B'; 'sub-BR275','A'; 'sub-BR276','B'; 'sub-BR496','B'; ...
    'sub-BR527','B'; 'sub-BR665','B'; 'sub-BR952','A'; 'sub-BR953','A'; 'sub-XX361','A'};

groupMap = containers.Map(group_assignments(:,1), group_assignments(:,2));

%% ------------------ Read TSV ------------------
opts = detectImportOptions(tsvPath, 'FileType', 'text', 'Delimiter', '\t');
T = readtable(tsvPath, opts);

% Ensure key columns are strings
T.participant_id = string(T.participant_id);
T.session        = string(T.session);

%% ------------------ Filter: participants + session ------------------
targetSession = "before_pre_fmri_arm_1";

keepP = ismember(T.participant_id, string(participant_ids_to_include));
keepS = (T.session == targetSession);
Tf = T(keepP & keepS, :);

% Add Group column
Tf.Group = strings(height(Tf),1);
for i = 1:height(Tf)
    pid = char(Tf.participant_id(i));
    if isKey(groupMap, pid)
        Tf.Group(i) = string(groupMap(pid));
    else
        Tf.Group(i) = "NA";
    end
end
Tf = Tf(Tf.Group ~= "NA", :);

% Medication count numeric (robust)
Tf.medhx_cmeds_count = toDoubleSafe(Tf.medhx_cmeds_count);

A = Tf(Tf.Group=="A", :); % Sham
B = Tf(Tf.Group=="B", :); % Active

fprintf('Session kept: %s\n', targetSession);
fprintf('N Sham (A)   = %d\n', height(A));
fprintf('N Active (B) = %d\n\n', height(B));

%% ============================================================
%  FIGURE 1: Medication count (mean±SEM + scatter) - COLORED
%% ============================================================
figure('Color','w'); hold on;

[mA, seA] = meanSem(A.medhx_cmeds_count);
[mB, seB] = meanSem(B.medhx_cmeds_count);

bar(1, mA, 0.6, 'FaceColor', colA, 'EdgeColor','none');
bar(2, mB, 0.6, 'FaceColor', colB, 'EdgeColor','none');

errorbar(1, mA, seA, 'k', 'LineStyle','none', 'LineWidth', 1.5);
errorbar(2, mB, seB, 'k', 'LineStyle','none', 'LineWidth', 1.5);

% Scatter (jitter)
rng(1);
jA = (rand(height(A),1)-0.5)*0.15;
jB = (rand(height(B),1)-0.5)*0.15;

scatter(1+jA, A.medhx_cmeds_count, 35, colA, 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor','k', 'MarkerEdgeAlpha', 0.20);
scatter(2+jB, B.medhx_cmeds_count, 35, colB, 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor','k', 'MarkerEdgeAlpha', 0.20);

set(gca,'XTick',[1 2],'XTickLabel',{'Sham (A)','Active (B)'});
ylabel('Number of current medications (past 3 months)');
title('Medication count by group (before\_pre\_fmri\_arm\_1)');
box off; grid on;

%% ============================================================
%  FIGURE 2: Distribution (histograms) - COLORED
%% ============================================================
figure('Color','w');

subplot(1,2,1);
histogram(A.medhx_cmeds_count, ...
    'BinMethod','integers', ...
    'FaceColor', colA, ...
    'EdgeColor','none');
title(sprintf('Sham (A), N=%d', height(A)));
xlabel('Medication count'); ylabel('Participants');
box off; grid on;

subplot(1,2,2);
histogram(B.medhx_cmeds_count, ...
    'BinMethod','integers', ...
    'FaceColor', colB, ...
    'EdgeColor','none');
title(sprintf('Active (B), N=%d', height(B)));
xlabel('Medication count'); ylabel('Participants');
box off; grid on;

sgtitle('Medication count distributions (before\_pre\_fmri\_arm\_1)');

%% ============================================================
%  FIGURE 3: Burden categories 0/1/2/3+ (stacked) - GROUP COLORED
%  Idea: each group uses SHADES of its group color (light->dark)
%% ============================================================
bins = [0 1 2 3]; % 3 means "3+"
pA = countBins(A.medhx_cmeds_count, bins);
pB = countBins(B.medhx_cmeds_count, bins);

figure('Color','w');
b = bar([pA; pB], 'stacked', 'BarWidth', 0.7);
set(gca,'XTick',[1 2], 'XTickLabel',{'Sham (A)','Active (B)'});
ylabel('Proportion of participants');
title('Medication burden categories (before\_pre\_fmri\_arm\_1)');
box off; grid on;

% Apply group-colored shades by setting CData per bar object (per segment)
% Each b(k) corresponds to category k (0,1,2,3+)
shades = [0.40 0.65 0.85 1.00]; % light -> full color
for k = 1:numel(b)
    b(k).FaceColor = 'flat';
    b(k).CData = [colA*shades(k); colB*shades(k)];
    b(k).EdgeColor = 'none';
end

legend({'0 meds','1 med','2 meds','3+ meds'}, 'Location','best');

%% ============================================================
%  FIGURE 4: Top medication names (from medhx_cmeds_name_#) - COLORED
%% ============================================================
nameVars = "medhx_cmeds_name_" + string(1:15);

namesA = extractMedNames(A, nameVars);
namesB = extractMedNames(B, nameVars);

namesA = normalizeMedStrings(namesA);
namesB = normalizeMedStrings(namesB);

topN = 15;
[labA, cntA] = topCounts(namesA, topN);
[labB, cntB] = topCounts(namesB, topN);

figure('Color','w');

subplot(1,2,1);
if isempty(cntA)
    text(0.1,0.5,'No medication names recorded','FontSize',11);
    axis off;
else
    barh(cntA, 'FaceColor', colA, 'EdgeColor','none');
    set(gca,'YTick',1:numel(labA),'YTickLabel',labA);
    xlabel('Count');
    title('Top medications – Sham (A)');
    xlim([0 12.1]);
    grid on; box off;
    grid minor
end

subplot(1,2,2);
if isempty(cntB)
    text(0.1,0.5,'No medication names recorded','FontSize',11);
    axis off;
else
    barh(cntB, 'FaceColor', colB, 'EdgeColor','none');
    set(gca,'YTick',1:numel(labB),'YTickLabel',labB);
    xlabel('Count');
    title('Top medications – Active (B)');
    xlim([0 12.1]);
    grid on; box off;
    grid minor
end

sgtitle('Most frequently reported medications (name fields)');

%% ============================================================
%  Summary table (prints in command window)
%% ============================================================
summaryTbl = table( ...
    ["A";"B"], ...
    [height(A); height(B)], ...
    [mA; mB], ...
    [seA; seB], ...
    'VariableNames', {'Group','N','Mean_MedCount','SEM_MedCount'} ...
);
disp(summaryTbl);

%% ===================== Helper functions ======================

function rgb = hex2rgb(hex)
    hex = char(hex);
    if hex(1)=='#', hex = hex(2:end); end
    rgb = [hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))]/255;
end

function x = toDoubleSafe(x)
    % Convert string/cell/char to double safely; keep numeric as-is
    if isnumeric(x)
        x = double(x);
        return;
    end
    if iscell(x)
        x = string(x);
    elseif ischar(x)
        x = string(x);
    end
    if isstring(x)
        x = str2double(x);
    else
        x = double(x);
    end
end

function [m, se] = meanSem(x)
    x = x(:);
    m = mean(x, 'omitnan');
    n = sum(~isnan(x));
    if n <= 1
        se = NaN;
    else
        se = std(x, 'omitnan') / sqrt(n);
    end
end

function props = countBins(x, bins)
    % bins = [0 1 2 3] where 3 means 3+
    x = x(:);
    n = sum(~isnan(x));
    if n==0
        props = zeros(1,numel(bins));
        return;
    end
    c0 = sum(x==0);
    c1 = sum(x==1);
    c2 = sum(x==2);
    c3 = sum(x>=3);
    props = [c0 c1 c2 c3] ./ n;
end

function names = extractMedNames(T, nameVars)
    names = strings(0,1);
    for v = 1:numel(nameVars)
        vn = nameVars(v);
        if ismember(vn, string(T.Properties.VariableNames))
            col = T.(vn);
            names = [names; string(col)];
        end
    end
    names = names(~ismissing(names));
end

function s = normalizeMedStrings(s)
    s = string(s);
    s = lower(strtrim(s));

    % common empties/NA
    bad = (s=="" | s=="na" | s=="n/a" | s=="none" | s=="no" | s=="-" );
    s(bad) = [];

    % optional cleanup of double spaces
    s = regexprep(s, '\s+', ' ');
end

function [labels, counts] = topCounts(names, topN)
    if isempty(names)
        labels = {};
        counts = [];
        return;
    end
    [u,~,idx] = unique(names);
    c = accumarray(idx, 1);
    [cSort, ord] = sort(c, 'descend');
    uSort = u(ord);

    k = min(topN, numel(uSort));
    labels = cellstr(uSort(1:k));
    counts = cSort(1:k);

    % reverse so largest appears at top in barh
    labels = flipud(labels(:));
    counts = flipud(counts(:));
end
