%% ============================================================
%  BCS (Brief Drug Craving Scale) - ONE FIGURE
%  Rows = Questions (bmcs_1..bmcs_4)
%  Cols = Timepoints (baseline, before_pre_fmri, day_after)
%  Each subplot: Bar (mean±SEM) + scatter, Sham vs Active
%  Y-lims: bmcs_1..3 -> [0 4], bmcs_4 -> [0 30]
%% ============================================================

clear; clc; close all;

%% ------------------ Path ------------------
tsvPath = '/Volumes/ExtremeSSDD/LIBR_tACS/Demographic/phenotype/brief_drug_craving_scale_bcs.tsv';

%% ------------------ Colors ------------------
colA = hex2rgb('#F8766D'); % Sham
colB = hex2rgb('#00BFC4'); % Active

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

%% ------------------ Selected timepoints ------------------
timepoints = { ...
    'baseline_arm_1'
    'before_pre_fmri_arm_1'
    'day_after_arm_1'
    };
tpLabels = {'Baseline','Before pre-fMRI','Day after'};

%% ------------------ BCS questions ------------------
bcs_vars  = {'bmcs_1','bmcs_2','bmcs_3','bmcs_4'};
bcs_names = { ...
    'BCS-1: Intensity (0–4)'
    'BCS-2: Frequency (0–4)'
    'BCS-3: Duration (0–4)'
    'BCS-4: Count (0–30)'
    };

%% ------------------ Read TSV ------------------
opts = detectImportOptions(tsvPath,'FileType','text','Delimiter','\t');
T = readtable(tsvPath,opts);
T.participant_id = string(T.participant_id);
T.session        = string(T.session);

% Filter participants
T = T(ismember(T.participant_id, string(participant_ids_to_include)), :);

% Add group
GA = string(group_assignments(:,1));
GB = string(group_assignments(:,2));
T.Group = strings(height(T),1);
for i = 1:height(T)
    idx = find(GA == T.participant_id(i),1);
    if ~isempty(idx), T.Group(i) = GB(idx); else, T.Group(i) = "NA"; end
end
T = T(T.Group~="NA",:);

%% ================== PLOT: ONE FIGURE (4x3) ==================
fig = figure('Color','w','Position',[100 100 800 600]);
%sgtitle('Brief Drug Craving Scale (BCS): Sham vs Active across timepoints','FontWeight','bold');

for r = 1:numel(bcs_vars)
    for c = 1:numel(timepoints)

        ax = subplot(numel(bcs_vars), numel(timepoints), (r-1)*numel(timepoints)+c);
        hold(ax,'on');

        tp = timepoints{c};
        idxT = T.session == tp;

        xA = toNumeric(T{idxT & T.Group=="A", bcs_vars{r}});
        xB = toNumeric(T{idxT & T.Group=="B", bcs_vars{r}});

        xA = xA(isfinite(xA));
        xB = xB(isfinite(xB));

        means = [mean(xA,'omitnan'), mean(xB,'omitnan')];
        sems  = [std(xA,'omitnan')/sqrt(max(numel(xA),1)), ...
                 std(xB,'omitnan')/sqrt(max(numel(xB),1))];

        % Bars
        b = bar([1 2], means, 'FaceColor','flat','EdgeColor','k');
        b.CData(1,:) = colA;
        b.CData(2,:) = colB;

        % Error bars
        errorbar([1 2], means, sems, 'k','LineStyle','none','LineWidth',1);

        % Scatter
        jitter = 0.08;
        scatter(1 + randn(size(xA))*jitter, xA, 18, colA, 'filled','MarkerFaceAlpha',0.6,'MarkerEdgeColor','k');
        scatter(2 + randn(size(xB))*jitter, xB, 18, colB, 'filled','MarkerFaceAlpha',0.6,'MarkerEdgeColor','k');

        set(ax,'XTick',[1 2],'XTickLabel',{'Sham','Active'});
        xtickangle(ax,0);
        grid(ax,'on'); box(ax,'off');

        % Row titles on leftmost column, col titles on top row
        if c == 1
            ylabel(ax, bcs_names{r}, 'FontSize',10, 'FontWeight','bold');
        end
        if r == 1
           % title(ax, tpLabels{c}, 'FontSize',11);
        end

        % Fixed y-limits as requested
        if r <= 3
            ylim(ax, [0 4]);
            yticks(ax, 0:1:4);
        else
            ylim(ax, [0 30]);
            yticks(ax, 0:5:30);
        end
    end
end

% Optional legend (single, outside)
%lg = legend({'Mean','Mean','SEM','Sham indiv','Active indiv'}, 'Location','southoutside');
%lg.Box = 'off';

%% ------------------ Save figure ------------------
outDir = fullfile(pwd,'figures');
if ~exist(outDir,'dir'); mkdir(outDir); end

figFile = fullfile(outDir,'BCS_Sham_vs_Active_3Timepoints');

exportgraphics(fig, [figFile '.png'], 'Resolution',300);
exportgraphics(fig, [figFile '.pdf'], 'ContentType','vector');

%%
%% ============================================================
%  REPORT: Descriptives + Group Comparisons (BCS)
%% ============================================================

reportRows = table();

for r = 1:numel(bcs_vars)
    for c = 1:numel(timepoints)

        tp = timepoints{c};
        idxT = T.session == tp;

        xA = toNumeric(T{idxT & T.Group=="A", bcs_vars{r}});
        xB = toNumeric(T{idxT & T.Group=="B", bcs_vars{r}});

        xA = xA(isfinite(xA));
        xB = xB(isfinite(xB));

        % Descriptives
        nA = numel(xA);   nB = numel(xB);
        mA = mean(xA,'omitnan');  mB = mean(xB,'omitnan');
        sA = std(xA,'omitnan');   sB = std(xB,'omitnan');
        semA = sA / sqrt(max(nA,1));
        semB = sB / sqrt(max(nB,1));

        % Stats (Welch t-test)
        p_t = NaN; tstat = NaN; df = NaN;
        if nA >= 3 && nB >= 3
            [~, p_t, ~, stats] = ttest2(xA, xB, 'Vartype','unequal');
            tstat = stats.tstat;
            df    = stats.df;
        end

        % Nonparametric check
        p_rs = NaN;
        if nA >= 3 && nB >= 3
            p_rs = ranksum(xA, xB);
        end

        % Effect size: Hedges g (Active − Sham)
        g = NaN;
        if nA >= 2 && nB >= 2
            sp = sqrt(((nA-1)*var(xA,1) + (nB-1)*var(xB,1)) / max(nA+nB-2,1));
            if sp > 0
                d = (mB - mA) / sp;
                J = 1 - (3/(4*(nA+nB)-9));
                g = J*d;
            end
        end

        % Store row
        newRow = table( ...
            string(bcs_vars{r}), string(bcs_names{r}), ...
            string(tpLabels{c}), ...
            nA, mA, semA, ...
            nB, mB, semB, ...
            tstat, df, p_t, p_rs, g, ...
            'VariableNames', { ...
                'BCS_Item','BCS_Label','Timepoint', ...
                'N_Sham','Mean_Sham','SEM_Sham', ...
                'N_Active','Mean_Active','SEM_Active', ...
                't','df','p_ttestWelch','p_ranksum','Hedges_g_ActiveMinusSham'} ...
            );

        reportRows = [reportRows; newRow]; %#ok<AGROW>
    end
end

% Round numeric values for readability
numVars = varfun(@isnumeric, reportRows, 'OutputFormat','uniform');
for v = find(numVars)
    reportRows.(reportRows.Properties.VariableNames{v}) = ...
        round(reportRows.(reportRows.Properties.VariableNames{v}), 4);
end

disp(' ');
disp('=== BCS Results: Sham vs Active ===');
disp(reportRows);

% Save CSV
outCSV = fullfile(outDir,'BCS_GroupComparison_Report.csv');
writetable(reportRows, outCSV);

fprintf('\nBCS results table saved to:\n%s\n', outCSV);

%% ------------------ Helpers ------------------
function x = toNumeric(x)
    if istable(x), x = table2array(x); end
    if iscell(x),  x = string(x); end
    if ischar(x),  x = string(x); end
    if isstring(x)
        x = str2double(strrep(x,"NA","NaN"));
    end
end

function rgb = hex2rgb(hexStr)
    if startsWith(hexStr,'#'), hexStr = hexStr(2:end); end
    rgb = reshape(sscanf(hexStr,'%2x'),1,3)/255;
end
