%% ============================================================
%  CEE DRUGS: After-post fMRI (arm_1) group comparison + report
%  - Reads: /Volumes/ExtremeSSDD/LIBR_tACS/Demographic/phenotype/cee_drugs.tsv
%  - Keeps only: session == "after_post_fmri_arm_1"
%  - Keeps only your participant list
%  - Adds Group (A=Sham, B=Active) from your mapping
%  - Robustly converts numeric columns (handles cell/string/"NA")
%  - Reports descriptives + group tests for selected subscales
%  - Saves CSV report
%% ============================================================

clear; clc;

%% ------------------ Paths ------------------
tsvPath = '/Volumes/ExtremeSSDD/LIBR_tACS/Demographic/phenotype/cee_drugs.tsv';
outDir  = fullfile(pwd, 'CEE_reports');
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% ------------------ Colors (optional) ------------------
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

%% ------------------ Read TSV (robust) ------------------
opts = detectImportOptions(tsvPath, 'FileType','text', 'Delimiter','\t');
% Force text to remain text; we will convert as needed
opts = setvaropts(opts, opts.VariableNames, 'WhitespaceRule','preserve');
opts = setvaropts(opts, opts.VariableNames, 'EmptyFieldRule','auto');
T = readtable(tsvPath, opts);

% Ensure required columns exist
reqCols = {'participant_id','session'};
for i = 1:numel(reqCols)
    if ~ismember(reqCols{i}, T.Properties.VariableNames)
        error('Missing required column "%s" in %s', reqCols{i}, tsvPath);
    end
end

% Normalize ID/session to string
T.participant_id = string(T.participant_id);
T.session        = string(T.session);

%% ------------------ Filter: participants + session ------------------
keepIDs = string(participant_ids_to_include(:));
T = T(ismember(T.participant_id, keepIDs), :);

targetSession = "after_post_fmri_arm_1";
T = T(T.session == targetSession, :);

fprintf('\nKept N=%d rows after filtering participants and session="%s"\n', height(T), targetSession);

%% ------------------ Add Group column ------------------
GA = string(group_assignments(:,1));
GB = string(group_assignments(:,2));

grp = strings(height(T),1);
for i = 1:height(T)
    idx = find(GA == T.participant_id(i), 1);
    if ~isempty(idx)
        grp(i) = GB(idx);   % "A" or "B"
    else
        grp(i) = "NA";
    end
end
T.Group = grp;

% Remove missing group rows if any
T = T(T.Group ~= "NA", :);
T.GroupLabel = categorical(T.Group, ["A","B"], ["Sham","Active"]);

%% ------------------ Subscales to report ------------------
subscales = { ...
    'success_control_score'        % Executive control / inhibition
    'unsuccess_control_score'      % Disinhibition / loss of control
    'reappraisal_score'            % Cognitive control / reappraisal
    'attentional_score'            % Attentional capture
    'ruminative_score'             % Perseveration/rumination
    'pos_appetitive_score'
    'neg_emotional_score'
    'interoceptive_score'
    'metacog_score'
    'verbalizing_score'
    'pos_memory_score'
    'neg_memory_score'
    'gen_memory_score'
    };

% Keep only columns that exist
subscales = subscales(ismember(subscales, T.Properties.VariableNames));
if isempty(subscales)
    error('None of the requested subscale columns were found in the TSV.');
end

%% ------------------ Descriptives + Tests ------------------
descRows = table();
testRows = table();

for s = 1:numel(subscales)
    vname = subscales{s};

    xA = toNumeric(T{T.Group=="A", vname});
    xB = toNumeric(T{T.Group=="B", vname});

    xA = xA(isfinite(xA));
    xB = xB(isfinite(xB));

    % Descriptives
    rowD = table( ...
        string(vname), ...
        numel(xA), mean(xA,'omitnan'), std(xA,'omitnan'), std(xA,'omitnan')/sqrt(max(numel(xA),1)), ...
        numel(xB), mean(xB,'omitnan'), std(xB,'omitnan'), std(xB,'omitnan')/sqrt(max(numel(xB),1)), ...
        'VariableNames', {'Measure','N_Sham','Mean_Sham','SD_Sham','SEM_Sham','N_Active','Mean_Active','SD_Active','SEM_Active'} ...
        );
    descRows = [descRows; rowD]; %#ok<AGROW>

    % Welch t-test
    p_t = NaN; tstat = NaN; df = NaN;
    if numel(xA) >= 3 && numel(xB) >= 3
        [~, p_t, ~, stats] = ttest2(xA, xB, 'Vartype','unequal');
        tstat = stats.tstat;
        df    = stats.df;
    end

    % Rank-sum
    p_w = NaN;
    if numel(xA) >= 3 && numel(xB) >= 3
        p_w = ranksum(xA, xB);
    end

    % Hedges' g (Active - Sham)
    g = NaN;
    if numel(xA) >= 2 && numel(xB) >= 2
        m1 = mean(xA,'omitnan'); m2 = mean(xB,'omitnan');
        v1 = var(xA,1,'omitnan'); v2 = var(xB,1,'omitnan');
        n1 = numel(xA); n2 = numel(xB);
        sp = sqrt(((n1-1)*v1 + (n2-1)*v2) / max(n1+n2-2,1));
        if sp > 0
            d  = (m2 - m1) / sp;
            J  = 1 - (3/(4*(n1+n2)-9));
            g  = J*d;
        end
    end

    rowT = table(string(vname), p_t, tstat, df, p_w, g, ...
        'VariableNames', {'Measure','p_ttestWelch','t','df','p_ranksum','Hedges_g_ActiveMinusSham'});
    testRows = [testRows; rowT]; %#ok<AGROW>
end

%% ------------------ Merge and format ------------------
Report = outerjoin(descRows, testRows, 'Keys','Measure', 'MergeKeys',true);

% Keep a nice order (executive-control relevant first)
priority = ["success_control_score","unsuccess_control_score","reappraisal_score","attentional_score","ruminative_score"];
[~, idxP] = ismember(priority, string(Report.Measure));
idxP = idxP(idxP>0);
idxRest = setdiff(1:height(Report), idxP, 'stable');
Report = Report([idxP(:); idxRest(:)], :);

% Optional rounding for readability
numVars = varfun(@isnumeric, Report, 'OutputFormat','uniform');
for j = find(numVars)
    Report.(Report.Properties.VariableNames{j}) = round(Report.(Report.Properties.VariableNames{j}), 4);
end

disp(' ');
disp('=== CEE After-post fMRI (arm_1): Group Report (Sham vs Active) ===');
disp(Report);

%% ------------------ Save ------------------
outCSV = fullfile(outDir, sprintf('CEE_after_post_fmri_arm1_GroupReport_%s.csv', datestr(now,'yyyymmdd_HHMMSS')));
writetable(Report, outCSV);
fprintf('\nSaved report to:\n%s\n', outCSV);

%% ========================= Helpers =========================

function rgb = hex2rgb(hexStr)
%HEX2RGB Convert hex color (e.g. '#F8766D') to [r g b] in 0-1
    hexStr = char(hexStr);
    if startsWith(hexStr,'#'); hexStr = hexStr(2:end); end
    if numel(hexStr) ~= 6; error('hex2rgb expects 6 hex digits.'); end
    r = hex2dec(hexStr(1:2));
    g = hex2dec(hexStr(3:4));
    b = hex2dec(hexStr(5:6));
    rgb = [r g b]/255;
end

function x = toNumeric(x)
%TONUMERIC Safely convert table column to double
% - Handles cell arrays, string arrays, char arrays
% - Converts "NA", "", etc. to NaN via str2double
    if istable(x)
        x = table2array(x);
    end
    if iscell(x)
        x = string(x);
    end
    if ischar(x)
        x = string(x);
    end
    if isstring(x)
        x = str2double(strrep(x, "NA", "NaN"));
    end
    if ~isnumeric(x)
        % last fallback
        try
            x = double(x);
        catch
            x = NaN(size(x));
        end
    end
end
