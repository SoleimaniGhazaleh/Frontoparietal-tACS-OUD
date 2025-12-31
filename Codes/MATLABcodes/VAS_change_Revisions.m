%% ============================================================
%  Reviewer Comment #9 — VAS craving patterns across timepoints
%  Goal:
%   (A) Quantify pre-task change (t2-t1) and post-task change (t4-t3)
%       separately for Sham (Group A) and Active (Group B)
%   (B) Test whether Sham shows a reduction during post task (t4-t3)
%       and whether this differs from Active (Group x Block interaction)
%   (C) Compare follow-up day t5 between groups (and also compare change vs baseline)
%
%  Data:
%    /Volumes/ExtremeSSDD/LIBR_tACS/VAS/Craving_Final_comparison.xlsx
%    columns: ID, Group, T0..T5
%
%  Notes:
%   - Uses Welch t-tests for between-group comparisons by default (robust to unequal variances)
%   - Uses sign-rank (nonparametric) as an optional sensitivity test
%   - Uses LME for within-subject Block effect (Pre vs Post task) and Group interaction
%% ============================================================

clear; clc; close all;

xlsxFile = '/Volumes/ExtremeSSDD/LIBR_tACS/VAS/Craving_Final_comparison.xlsx';

% Colors (optional for plots)
colA = [248 118 109]/255; % Sham
colB = [0 191 196]/255;   % Active


outDir = '/Volumes/ExtremeSSDD/LIBR_tACS/VAS/Reviewer9_Results/';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ------------------ Read data ------------------
T = readtable(xlsxFile);

reqVars = {'ID','Group','T0','T1','T2','T3','T4','T5'};
missingVars = setdiff(reqVars, T.Properties.VariableNames);
if ~isempty(missingVars)
    error('Missing columns in Excel: %s', strjoin(missingVars, ', '));
end

T.ID    = string(T.ID);
T.Group = categorical(string(T.Group));

% If your groups are coded differently, set order explicitly:
% T.Group = reordercats(T.Group, {'A','B'}); % A=Sham, B=Active
cats = categories(T.Group);
if numel(cats) ~= 2
    error('Expected exactly 2 groups. Found: %s', strjoin(cats, ', '));
end
grpA = cats{1}; % Sham
grpB = cats{2}; % Active

%% ------------------ Compute task changes ------------------
% Pre-stimulation task change: t2 - t1
% Post-stimulation task change: t4 - t3
T.deltaTask_Pre  = T.T2 - T.T1;
T.deltaTask_Post = T.T4 - T.T3;

% Optional: overall pre-to-post fMRI cue measure (if needed for other comment)
T.deltaCue_Post  = T.T4 - T.T2;

% Follow-up day comparisons
T.deltaFollowup_fromT0 = T.T5 - T.T0;
T.deltaFollowup_fromT1 = T.T5 - T.T1;  % baseline (before pre-task)
T.deltaFollowup_fromT2 = T.T5 - T.T2;  % baseline (after pre-task)

%% ------------------ Helper: descriptive table ------------------
mk_desc = @(x) [mean(x,'omitnan'), std(x,'omitnan'), sum(~isnan(x))];

desc = table();
for g = 1:2
    G = cats{g};
    idx = (T.Group == G);

    [m1,s1,n1] = deal(mk_desc(T.deltaTask_Pre(idx)));
    [m2,s2,n2] = deal(mk_desc(T.deltaTask_Post(idx)));
    [m5,s5,n5] = deal(mk_desc(T.T5(idx)));

    desc = [desc; table( ...
        string(G), n1, m1, s1, ...
        n2, m2, s2, ...
        n5, m5, s5, ...
        'VariableNames', {'Group','N_preTask','Mean_preTask','SD_preTask', ...
                          'N_postTask','Mean_postTask','SD_postTask', ...
                          'N_T5','Mean_T5','SD_T5'})];
end

disp('=== Descriptives by group ===');
disp(desc);

%% ------------------ (A) Within-group: is Sham post-task change < 0 ? ------------------
% One-sample t-test against 0: does deltaTask_Post differ from 0?
xA = T.deltaTask_Post(T.Group==grpA);
xB = T.deltaTask_Post(T.Group==grpB);

xA = xA(~isnan(xA));
xB = xB(~isnan(xB));

fprintf('\n=== Within-group tests: post-task change (T4-T3) vs 0 ===\n');
[~,pA,~,stA] = ttest(xA, 0);  % two-sided
[~,pB,~,stB] = ttest(xB, 0);

fprintf('Sham (%s): mean=%.3f, SD=%.3f, N=%d, t(%d)=%.3f, p=%.4f\n', ...
    string(grpA), mean(xA), std(xA), numel(xA), stA.df, stA.tstat, pA);
fprintf('Active (%s): mean=%.3f, SD=%.3f, N=%d, t(%d)=%.3f, p=%.4f\n', ...
    string(grpB), mean(xB), std(xB), numel(xB), stB.df, stB.tstat, pB);

% Optional nonparametric sensitivity (signrank)
pA_sr = signrank(xA, 0);
pB_sr = signrank(xB, 0);
fprintf('Sign-rank p (Sham) = %.4f | Sign-rank p (Active) = %.4f\n', pA_sr, pB_sr);

%% ------------------ (B) Between-group: does Sham differ from Active on post-task change? ------------------
fprintf('\n=== Between-group test: post-task change (T4-T3) Sham vs Active ===\n');
[~,pWelch,~,stWelch] = ttest2(xA, xB, 'Vartype','unequal'); % Welch
fprintf('Welch t(%0.2f)=%.3f, p=%.4f | mean(Sham)=%.3f, mean(Active)=%.3f\n', ...
    stWelch.df, stWelch.tstat, pWelch, mean(xA), mean(xB));

%% ------------------ (C) LME: deltaTask ~ Group*Block + (1|ID) ------------------
% Build long table with two rows per subject: Pre vs Post task changes
TL = table();
TL.ID = [T.ID; T.ID];
TL.Group = [T.Group; T.Group];
TL.Block = [repmat(categorical("Pre"), height(T), 1); repmat(categorical("Post"), height(T), 1)];
TL.deltaTask = [T.deltaTask_Pre; T.deltaTask_Post];

% Drop missing
TL = TL(~isnan(TL.deltaTask), :);

% Fit LME
fprintf('\n=== LME: deltaTask ~ Group*Block + (1|ID) ===\n');
lme = fitlme(TL, 'deltaTask ~ Group*Block + (1|ID)');
disp(anova(lme));
disp(lme.Coefficients);

%% ------------------ (D) Follow-up day (T5): Sham vs Active ------------------
fprintf('\n=== Follow-up day (T5): between-group comparisons ===\n');

t5A = T.T5(T.Group==grpA); t5A = t5A(~isnan(t5A));
t5B = T.T5(T.Group==grpB); t5B = t5B(~isnan(t5B));

[~,pT5,~,stT5] = ttest2(t5A, t5B, 'Vartype','unequal'); % Welch
fprintf('T5 Welch t(%0.2f)=%.3f, p=%.4f | mean(Sham)=%.3f, mean(Active)=%.3f\n', ...
    stT5.df, stT5.tstat, pT5, mean(t5A), mean(t5B));

% Also compare change-from-baseline to T5 (choose baseline definition)
dA = T.deltaFollowup_fromT1(T.Group==grpA); dA = dA(~isnan(dA));
dB = T.deltaFollowup_fromT1(T.Group==grpB); dB = dB(~isnan(dB));

[~,pD,~,stD] = ttest2(dA, dB, 'Vartype','unequal'); % Welch
fprintf('T5-T1 Welch t(%0.2f)=%.3f, p=%.4f | mean(Sham)=%.3f, mean(Active)=%.3f\n', ...
    stD.df, stD.tstat, pD, mean(dA), mean(dB));

%% ------------------ (E) Optional: quick visualization (bar + scatter) ------------------
% Mean±SEM for deltaTask_Pre and deltaTask_Post by group
figure('Color','w','Position',[100 100 900 360]);

% Prepare
valsA = [T.deltaTask_Pre(T.Group==grpA), T.deltaTask_Post(T.Group==grpA)];
valsB = [T.deltaTask_Pre(T.Group==grpB), T.deltaTask_Post(T.Group==grpB)];

mA = mean(valsA,1,'omitnan');  sA = std(valsA,0,1,'omitnan') ./ sqrt(sum(~isnan(valsA),1));
mB = mean(valsB,1,'omitnan');  sB = std(valsB,0,1,'omitnan') ./ sqrt(sum(~isnan(valsB),1));

x = [1 2];
barWidth = 0.35;

% Bars
b1 = bar(x - barWidth/2, mA, barWidth); hold on;
b2 = bar(x + barWidth/2, mB, barWidth);
b1.FaceColor = colA; b1.EdgeColor = 'k';
b2.FaceColor = colB; b2.EdgeColor = 'k';

% Errors
errorbar(x - barWidth/2, mA, sA, 'k', 'LineStyle','none', 'LineWidth',1);
errorbar(x + barWidth/2, mB, sB, 'k', 'LineStyle','none', 'LineWidth',1);

% Scatter
j = 0.08;
scatter((x(1)-barWidth/2) + (rand(sum(~isnan(valsA(:,1))),1)-0.5)*2*j, valsA(~isnan(valsA(:,1)),1), 30, colA, 'filled', 'MarkerFaceAlpha',0.6,'MarkerEdgeColor','k');
scatter((x(2)-barWidth/2) + (rand(sum(~isnan(valsA(:,2))),1)-0.5)*2*j, valsA(~isnan(valsA(:,2)),2), 30, colA, 'filled', 'MarkerFaceAlpha',0.6,'MarkerEdgeColor','k');

scatter((x(1)+barWidth/2) + (rand(sum(~isnan(valsB(:,1))),1)-0.5)*2*j, valsB(~isnan(valsB(:,1)),1), 30, colB, 'filled', 'MarkerFaceAlpha',0.6,'MarkerEdgeColor','k');
scatter((x(2)+barWidth/2) + (rand(sum(~isnan(valsB(:,2))),1)-0.5)*2*j, valsB(~isnan(valsB(:,2)),2), 30, colB, 'filled', 'MarkerFaceAlpha',0.6,'MarkerEdgeColor','k');

% Cosmetics
xticks(x); xticklabels({'Pre task (T2-T1)','Post task (T4-T3)'});
ylabel('\Delta Craving (VAS)');
yline(0,'k--');
legend({'Sham','Active'}, 'Location','best');
title('Task-induced craving change by group (Pre vs Post stimulation)');
box on; grid on;

%% ------------------ Save outputs (optional) ------------------
writetable(desc, fullfile(outDir, 'Reviewer9_Descriptives_byGroup.csv'));
writetable(TL,   fullfile(outDir, 'Reviewer9_LongTable_forLME.csv'));
fprintf('\nSaved reviewer tables to: %s\n', outDir);
