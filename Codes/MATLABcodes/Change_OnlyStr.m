%% ============================================================
%  DELTA (Post - Pre) PLOT: ONE AXIS, TWO BARS (Sham vs Active)
%  - Directed edges (i->j)
%  - Bar (mean±SEM) + scatter
%  - Shared y-axis
%% ============================================================

clear; clc;

%% -------------------------
% USER OPTIONS
%% -------------------------
USE_ABS  = true;     % must match how Zpre/Zpost were computed
PAD_FRAC = 0.10;
FIG_W = 300;
FIG_H = 250;
SAVE_DPI = 300;

%% -------------------------
% Load data (expects Z)
%% -------------------------
p1 = '/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_Striatum/resultsROI_Condition001.mat'; % opioid_pre
p2 = '/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_Striatum/resultsROI_Condition002.mat'; % neutral_pre
p3 = '/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_Striatum/resultsROI_Condition003.mat'; % opioid_post
p4 = '/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_Striatum/resultsROI_Condition004.mat'; % neutral_post

opioid_pre  = load(p1).Z;
neutral_pre = load(p2).Z;
opioid_post = load(p3).Z;
neutral_post= load(p4).Z;

%% -------------------------
% Groups
%% -------------------------
groupA = [0 0 0 1 0 0 0 0 0 1 1 1 0 0 1 0 1 1 0 1 1 0 1 1 0 1 0 1 1 0 1 0 0 0 1 1 1 1 1 0 1 0 1 0 0 1 1 1 1 1 0 0 1 0 0 0 1 1 0 1];
groupB = [1 1 1 0 1 1 1 1 1 0 0 0 1 1 0 1 0 0 1 0 0 1 0 0 1 0 1 0 0 1 0 1 1 1 0 0 0 0 0 1 0 1 0 1 1 0 0 0 0 0 1 1 0 1 1 1 0 0 1 0];

idxA = find(groupA==1);
idxB = find(groupB==1);

%% -------------------------
% ROI labels
%% -------------------------
roi = { ...
'Right Parietal)'
'Right Frontal'
'STR1_r'
'STR1_l'
'STR2_r'
'STR2_l'
'STR3_r'
'STR3_l'
'Left Parietal'
'Left Frontal'};

%% -------------------------
% Colors
%% -------------------------
colA = hex2rgb('#F8766D'); % sham
colB = hex2rgb('#00BFC4'); % active

%% -------------------------
% Compute Δ(Post - Pre)
%% -------------------------
Zpre_raw  = opioid_pre  - neutral_pre;
Zpost_raw = opioid_post - neutral_post;

if USE_ABS
    Zpre  = abs(Zpre_raw);
    Zpost = abs(Zpost_raw);
    ylab = 'Δ|gPPI Z| (Post − Pre)';
else
    Zpre  = Zpre_raw;
    Zpost = Zpost_raw;
    ylab = 'ΔgPPI Z (Post − Pre)';
end

dZ = Zpost - Zpre;  % 14x14x60

%% -------------------------
% Output directory
%% -------------------------
outDir = fullfile(pwd, 'Figures_Delta_TwoBars_Directed_notABS_OnlyStr');
if ~exist(outDir,'dir'); mkdir(outDir); end

%% -------------------------
% Directed edges (i ≠ j)
%% -------------------------
[iList, jList] = find(~eye(10));
nEdges = numel(iList);

fprintf('Generating %d delta figures (two bars)...\n', nEdges);

%% -------------------------
% Loop edges
%% -------------------------
for e = 1:nEdges
    i = iList(e);
    j = jList(e);

    dA = squeeze(dZ(i,j,idxA)); % sham
    dB = squeeze(dZ(i,j,idxB)); % active

    % Shared y-limits
    yAll = [dA(:); dB(:)];
    yAll = yAll(isfinite(yAll));
    r = range(yAll);
    pad = max(PAD_FRAC*r, 0.1);
    yl = [min(yAll)-pad, max(yAll)+pad];

    % Figure
    f = figure('Color','w','Position',[100 100 FIG_W FIG_H]);

    plot_two_bar_scatter(dA, dB, colA, colB);
    ylim(yl);
    ylabel(ylab);
    title(sprintf('%s  →  %s', roi{i}, roi{j}), 'Interpreter','none');
    set(gca,'FontSize',11);

    if ~USE_ABS
        yline(0,'--','Color',[0 0 0 0.35],'LineWidth',1);
    end

    fname = sprintf('Delta_%02d_to_%02d__%s__to__%s.png', ...
        i, j, safe_name(roi{i}), safe_name(roi{j}));
    exportgraphics(f, fullfile(outDir, fname), 'Resolution', SAVE_DPI);

    drawnow;
end

fprintf('DONE.\nSaved to:\n%s\n', outDir);

%% =========================
% Helper functions
%% =========================

function plot_two_bar_scatter(yA, yB, colA, colB)
    yA = yA(isfinite(yA));
    yB = yB(isfinite(yB));

    m  = [mean(yA,'omitnan'), mean(yB,'omitnan')];
    se = [std(yA,'omitnan')/sqrt(numel(yA)), ...
          std(yB,'omitnan')/sqrt(numel(yB))];

    x = [1 2];

    b = bar(x, m, 0.55);
    b.FaceColor = 'flat';
    b.CData = [colA; colB];
    b.EdgeColor = 'none';
    hold on;

    errorbar(x, m, se, 'k.', 'LineWidth', 1.2, 'CapSize', 10);

    jit = 0.10;
    scatter(1 + (rand(numel(yA),1)-0.5)*2*jit, yA, 30, ...
        'MarkerFaceColor', colA, 'MarkerEdgeColor','k','MarkerFaceAlpha',0.85);
    scatter(2 + (rand(numel(yB),1)-0.5)*2*jit, yB, 30, ...
        'MarkerFaceColor', colB, 'MarkerEdgeColor','k','MarkerFaceAlpha',0.85);

    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'Sham','Active'});
    grid on;
    box off;
end

function rgb = hex2rgb(hex)
    hex = strrep(hex,'#','');
    rgb = [hex2dec(hex(1:2)) hex2dec(hex(3:4)) hex2dec(hex(5:6))] ./ 255;
end

function s = safe_name(str)
    s = regexprep(str,'[^\w\-]+','_');
    s = regexprep(s,'_+','_');
    s = regexprep(s,'^_|_$','');
end
