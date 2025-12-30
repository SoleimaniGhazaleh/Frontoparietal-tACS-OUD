%% ============================================================
%  FULL PIPELINE: ROI-to-ROI gPPI (14x14x60) - DIRECTED EDGES
%  - Load 4 CONN conditions (Opioid/Neutral x Pre/Post)
%  - Compute contrasts: (Opioid - Neutral) at Pre and at Post
%  - Split by Group A (Sham) vs Group B (Active)
%  - For EACH ROI pair (DIRECTED; i->j and j->i both included), generate:
%       Bar (mean±SEM) + scatter + paired spaghetti (Pre->Post)
%    in two panels: Group A and Group B
%  - SAME y-limits for the two subplots within each figure
%  - Save one PNG per directed ROI-pair
%  - DOES NOT CLOSE FIGURES (per your preference)
%
%  Colors:
%    Group A (Sham)  = #F8766D
%    Group B (Active)= #00BFC4
%% ============================================================

clear; clc;

%% -------------------------
% USER OPTIONS
%% -------------------------
USE_ABS = false;     % keep your current choice (magnitude-only). Set false to keep sign.
PAD_FRAC = 0.10;    % padding fraction for shared y-limits (10%)

%% -------------------------
% Paths to CONN ROI results
%% -------------------------
p1 = '/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_Fox/resultsROI_Condition001.mat'; % opioid_pre
p2 = '/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_Fox/resultsROI_Condition002.mat'; % neutral_pre
p3 = '/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_Fox/resultsROI_Condition003.mat'; % opioid_post
p4 = '/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_Fox/resultsROI_Condition004.mat'; % neutral_post

assert(exist(p1,'file')==2, 'Missing file: %s', p1);
assert(exist(p2,'file')==2, 'Missing file: %s', p2);
assert(exist(p3,'file')==2, 'Missing file: %s', p3);
assert(exist(p4,'file')==2, 'Missing file: %s', p4);

%% -------------------------
% Load (expects field "Z")
%% -------------------------
Cond1 = load(p1); opioid_pre  = Cond1.Z;
Cond2 = load(p2); neutral_pre = Cond2.Z;
Cond3 = load(p3); opioid_post = Cond3.Z;
Cond4 = load(p4); neutral_post= Cond4.Z;

%% -------------------------
% Sanity checks
%% -------------------------
assert(isequal(size(opioid_pre),  [4 4 60]), 'opioid_pre is not 14x14x60');
assert(isequal(size(neutral_pre), [4 4 60]), 'neutral_pre is not 14x14x60');
assert(isequal(size(opioid_post), [4 4 60]), 'opioid_post is not 14x14x60');
assert(isequal(size(neutral_post),[4 4 60]), 'neutral_post is not 14x14x60');

%% -------------------------
% Groups (A=sham, B=active)
%% -------------------------
groupA = [0 0 0 1 0 0 0 0 0 1 1 1 0 0 1 0 1 1 0 1 1 0 1 1 0 1 0 1 1 0 1 0 0 0 1 1 1 1 1 0 1 0 1 0 0 1 1 1 1 1 0 0 1 0 0 0 1 1 0 1];
groupB = [1 1 1 0 1 1 1 1 1 0 0 0 1 1 0 1 0 0 1 0 0 1 0 0 1 0 1 0 0 1 0 1 1 1 0 0 0 0 0 1 0 1 0 1 1 0 0 0 0 0 1 1 0 1 1 1 0 0 1 0];

assert(numel(groupA)==60 && numel(groupB)==60, 'Group vectors must be length 60');
assert(all((groupA+groupB)==1), 'Groups must be complementary (A+B==1 for all subjects)');

idxA = find(groupA==1);
idxB = find(groupB==1);

%% -------------------------
% ROI labels (order must match 14x14)
%% -------------------------
roi = { ...
'R-Parietal'
'R-Frontal'
'L-Parietal'
'L-Frontal'};

assert(numel(roi)==4, 'ROI label list must have 4 items');

%% -------------------------
% Colors (requested)
%% -------------------------
colA = hex2rgb('#F8766D'); % sham
colB = hex2rgb('#00BFC4'); % active

%% ------------------------------------------------------------
% Compute contrasts: (Opioid - Neutral) for Pre and for Post
%% ------------------------------------------------------------
Zpre_raw  = opioid_pre  - neutral_pre;   % 14x14x60
Zpost_raw = opioid_post - neutral_post;  % 14x14x60

if USE_ABS
    Zpre  = abs(Zpre_raw);
    Zpost = abs(Zpost_raw);
    ylab = '|gPPI Z| (|Opioid − Neutral|)';
else
    Zpre  = Zpre_raw;
    Zpost = Zpost_raw;
    ylab = 'gPPI Z (Opioid − Neutral)';
end

%% ------------------------------------------------------------
% Plot settings + output directory
%% ------------------------------------------------------------
tag = ternary(USE_ABS, 'ABS', 'SIGNED');
outDir = fullfile(pwd, ['Figures_OpioidMinusNeutral_PrePost_' tag '_DIRECTED_SHAREDYLIM_OnlyFrontoParietal']);
if ~exist(outDir, 'dir'); mkdir(outDir); end

% ALL directed edges: i ~= j (includes i->j and j->i)
[iList, jList] = find(~eye(4));
nEdges = numel(iList);

fprintf('Generating %d directed edge figures (i->j, i~=j)...\n', nEdges);

%% ------------------------------------------------------------
% Loop over ALL directed ROI pairs and save figures
%% ------------------------------------------------------------
for e = 1:nEdges
    i = iList(e);
    j = jList(e);

    % Group A (Sham)
    preA  = squeeze(Zpre(i,j,idxA));
    postA = squeeze(Zpost(i,j,idxA));

    % Group B (Active)
    preB  = squeeze(Zpre(i,j,idxB));
    postB = squeeze(Zpost(i,j,idxB));

    % Compute shared y-limits for THIS figure (both panels together)
    yAll = [preA(:); postA(:); preB(:); postB(:)];
    yAll = yAll(isfinite(yAll));

    if isempty(yAll)
        yl = [-1 1]; % fallback
    else
        r = max(yAll) - min(yAll);
        pad = PAD_FRAC * r;
        if pad == 0
            pad = 0.1; % fallback padding when all values identical
        end
        yl = [min(yAll)-pad, max(yAll)+pad];
    end

    % Make figure (do NOT close)
    f = figure('Color','w','Position',[100 100 600 250]);

    % Panel A
    ax1 = subplot(1,2,1);
    plot_bar_scatter_spaghetti(preA, postA, colA, USE_ABS);
    title(sprintf('%s  →  %s', roi{i}, roi{j}), 'Interpreter','none');
    ylabel(ylab);
    ylim(yl);
    set(gca,'FontSize',11);

    % Panel B
    ax2 = subplot(1,2,2);
    plot_bar_scatter_spaghetti(preB, postB, colB, USE_ABS);
    title(sprintf('%s  →  %s', roi{i}, roi{j}), 'Interpreter','none');
    ylabel(ylab);
    ylim(yl);
    set(gca,'FontSize',11);

    % Optional: also link y-axes (keeps them synced if edited)
    linkaxes([ax1 ax2], 'y');

    sgtitle(sprintf('Opioid − Neutral (Pre vs Post) | Directed Edge %d/%d', e, nEdges));

    % Save
    fname = sprintf('Edge_%02d_to_%02d__%s__to__%s.png', i, j, safe_name(roi{i}), safe_name(roi{j}));
    exportgraphics(f, fullfile(outDir, fname), 'Resolution', 200);
    

    drawnow;

    if mod(e,10)==0
        fprintf('  Saved %d/%d\n', e, nEdges);
    end
end

fprintf('DONE.\nSaved figures to:\n%s\n', outDir);

%% =========================
% Helper functions (local)
%% =========================

function plot_bar_scatter_spaghetti(pre, post, col, use_abs)
    % Bar (mean±SEM) + scatter + paired spaghetti

    ok = isfinite(pre) & isfinite(post);
    pre  = pre(ok);
    post = post(ok);

    n = numel(pre);
    if n < 2
        cla; text(0.5,0.5,'Not enough data','HorizontalAlignment','center');
        return;
    end

    m   = [mean(pre,'omitnan'),  mean(post,'omitnan')];
    sem = [std(pre,'omitnan')/sqrt(n), std(post,'omitnan')/sqrt(n)];

    x = [1 2];

    b = bar(x, m, 0.55);
    b.FaceColor = 'flat';
    b.CData = [col; col];
    b.EdgeColor = 'none';
    hold on;

    errorbar(x, m, sem, 'k.', 'LineWidth', 1.2, 'CapSize', 10);

    jitter = 0.08;
    x1 = 1 + (rand(n,1)-0.5)*2*jitter;
    x2 = 2 + (rand(n,1)-0.5)*2*jitter;

    for s = 1:n
        plot([x1(s) x2(s)], [pre(s) post(s)], '-', 'Color', [0 0 0 0.18], 'LineWidth', 0.8);
    end

    scatter(x1, pre,  28, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.85);
    scatter(x2, post, 28, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.85);

    % Zero line only makes sense for signed data
    if ~use_abs
        yline(0, '--', 'Color', [0 0 0 0.35], 'LineWidth', 1);
    end

    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'Pre','Post'});
    grid on;
    box off;
end

function rgb = hex2rgb(hex)
    hex = strrep(hex,'#','');
    rgb = [hex2dec(hex(1:2)) hex2dec(hex(3:4)) hex2dec(hex(5:6))] ./ 255;
end

function s = safe_name(str)
    s = regexprep(str, '[^\w\-]+', '_');
    s = regexprep(s, '_+', '_');
    s = regexprep(s, '^_|_$', '');
    if numel(s) > 80
        s = s(1:80);
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
