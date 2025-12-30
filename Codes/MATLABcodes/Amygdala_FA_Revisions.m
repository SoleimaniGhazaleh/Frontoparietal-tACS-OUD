%% =======================
%  BNA 246-Region Bar Plots (Mean ± SE)
%  - Align Pre/Post by ID (robust)
%  - Plot 246 bars in ONE figure for:
%       Group A - Pre
%       Group A - Post
%       Group B - Pre
%       Group B - Post
%    (i.e., Pre and Post are plotted SEPARATELY per group)
%  - NaN-safe (per region), SE uses N_nonNaN per region
%  - Small errorbar caps + thicker bars
%% =======================

clear; clc;
close all

%% ------------------ Paths ------------------
preFile  = '/Users/ghazaleh/Documents/1.PROJECTS/PUBLICATIONS/MolecularPsychiatry/MolecularPsy/Review_round1/BNA_Activation/Final60_Pre.xlsx';
postFile = '/Users/ghazaleh/Documents/1.PROJECTS/PUBLICATIONS/MolecularPsychiatry/MolecularPsy/Review_round1/BNA_Activation/Final60_Post.xlsx';

%% ------------------ Colors ------------------
colA = [248 118 109]/255; % #F8766D  Group A = Sham
colB = [0 191 196]/255;   % #00BFC4  Group B = Active

%% ------------------ Plot settings ------------------
capSize   = 2;        % smaller caps (try 2-6)
xtickStep = 10;       % show every 10th region tick
figW = 250; figH = 300;
barW = 0.90;          % bar thickness (must be <= 1)

%% ------------------ Load data ------------------
Tpre  = readtable(preFile);
Tpost = readtable(postFile);

% Ensure comparable ID types
Tpre.ID  = string(Tpre.ID);
Tpost.ID = string(Tpost.ID);

% Align rows by ID (inner join keeps only overlapping IDs)
T = innerjoin(Tpre, Tpost, 'Keys','ID', ...
    'LeftVariables', Tpre.Properties.VariableNames, ...
    'RightVariables', Tpost.Properties.VariableNames);

% Get Group labels (prefer Pre file if duplicated by join)
if ismember('Group_Tpre', T.Properties.VariableNames)
    Group = string(T.Group_Tpre);
elseif ismember('Group', T.Properties.VariableNames)
    Group = string(T.Group);
else
    error('Group column not found after join. Check table variable names.');
end

% Extract BNA matrices
preVars  = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'Pre_'));
postVars = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'Post_'));

Xpre  = T{:, preVars};    % [Nsub x 246]
Xpost = T{:, postVars};   % [Nsub x 246]

% Sanity checks
if ~isequal(size(Xpre), size(Xpost))
    error('Pre/Post matrices are different sizes after join: Pre=%s, Post=%s', ...
        mat2str(size(Xpre)), mat2str(size(Xpost)));
end
if size(Xpre,2) ~= 246
    warning('Expected 246 regions. Found %d columns (Pre/Post).', size(Xpre,2));
end

fprintf('Aligned subjects (intersection of IDs): %d\n', size(Xpre,1));

%% ------------------ Group indices ------------------
idxA = (Group == "A");
idxB = (Group == "B");

fprintf('Group A (Sham): %d subjects\n', sum(idxA));
fprintf('Group B (Active): %d subjects\n', sum(idxB));

%% ------------------ Helper fns: NaN-safe mean & SE per region ------------------
nanMean = @(M) mean(M, 1, 'omitnan'); % 1 x 246
nanSE   = @(M) std(M, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(M),1)); % 1 x 246

%% =========================
%  GROUP A: Pre and Post (SEPARATE FIGURES)
%% =========================
preA  = Xpre(idxA,:);
postA = Xpost(idxA,:);

%% =========================
%  GROUP B: Pre and Post (SEPARATE FIGURES)
%% =========================
preB  = Xpre(idxB,:);
postB = Xpost(idxB,:);
%% ==========================================
%  Paired bar + scatter for selected BNA regions (211–214)
%  Separate panels for Group A (Sham) and Group B (Active)
%% ==========================================

regions = [211 212 213 214];

% Figure size (change if you want)
figW2 = 250;
figH2 = 300;

% Scatter settings
ms = 28;              % marker size
jitterAmt = 0.08;     % jitter around x=1 and x=2
lineAlpha = 0.25;

for rr = 1:numel(regions)
    reg = regions(rr);

    % Extract subject-level paired values for each group (Nsub x 1)
    yPreA  = preA(:,reg);
    yPostA = postA(:,reg);

    yPreB  = preB(:,reg);
    yPostB = postB(:,reg);

    % Means + SEM (NaN-safe)
    mA  = [mean(yPreA,'omitnan')  mean(yPostA,'omitnan')];
    seA = [std(yPreA,'omitnan')/sqrt(sum(~isnan(yPreA)))  std(yPostA,'omitnan')/sqrt(sum(~isnan(yPostA)))];

    mB  = [mean(yPreB,'omitnan')  mean(yPostB,'omitnan')];
    seB = [std(yPreB,'omitnan')/sqrt(sum(~isnan(yPreB)))  std(yPostB,'omitnan')/sqrt(sum(~isnan(yPostB)))];

    % Create figure: 1 row, 2 columns (Sham vs Active)
    figure('Color','w', 'Position',[200 200 figW2 figH2], ...
        'Name', sprintf('BNA Region %d (Pre/Post)', reg));

    x = [1 2];

    %% ------------------ Group A (Sham) ------------------
    subplot(1,2,1); hold on;

    bar(x, mA, 'BarWidth', 0.65, 'EdgeColor','none', 'FaceColor', colA); % colB = Sham color
    errorbar(x, mA, seA, 'k', 'LineStyle','none', 'LineWidth', 1.2, 'CapSize', 10);

    % paired lines + scatter (only where both are non-NaN)
    okA = ~isnan(yPreA) & ~isnan(yPostA);
    yPreA2 = yPreA(okA);
    yPostA2 = yPostA(okA);

    for i = 1:numel(yPreA2)
        x1 = 1 + (rand-0.5)*2*jitterAmt;
        x2 = 2 + (rand-0.5)*2*jitterAmt;
        plot([x1 x2], [yPreA2(i) yPostA2(i)], '-', 'Color', [0 0 0 lineAlpha], 'LineWidth', 0.8);
        scatter(x1, yPreA2(i), ms, colA, 'filled', 'MarkerFaceAlpha', 0.75,'MarkerEdgeColor','k');
        scatter(x2, yPostA2(i), ms, colA, 'filled', 'MarkerFaceAlpha', 0.75,'MarkerEdgeColor','k');
    end

    set(gca,'XTick',x,'XTickLabel',{'Pre','Post'});
    title(sprintf('Sham (A) | Region %d', reg), 'Interpreter','none');
    ylabel('Value');
    xlim([0.5 2.5]);
    box off;

    %% ------------------ Group B (Active) ------------------
    subplot(1,2,2); hold on;

    bar(x, mB, 'BarWidth', 0.65, 'EdgeColor','none', 'FaceColor', colB); % colA = Active color
    errorbar(x, mB, seB, 'k', 'LineStyle','none', 'LineWidth', 1.2, 'CapSize', 10);

    okB = ~isnan(yPreB) & ~isnan(yPostB);
    yPreB2 = yPreB(okB);
    yPostB2 = yPostB(okB);

    for i = 1:numel(yPreB2)
        x1 = 1 + (rand-0.5)*2*jitterAmt;
        x2 = 2 + (rand-0.5)*2*jitterAmt;
        plot([x1 x2], [yPreB2(i) yPostB2(i)], '-', 'Color', [0 0 0 lineAlpha], 'LineWidth', 0.8);
        scatter(x1, yPreB2(i), ms, colB, 'filled', 'MarkerFaceAlpha', 0.75,'MarkerEdgeColor','k');
        scatter(x2, yPostB2(i), ms, colB, 'filled', 'MarkerFaceAlpha', 0.75,'MarkerEdgeColor','k');
    end

    set(gca,'XTick',x,'XTickLabel',{'Pre','Post'});
    title(sprintf('Active (B) | Region %d', reg), 'Interpreter','none');
    ylabel('Value');
    xlim([0.5 2.5]);
    box off;

    % ---- Match y-limits across the two panels for this region ----
    ax1 = subplot(1,2,1);
    ax2 = subplot(1,2,2);
    yAll = [ax1.YLim; ax2.YLim];
    yMin = min(yAll(:)); yMax = max(yAll(:));
    ax1.YLim = [yMin yMax];
    ax2.YLim = [yMin yMax];
end
