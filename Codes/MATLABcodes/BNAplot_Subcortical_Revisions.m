%% ==========================================
%  BNA Regions 211–246 (Mean ± SE)
%  Separate Pre/Post per Group
%% ==========================================

clear; clc; close all;

%% ------------------ Paths ------------------
preFile  = '/Users/ghazaleh/Documents/1.PROJECTS/PUBLICATIONS/MolecularPsychiatry/MolecularPsy/Review_round1/BNA_Activation/Final60_Pre.xlsx';
postFile = '/Users/ghazaleh/Documents/1.PROJECTS/PUBLICATIONS/MolecularPsychiatry/MolecularPsy/Review_round1/BNA_Activation/Final60_Post.xlsx';

%% ------------------ Colors ------------------
colB = [248 118 109]/255; % Group A = Sham
colA = [0 191 196]/255;   % Group B = Active

%% ------------------ Plot settings ------------------
capSize   = 4;
xtickStep = 1;        % show every region label (211–246)
figW = 350; figH = 250;
barW = 0.85;

regionIdx = 175:188;  % <<< ONLY THESE REGIONS

%% ------------------ Load & align ------------------
Tpre  = readtable(preFile);
Tpost = readtable(postFile);

Tpre.ID  = string(Tpre.ID);
Tpost.ID = string(Tpost.ID);

T = innerjoin(Tpre, Tpost, 'Keys','ID', ...
    'LeftVariables', Tpre.Properties.VariableNames, ...
    'RightVariables', Tpost.Properties.VariableNames);

if ismember('Group_Tpre', T.Properties.VariableNames)
    Group = string(T.Group_Tpre);
else
    Group = string(T.Group);
end

preVars  = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'Pre_'));
postVars = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'Post_'));

Xpre  = T{:, preVars};
Xpost = T{:, postVars};

%% ------------------ Group indices ------------------
idxA = (Group == "A");
idxB = (Group == "B");

%% ------------------ Helper fns ------------------
nanMean = @(M) mean(M, 1, 'omitnan');
nanSE   = @(M) std(M, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(M),1));

%% =========================
%  GROUP A (regions 211–246)
%% =========================
preA  = Xpre(idxA, regionIdx);
postA = Xpost(idxA, regionIdx);

mA_pre   = nanMean(preA);
mA_post  = nanMean(postA);
seA_pre  = nanSE(preA);
seA_post = nanSE(postA);

% Group A (Sham)
plot_subset(regionIdx, mA_pre, seA_pre, preA, colA, ...
    'Group A (Sham) – Pre (Regions 211–230)', capSize, figW, figH, barW);

plot_subset(regionIdx, mA_post, seA_post, postA, colA, ...
    'Group A (Sham) – Post (Regions 211–230)', capSize, figW, figH, barW);

%% =========================
%  GROUP B (regions 211–246)
%% =========================
preB  = Xpre(idxB, regionIdx);
postB = Xpost(idxB, regionIdx);

mB_pre   = nanMean(preB);
mB_post  = nanMean(postB);
seB_pre  = nanSE(preB);
seB_post = nanSE(postB);

% Group B (Active)
plot_subset(regionIdx, mB_pre, seB_pre, preB, colB, ...
    'Group B (Active) – Pre (Regions 211–230)', capSize, figW, figH, barW);

plot_subset(regionIdx, mB_post, seB_post, postB, colB, ...
    'Group B (Active) – Post (Regions 211–230)', capSize, figW, figH, barW);

%% =========================
%  Δ(Post − Pre), regions 211–246
%% =========================

% ---- Group A ----
dA = abs(postA) - abs(preA);   % [NsubA x regions]
mA_d  = nanMean(dA);
seA_d = nanSE(dA);


% ---- Group B ----
dB = -abs(postB) + abs(preB);   % [NsubB x regions]
mB_d  = nanMean(dB);
seB_d = nanSE(dB);

% Δ(Post − Pre)
plot_subset(regionIdx, mA_d, seA_d, dA, colA, ...
    'Group A (Sham) – Δ(Post − Pre) (Regions 211–230)', capSize, figW, figH, barW);

plot_subset(regionIdx, mB_d, seB_d, dB, colB, ...
    'Group B (Active) – Δ(Post − Pre) (Regions 211–230)', capSize, figW, figH, barW);

%% =========================
%  Local plotting function
%% =========================
function plot_subset(regionIdx, m, se, rawData, col, figTitle, capSize, figW, figH, barW)

    x = regionIdx;                 % true BNA indices
    nSub = size(rawData, 1);

    figure('Color','w','Position',[100 100 figW figH]); hold on;

    % ---- Bars ----
    b = bar(x, m(:), 'BarWidth', barW, 'EdgeColor','none');
    b.FaceColor = col;

    % ---- Error bars ----
    errorbar(x, m(:), se(:), 'k', ...
        'LineStyle','none', 'LineWidth',0.8, 'CapSize',capSize);

    % ---- Scatter (subjects) ----
    jitterAmount = 0.18;  % tweak if you want tighter/wider
    ms = 14;              % marker size

    for r = 1:numel(x)
        y = rawData(:, r);
        ok = ~isnan(y);
        xj = x(r) + (rand(sum(ok),1)-0.5) * jitterAmount;
        scatter(xj, y(ok), ms, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 0.45);
    end

    % ---- Formatting ----
    title(figTitle, 'FontWeight','bold');
    ylabel('BNA Activation');
    xlabel('BNA Region Index');

    set(gca,'FontSize',11);
    xticks(x);
    xlim([min(x)-1 max(x)+1]);
    ylim([-0.45 0.26]);

    grid on; grid minor;
    box off;
end
