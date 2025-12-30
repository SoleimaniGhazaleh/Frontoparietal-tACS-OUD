function plot_striatum_prepost_bars_scatter_AND_REPORT_STATS()
    clear; clc; close all;

    %% -------- Load data --------
    load('/Volumes/ExtremeSSDD/LIBR_tACS/FurtherAnalysis/Stritum_PrePost.mat'); % Pre, Post (6x60)
    assert(exist('Pre','var')==1 && exist('Post','var')==1, 'Pre/Post not found in .mat file.');
    assert(all(size(Pre)==size(Post)), 'Pre and Post must have the same size.');
    assert(size(Pre,1)==6 && size(Pre,2)==60, 'Expected Pre/Post size = 6x60.');

    %% -------- Colors --------
    colA = [248 118 109]/255; % Sham (Group A)
    colB = [0 191 196]/255;   % Active (Group B)

    %% -------- Group labels (fixed order; length=60) --------
    % A = Sham, B = Active
    grp = { ...
        'B','B','B','A','B','B','B','B','B','A','A','A','B','A','B','A','B','A','B','A', ...
        'A','B','A','A','B','A','B','A','A','B','A','B','B','A','A','A','A','A','B','A', ...
        'B','A','B','B','A','B','A','A','A','A','B','B','A','B','B','B','B','A','A','A' ...
    };
    assert(numel(grp)==60, 'Group label list must be length 60.');

    isSham   = strcmp(grp,'A');
    isActive = strcmp(grp,'B');

    %% -------- Outlier settings --------
    removeOutliers = true;

    % Outlier metric within each group:
    % 'delta' (recommended), 'pre', or 'post'
    outlierMetric = 'delta';
    thrMAD = 3; % median ± thrMAD*MAD

    %% -------- Figure settings --------
   figW = 250; figH = 300;
    barW  = 0.6;
    jitterAmt = 0.10;
    ms = 28;

    %% -------- Stats table container --------
    % Columns: Row, Group, Timepoint, N, Mean, SD
    statsRow = 0;
Stats = table('Size',[0 6], ...
    'VariableTypes',{'double','string','string','double','double','double'}, ...
    'VariableNames',{'RowIdx','Group','Timepoint','N','Mean','SD'});


    %% -------- Loop rows --------
    for r = 1:size(Pre,1)

        preRow  = Pre(r,:);
        postRow = Post(r,:);

        %% Remove outliers (paired) within each group for this row
        if removeOutliers
            [preRow, postRow] = remove_group_outliers_paired(preRow, postRow, isSham,   outlierMetric, thrMAD, r, "Sham(A)");
            [preRow, postRow] = remove_group_outliers_paired(preRow, postRow, isActive, outlierMetric, thrMAD, r, "Active(B)");
        end

        %% Split by group after removal
        preA  = preRow(isSham);    postA = postRow(isSham);
        preB  = preRow(isActive);  postB = postRow(isActive);

        %% Compute mean/SD/N after removal
        % Sham(A)
        [m_preA, sd_preA, n_preA]   = mean_sd_n(preA);
        [m_postA, sd_postA, n_postA]= mean_sd_n(postA);

        % Active(B)
        [m_preB, sd_preB, n_preB]   = mean_sd_n(preB);
        [m_postB, sd_postB, n_postB]= mean_sd_n(postB);

        % Append to stats table
        Stats = [Stats; ...
            make_stat_row(r,"Sham","Pre", n_preA,  m_preA,  sd_preA); ...
            make_stat_row(r,"Sham","Post",n_postA, m_postA, sd_postA); ...
            make_stat_row(r,"Active","Pre", n_preB,  m_preB,  sd_preB); ...
            make_stat_row(r,"Active","Post",n_postB, m_postB, sd_postB)];

        %% Plot bars + scatter (paired)
        mA  = [m_preA  m_postA];
        seA = [sd_preA/sqrt(max(n_preA,1))  sd_postA/sqrt(max(n_postA,1))];

        mB  = [m_preB  m_postB];
        seB = [sd_preB/sqrt(max(n_preB,1))  sd_postB/sqrt(max(n_postB,1))];

        figure('Color','w', ...
            'Name',sprintf('Striatum Row %d', r), ...
            'Position',[100 100 figW figH]);

        x = [1 2];

        % ---- Sham ----
        subplot(1,2,1); hold on;
        bar(x, mA, 'BarWidth', barW, 'EdgeColor','none', 'FaceColor', colA);
        errorbar(x, mA, seA, 'k', 'LineStyle','none', 'LineWidth', 1.2, 'CapSize', 10);

        for i = 1:numel(preA)
            if isnan(preA(i)) || isnan(postA(i)), continue; end
            x1 = 1 + (rand-0.5)*2*jitterAmt;
            x2 = 2 + (rand-0.5)*2*jitterAmt;
            plot([x1 x2], [preA(i) postA(i)], '-', 'Color', [0 0 0 0.25], 'LineWidth', 0.8);
            scatter(x1, preA(i), ms, colA, 'filled', 'MarkerFaceAlpha', 0.75,'MarkerEdgeColor','k');
            scatter(x2, postA(i), ms, colA, 'filled', 'MarkerFaceAlpha', 0.75,'MarkerEdgeColor','k');
        end
        set(gca,'XTick',x,'XTickLabel',{'Pre','Post'});
        title(sprintf('Sham (A) | Row %d', r));
        ylabel('Value');
        xlim([0.5 2.5]); box off;

        % ---- Active ----
        subplot(1,2,2); hold on;
        bar(x, mB, 'BarWidth', barW, 'EdgeColor','none', 'FaceColor', colB);
        errorbar(x, mB, seB, 'k', 'LineStyle','none', 'LineWidth', 1.2, 'CapSize', 10);

        for i = 1:numel(preB)
            if isnan(preB(i)) || isnan(postB(i)), continue; end
            x1 = 1 + (rand-0.5)*2*jitterAmt;
            x2 = 2 + (rand-0.5)*2*jitterAmt;
            plot([x1 x2], [preB(i) postB(i)], '-', 'Color', [0 0 0 0.25], 'LineWidth', 0.8);
            scatter(x1, preB(i), ms, colB, 'filled', 'MarkerFaceAlpha', 0.75,'MarkerEdgeColor','k');
            scatter(x2, postB(i), ms, colB, 'filled', 'MarkerFaceAlpha', 0.75,'MarkerEdgeColor','k');
        end
        set(gca,'XTick',x,'XTickLabel',{'Pre','Post'});
        title(sprintf('Active (B) | Row %d', r));
        ylabel('Value');
        xlim([0.5 2.5]); box off;

        % Match y-lims
        ax1 = subplot(1,2,1);
        ax2 = subplot(1,2,2);
        yAll = [ax1.YLim; ax2.YLim];
        yMin = min(yAll(:)); yMax = max(yAll(:));
        ax1.YLim = [yMin yMax];
        ax2.YLim = [yMin yMax];
    end

    %% -------- Print summary table in Command Window --------
    disp('=== Mean (SD) after outlier removal (per row, group, timepoint) ===');
    disp(Stats);

    %% -------- Save CSV --------
    outCSV = '/Volumes/ExtremeSSDD/LIBR_tACS/FurtherAnalysis/Striatum_PrePost_Stats_AfterOutlierRemoval.csv';
    writetable(Stats, outCSV);
    fprintf('Saved stats table to:\n%s\n', outCSV);
end

%% =========================================================
% Helper: remove paired outliers within one group using MAD
%% =========================================================
function [preRowOut, postRowOut] = remove_group_outliers_paired(preRow, postRow, groupMask, metric, thrMAD, rowIdx, groupLabel)
    preRowOut  = preRow;
    postRowOut = postRow;

    idxG = find(groupMask);
    preG  = preRow(idxG);
    postG = postRow(idxG);

    switch lower(metric)
        case 'delta'
            x = postG - preG;
        case 'pre'
            x = preG;
        case 'post'
            x = postG;
        otherwise
            error('metric must be one of: delta | pre | post');
    end

    medx = median(x,'omitnan');
    madx = mad(x,1);

    if madx <= 0 || isnan(madx)
        return;
    end

    isOutLocal = abs(x - medx) > thrMAD*madx;

    if any(isOutLocal)
        outIdxOriginal = idxG(isOutLocal);
        fprintf('[Row %d] Removing outlier(s) in %s: %s | metric=%s | thr=%g*MAD\n', ...
            rowIdx, groupLabel, mat2str(outIdxOriginal), metric, thrMAD);

        preRowOut(outIdxOriginal)  = NaN;
        postRowOut(outIdxOriginal) = NaN;
    end
end

%% =========================================================
% Helper: mean, SD, N (ignoring NaN)
%% =========================================================
function [m, sd, n] = mean_sd_n(x)
    x = x(~isnan(x));
    n = numel(x);
    if n==0
        m = NaN; sd = NaN;
    elseif n==1
        m = x; sd = NaN;
    else
        m  = mean(x);
        sd = std(x);
    end
end

%% =========================================================
% Helper: one stats row as a table row
%% =========================================================
function T = make_stat_row(rowIdx, groupName, timeName, n, m, sd)
    T = table(rowIdx, string(groupName), string(timeName), n, m, sd, ...
        'VariableNames', {'RowIdx','Group','Timepoint','N','Mean','SD'});
end

