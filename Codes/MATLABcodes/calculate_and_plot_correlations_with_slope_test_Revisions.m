function calculate_and_plot_correlations_with_slope_test_Revision()
    % ---------------------------
    % Load VAS
    % ---------------------------
    filename = '/Volumes/ExtremeSSDD/LIBR_tACS/VAS/Craving_Final.xlsx';
    data = readtable(filename);

    % ---------------------------
    % Load gPPI (CONN) matrices
    % Z is [nROI x nROI x nSubjects] in CONN ROI results
    % ---------------------------
    load('/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_BNA/resultsROI_Condition001.mat'); Opioid_Pre  = Z;
    load('/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_BNA/resultsROI_Condition002.mat'); Neutral_Pre = Z;
    load('/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_BNA/resultsROI_Condition003.mat'); Opioid_Post = Z;
    load('/Volumes/ExtremeSSDD/LIBR_tACS_CONN/conn_project01_LIBR_tACS/results/firstlevel/gPPI_BNA/resultsROI_Condition004.mat'); Neutral_Post= Z;

    OvN_Pre    = Opioid_Pre  - Neutral_Pre;
    OvN_Post   = Opioid_Post - Neutral_Post;
    OvN_Change = OvN_Post - OvN_Pre;   % [ROI x ROI x subj]
    Conn_var   = OvN_Change;

    % ---------------------------
    % Group info (match your coding)
    % ---------------------------
    % data.Group contains 'A' (Sham) and 'B' (Active)
    groups = {'A','B'};
    groupLabels = {'Sham','Active'};
    groupColorsHex = {'#F8766D', '#00BFC4'};
    groupColors = cellfun(@hex2rgb, groupColorsHex, 'UniformOutput', false);

    % ---------------------------
    % ROI labels (same as yours)
    % ---------------------------
    regionsOfInterest = [13, 14, 211, 212, 213, 214, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258];
    regionNames = {
        'Left VMPFC', 'Right VMPFC', 'Left medial Amyg', 'Right medial Amyg', 'Left lateral Amyg',...
        'Right lateral Amyg', 'Bilateral VMPFC', 'Right Parietal', 'Right Frontal', 'Left Parietal',...
        'Left Frontal', 'rVS', 'Striatum_subregions.STR1_r', 'Striatum_subregions.STR1_l', ...
        'Striatum_subregions.STR2_r', 'Striatum_subregions.STR2_l', 'Striatum_subregions.STR3_r',...
        'Striatum_subregions.STR3_l'
    };

    % ROI pairs (same as yours)
    roiPairs = [
        249, 256;
        247, 250;
        248, 13;
        248, 212;
        248, 251;
        248, 258;
        251, 211;
        251, 258;
        252, 257;
        252, 258;
        253, 257;
        254, 214
    ];

    % Which VAS column to use (you used 'Change')
    timePoint = 'Change';  % must be a column name in the Craving_Final.xlsx

    % ---------------------------
    % IMPORTANT: subject alignment
    % ---------------------------
    % Assumption: CONN subject dimension matches the subject order in your
    % Craving_Final.xlsx. If not, you MUST add an explicit mapping here.
    nConnSubj = size(Conn_var, 3);
    if height(data) ~= nConnSubj
        warning('Row count in VAS table (%d) does not match CONN subjects (%d). You may need an ID->index mapping.', height(data), nConnSubj);
    end

    % Build group categorical for regression
    % Baseline group will be Sham ('A') unless you reorder categories
    GroupCat = categorical(data.Group, groups, groupLabels); % A->Sham, B->Active

    % Output table to save interaction tests
    results = table();

    for k = 1:size(roiPairs,1)
        roi1 = roiPairs(k,1);
        roi2 = roiPairs(k,2);

        % Extract connectivity for ALL subjects in the same row order as "data"
        connAll = squeeze(Conn_var(roi1, roi2, :));  % [nSubj x 1]
        vasAll  = data.(timePoint);

        % Remove NaNs (both conn and vas)
        valid = ~isnan(connAll) & ~isnan(vasAll) & ~isundefined(GroupCat);
        connAll = connAll(valid);
        vasAll  = vasAll(valid);
        grpAll  = GroupCat(valid);

        % Need at least some subjects per group
        if numel(connAll) < 6 || numel(unique(grpAll)) < 2
            fprintf('Skipping ROI %d-%d: not enough valid data or only one group present.\n', roi1, roi2);
            continue;
        end

        % ---------------------------
        % Reviewer-requested test:
        % One model with interaction (difference in slopes)
        % ---------------------------
        T = table(connAll, grpAll, vasAll, 'VariableNames', {'Conn','Group','VAS'});
        mdl = fitlm(T, 'VAS ~ Conn*Group');  % includes Conn, Group, Conn:Group

        % Get interaction p-value
        coef = mdl.Coefficients;
        % Term name will look like 'Conn:Group_Active' depending on MATLAB version
        termNames = coef.Properties.RowNames;
        intIdx = find(contains(termNames,'Conn:Group'), 1);
        if isempty(intIdx)
            error('Could not find interaction term in model coefficients. Terms: %s', strjoin(termNames', ', '));
        end
        pInteraction = coef.pValue(intIdx);
        betaInteraction = coef.Estimate(intIdx);

        % Optional: within-group correlations/slopes (to match your plot labels)
        rVals = nan(2,1); pVals = nan(2,1);
        slopes = nan(2,1);

        % Plot
        figure; hold on;

        for g = 1:2
            thisLabel = groupLabels{g};
            thisColor = groupColors{g};

            idxG = (grpAll == thisLabel);
            x = connAll(idxG);
            y = vasAll(idxG);

            if numel(x) > 1
                [rVals(g), pVals(g)] = corr(x(:), y(:), 'Type','Pearson');

                % Within-group slope (simple lm)
                mdlG = fitlm(x, y);
                slopes(g) = mdlG.Coefficients.Estimate(2);

                scatter(x, y, 36, 'MarkerFaceColor', thisColor, 'MarkerEdgeColor', thisColor);

                x_fit = linspace(min(x), max(x), 100)';
                [y_fit, y_ci] = predict(mdlG, x_fit);

                plot(x_fit, y_fit, 'Color', thisColor, 'LineWidth', 1.5);
                fill([x_fit; flipud(x_fit)], [y_ci(:,1); flipud(y_ci(:,2))], thisColor, ...
                    'FaceAlpha', 0.10, 'EdgeColor','none');
            end
        end

        grid on; grid minor;

        roi1Name = regionNames{regionsOfInterest == roi1};
        roi2Name = regionNames{regionsOfInterest == roi2};

        xlabel(sprintf('Δ gPPI (Opioid–Neutral) Post–Pre: %s ↔ %s', roi1Name, roi2Name), 'Interpreter','none');
        ylabel(sprintf('VAS %s', timePoint), 'Interpreter','none');

        title(sprintf('ROI %d–%d | Slope difference test: p_{int}=%.3g', roi1, roi2, pInteraction), 'Interpreter','none');

        % Text box: within-group r/p + interaction p
        xL = xlim; yL = ylim;
        xText = xL(1) + 0.03*(xL(2)-xL(1));
        yText = yL(2) - 0.05*(yL(2)-yL(1));
        text(xText, yText, sprintf('Interaction (Conn×Group): \\beta=%.3f, p=%.3g', betaInteraction, pInteraction), ...
            'FontSize', 11, 'BackgroundColor','white');

        if ~isnan(rVals(1))
            text(xText, yText - 0.08*(yL(2)-yL(1)), sprintf('Sham: r=%.2f, p=%.3g', rVals(1), pVals(1)), ...
                'Color', groupColors{1}, 'FontSize', 11, 'BackgroundColor','white');
        end
        if ~isnan(rVals(2))
            text(xText, yText - 0.16*(yL(2)-yL(1)), sprintf('Active: r=%.2f, p=%.3g', rVals(2), pVals(2)), ...
                'Color', groupColors{2}, 'FontSize', 11, 'BackgroundColor','white');
        end

        hold off;

        % Save results row
        newRow = table(roi1, roi2, string(roi1Name), string(roi2Name), ...
            betaInteraction, pInteraction, ...
            slopes(1), slopes(2), rVals(1), pVals(1), rVals(2), pVals(2), ...
            'VariableNames', {'ROI1','ROI2','ROI1_Name','ROI2_Name', ...
                              'Beta_Interaction','P_Interaction', ...
                              'Slope_Sham','Slope_Active', ...
                              'r_Sham','p_Sham','r_Active','p_Active'});
        results = [results; newRow]; %#ok<AGROW>
    end

    % Write out a summary table (recommended for your revision)
    outCSV = fullfile(pwd, sprintf('SlopeDifference_ConnxGroup_%s.csv', timePoint));
    writetable(results, outCSV);
    fprintf('\nSaved interaction-test summary to: %s\n', outCSV);
end


function rgb = hex2rgb(hex)
    hex = char(hex);
    if hex(1) == '#', hex = hex(2:end); end
    rgb = reshape(sscanf(hex, '%2x')/255, 1, 3);
end
