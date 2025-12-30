%% ============================================================
%  plot_YellowBorder_timeseries_perSubject.m
%
%  Uses per-subject MAT files created earlier:
%    SubjectXX_Pre[_ses-*]_YellowBorder.mat
%    SubjectXX_Post[_ses-*]_YellowBorder.mat
%
%  Each MAT should contain struct S with:
%    S.Onsets_sec, S.Durations_sec
%
%  Outputs:
%    - One plot per subject per timepoint
% ============================================================

clear; clc;

%% ---------------- USER SETTINGS ----------------
inDir  = '/Volumes/ExtremeSSDD/LIBR_tACS/InsideScannerData/YellowBorder_PerSubject'; % where SubjectXX_*.mat files are
outDir = fullfile(inDir, 'Plots_TimeSeries');
if ~isfolder(outDir), mkdir(outDir); end

% Choose time resolution for plotting (seconds)
dt = 0.05;  % 50 ms resolution looks smooth; you can set dt=TR if you want TR grid

% Optional: force a fixed plot window (seconds). Leave [] to auto per subject
xlimSec = [];  % e.g., [0 600]

% Save?
savePNGs = true;
dpi = 200;

% Timepoints to plot
timepoints = {'Pre','Post'};

%% ---------------- Find subject files ----------------
matFiles = dir(fullfile(inDir, 'Subject*_*.mat'));
if isempty(matFiles)
    error('No Subject*.mat files found in %s', inDir);
end

% Group by SubjectXX and Timepoint
for t = 1:numel(timepoints)
    tp = timepoints{t};

    tpFiles = dir(fullfile(inDir, sprintf('Subject*_%s*_YellowBorder.mat', tp)));
    if isempty(tpFiles)
        warning('No files found for %s', tp);
        continue
    end

    for f = 1:numel(tpFiles)
        filePath = fullfile(tpFiles(f).folder, tpFiles(f).name);
        D = load(filePath);

        if ~isfield(D,'S')
            warning('Skipping (no struct S): %s', filePath);
            continue
        end
        S = D.S;

        on  = S.Onsets_sec(:);
        dur = S.Durations_sec(:);

        % Remove invalid
        good = isfinite(on) & isfinite(dur) & dur>=0;
        on = on(good);
        dur = dur(good);

        % If empty, still plot an empty panel
        if isempty(on)
            tEnd = 10;
            tVec = 0:dt:tEnd;
            y = zeros(size(tVec));
        else
            tStart = max(0, min(on) - 2);
            tEnd   = max(on + max(dur,0.1)) + 2;
            tVec   = tStart:dt:tEnd;
            y      = buildBinaryTS(tVec, on, dur);
        end

        % ---- Plot ----
        fig = figure('Color','w', 'Position',[100 100 1100 250]); %#ok<NASGU>
        plot(tVec, y, 'LineWidth', 1.5);
        ylim([-0.1 1.1]);
        yticks([0 1]);
        yticklabels({'0','1'});
        ylabel('Yellow border (0/1)');
        xlabel('Time (s)');

        ttl = sprintf('%s | %s | %s', S.SubjectCode, S.Timepoint, S.SessionID);
        title(ttl, 'Interpreter','none');
        grid on; box on;

        if ~isempty(xlimSec)
            xlim(xlimSec);
        end

        % Add onset markers as vertical lines (helps readability)
        hold on;
        for i = 1:numel(on)
            xline(on(i), '--');
        end

        % Save
        if savePNGs
            [~,base,~] = fileparts(filePath);
            outPng = fullfile(outDir, [base '_TimeSeries.png']);
            exportgraphics(gcf, outPng, 'Resolution', dpi);
        end

        close(gcf);
    end
end

fprintf('Done. Plots saved in:\n  %s\n', outDir);

%% ===================== Helper =====================
function y = buildBinaryTS(tVec, on, dur)
% Build y(t)=1 during [onset, onset+duration], else 0
y = zeros(size(tVec));
for k = 1:numel(on)
    t0 = on(k);
    t1 = on(k) + dur(k);
    if t1 <= t0
        % if duration is 0, make it a single-sample pulse
        [~,ix] = min(abs(tVec - t0));
        y(ix) = 1;
    else
        y = y | (tVec >= t0 & tVec <= t1);
    end
end
y = double(y);
end
