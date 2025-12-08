%% ============================================================
%                  CLEAN START
% ============================================================
restoredefaultpath;
rehash toolboxcache;

addpath(genpath('/Users/alexlipman/Documents/github/isetbio'));
addpath(genpath('/Users/alexlipman/Documents/github/isetcam'));
ieInit;

%% ============================================================
%          DIRECTORY OF AUTO-RUN PNG EXPORTS
% ============================================================
imgDir = '../../Downloads/inputs';   % <-- ALL PNGs go here

%% ============================================================
%                  PARAMETERS
% ============================================================
TARGET_FRAMES = 80;     % final frame count per pitch
targetN       = 128;     % resize image to NxN
nSub          = 4;       % motion-blur subframes

%% ============================================================
%                  ISETBIO OBJECTS
% ============================================================
oi = oiCreate;
cm = cMosaic('sizeDegs', [0.5 0.5]);
cm.integrationTime = 0.03;
cm.noiseFlag = 'none';

%% ============================================================
%                  MASTER DATA TABLE
% ============================================================
masterTable = table();

%% ============================================================
%             FIND & PARSE ALL PNG FRAME NAMES
% ============================================================
files = dir(fullfile(imgDir, '*.png'));

if isempty(files)
    error('No PNG files found in %s', imgDir);
end

fprintf('Found %d PNG frames.\n', numel(files));

% Example filename:
% pitchType_angle_speed_spin_run3_frame012.png
%
% Regex groups:
%   1) pitch type
%   2) angle (with deg)
%   3) speed
%   4) spin
%   5) run #
%   6) frame #
pattern = '([A-Za-z0-9]+)_([-\+A-Za-z0-9]+deg)_([A-Za-z0-9\+\-]+)_([A-Za-z0-9\+\-]+)_run(\d+)_frame(\d+)\.png';

pitchRuns = containers.Map;

%% === Build structure of runs: runKey → list of filenames ===
for i = 1:numel(files)
    fname = files(i).name;
    tok = regexp(fname, pattern, 'tokens', 'once');

    if isempty(tok)
        fprintf('Skipping unrecognized file: %s\n', fname);
        continue;
    end

    pitchType   = tok{1};
    angleStr    = tok{2};
    speedStr    = tok{3};
    spinStr     = tok{4};
    runID       = sprintf('run%s', tok{5});
    frameNum    = str2double(tok{6});

    % Unique group key:
    runKey = sprintf('%s__%s__%s__%s__%s', ...
        pitchType, angleStr, speedStr, spinStr, runID);

    if ~isKey(pitchRuns, runKey)
        pitchRuns(runKey) = struct( ...
            'pitch', pitchType, ...
            'angle', angleStr, ...
            'speed', speedStr, ...
            'spin', spinStr, ...
            'runID', runID, ...
            'frames', [], ...
            'fullpaths', [] );
    end

    entry = pitchRuns(runKey);
    entry.frames(end+1) = frameNum;
    entry.fullpaths{end+1} = fullfile(imgDir, fname);
    pitchRuns(runKey) = entry;
end

allKeys = pitchRuns.keys;
N = 3;   % <-- test first 3 runs
% allKeys = allKeys(1:min(N, numel(allKeys)));
fprintf('Detected %d unique pitch runs.\n', numel(allKeys));

%% ============================================================
%        PROCESS EACH RUN AS IF IT WERE A "PITCH FOLDER"
% ============================================================
for k = 1:numel(allKeys)
    key = allKeys{k};
    entry = pitchRuns(key);

    pitchName = sprintf('%s_%s_%s_%s_%s', ...
        entry.pitch, entry.angle, entry.speed, entry.spin, entry.runID);

    fprintf('\n\n########################################\n');
    fprintf('### PROCESSING: %s\n', pitchName);
    fprintf('########################################\n');

    % Sort frames in temporal order
    [sortedFrames, sortIdx] = sort(entry.frames);
    frameFiles = entry.fullpaths(sortIdx);

    nFramesActual = numel(frameFiles);

    %% -------- Normalize frame count ----------
    frameIdx = linspace(1, nFramesActual, TARGET_FRAMES);
    frameIdx = round(frameIdx);
    frameIdx = max(frameIdx, 1);
    frameIdx = min(frameIdx, nFramesActual);

    absorptions = cell(1, TARGET_FRAMES);

    %% ============================================================
    %         PROCESS EACH NORMALIZED FRAME
    % ============================================================
    for fi = 1:TARGET_FRAMES
        srcIdx = frameIdx(fi);

        fprintf('Frame %d/%d  (source %d/%d)\n', ...
            fi, TARGET_FRAMES, srcIdx, nFramesActual);

        % --- Motion blur accumulation ---
        accum = zeros(targetN, targetN, 3);

        for s = 0:(nSub - 1)
            idx = min(srcIdx + s, nFramesActual);

            imgTemp = imread(frameFiles{idx});
            imgTemp = imresize(imgTemp, [targetN targetN]);

            accum = accum + double(imgTemp);
        end

        imgBlur = uint8(accum / nSub);

        % --- Image → Scene → OI → Cones ---
        sc = sceneFromFile(imgBlur, 'rgb', [], displayCreate('LCD-Apple'));
        sc = sceneSet(sc, 'fov', 2.0);

        oi = oiCompute(oi, sc);
        absorptions{fi} = cm.compute(oi);
    end

    %% ============================================================
    %                    PCA COMPRESSION (5 PCs)
    % ============================================================
    
    % Build matrix: [TARGET_FRAMES × nPixels]
    M = zeros(TARGET_FRAMES, numel(absorptions{1}));
    for ii = 1:TARGET_FRAMES
        M(ii,:) = absorptions{ii}(:)';
    end
    
    % PCA
    [coeff, score, ~, ~, explained] = pca(M, 'NumComponents', 5);
    
    % score → each row = frame, each column = PC#
    PC1 = score(:,1);
    PC2 = score(:,2);
    PC3 = score(:,3);
    PC4 = score(:,4);
    PC5 = score(:,5);
    
    fprintf('PCA explained variance (first 5 PCs): %.2f%%\n', sum(explained(1:5)));


    %% ============================================================
    %            METRICS FOR THIS RUN
    % ============================================================
    totalActivity = zeros(TARGET_FRAMES, 1);
    spreadActivity = zeros(TARGET_FRAMES, 1);
    meanActivity = zeros(TARGET_FRAMES, 1);

    for iFrame = 1:TARGET_FRAMES
        A = absorptions{iFrame};
        totalActivity(iFrame) = sum(A(:));
        spreadActivity(iFrame) = std(A(:));
        meanActivity(iFrame) = mean(A(:));
    end

        %% ============================================================
    %            BUILD TABLE ROWS (resolution-robust)
    % ============================================================
    rows = {};

    for iFrame = 1:TARGET_FRAMES
        A = absorptions{iFrame};
        flat = A(:)';

        % ---- Get size of cone mosaic frame ----
        [H, W] = size(A);

        % ---- Center of the mosaic ----
        cx = round(W/2);
        cy = round(H/2);

        % ---- Fovea: ~20×20 or smaller if needed ----
        foveaHalf = 10;  % desired half-size (20×20)
        fRowStart = max(cy - foveaHalf + 1, 1);
        fRowEnd   = min(cy + foveaHalf, H);
        fColStart = max(cx - foveaHalf + 1, 1);
        fColEnd   = min(cx + foveaHalf, W);

        % ---- Mid: bigger square around fovea ----
        midHalf = round(min(H, W) * 0.25);  % ~50% of min dimension
        mRowStart = max(cy - midHalf + 1, 1);
        mRowEnd   = min(cy + midHalf, H);
        mColStart = max(cx - midHalf + 1, 1);
        mColEnd   = min(cx + midHalf, W);

        % ---- Build masks instead of slicing directly ----
        maskFovea = false(H, W);
        maskMid   = false(H, W);

        % fovea mask
        maskFovea(fRowStart:fRowEnd, fColStart:fColEnd) = true;

        % mid mask (big square)
        maskMid(mRowStart:mRowEnd, mColStart:mColEnd) = true;
        % remove fovea from mid
        maskMid(maskFovea) = false;

        % far = everything else
        maskFar = ~(maskFovea | maskMid);

        % ---- Region totals ----
        fovea_total = sum(A(maskFovea));
        mid_total   = sum(A(maskMid));
        far_total   = sum(A(maskFar));

        % ---- Temporal stuff (same as before) ----
        if iFrame == 1
            delta_total = NaN;
            corr_prev = NaN;
        else
            delta_total = totalActivity(iFrame) - totalActivity(iFrame-1);
            corr_prev = corr(flat', absorptions{iFrame-1}(:));
        end

        % ---- Build table row (unchanged schema) ----
        row = table( ...
            string(entry.pitch), string(entry.angle), ...
            string(entry.speed), string(entry.spin), ...
            string(entry.runID), ...
            iFrame, frameIdx(iFrame), ...
            totalActivity(iFrame), spreadActivity(iFrame), meanActivity(iFrame), ...
            fovea_total, mid_total, far_total, ...
            delta_total, corr_prev, ...
            PC1(iFrame), PC2(iFrame), PC3(iFrame), PC4(iFrame), PC5(iFrame), ...
            'VariableNames', { ...
                'pitch','angle','speed','spin','run', ...
                'frame','src_frame', ...
                'total','spread','meanSignal', ...
                'fovea_total','mid_total','far_total', ...
                'delta_total','corr_prev', ...
                'PC1','PC2','PC3','PC4','PC5' ...
            });


        rows{end+1} = row;
    end


    pitchTable = vertcat(rows{:});
    masterTable = [masterTable; pitchTable];
end

%% ============================================================
%                  SAVE OUTPUT
% ============================================================
fprintf('\nSaving final dataset...\n');

writetable(masterTable, 'allPitchCones.csv');
save('allPitchCones.mat', 'masterTable', '-v7.3');

fprintf('DONE. Saved %d total rows.\n', height(masterTable));
