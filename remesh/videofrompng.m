

% pngs_to_video_padded.m
clear; clc;
% --- Settings ---
inputFolder = './render/multicomp';
outputVideo = './render/multicomp';
frameRate   = 10;
bgColor     = 255;   % 255 = white, 0 = black

files = dir(fullfile(inputFolder, '*.png'));
if isempty(files)
    error('No PNG files found in %s', inputFolder);
end

% Sort by numeric filename: 1.png, 2.png, 10.png, ...
nums = nan(numel(files),1);
for k = 1:numel(files)
    [~, name, ~] = fileparts(files(k).name);
    nums(k) = str2double(name);
end

if any(isnan(nums))
    error('All filenames must be integers like 1.png, 2.png, 15.png');
end

[~, idx] = sort(nums);
files = files(idx);

% First pass: find max frame size
maxH = 0;
maxW = 0;

for k = 1:numel(files)
    img = imread(fullfile(files(k).folder, files(k).name));
    [h, w, ~] = size(img);
    maxH = max(maxH, h);
    maxW = max(maxW, w);
end

% Some video profiles prefer even dimensions
maxH = maxH + mod(maxH, 2);
maxW = maxW + mod(maxW, 2);

vw = VideoWriter(outputVideo, 'Motion JPEG AVI');
vw.FrameRate = frameRate;
open(vw);

for k = 1:numel(files)
    img = imread(fullfile(files(k).folder, files(k).name));

    if ndims(img) == 2
        img = repmat(img, [1 1 3]);
    end

    [h, w, ~] = size(img);

    padTop  = floor((maxH - h)/2);
    padLeft = floor((maxW - w)/2);

    padded = uint8(bgColor * ones(maxH, maxW, 3));
    padded(padTop+1:padTop+h, padLeft+1:padLeft+w, :) = img;

    writeVideo(vw, padded);
end

close(vw);
fprintf('Video written to: %s\n', outputVideo);