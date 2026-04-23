clear; clc;

% --- Parameters ---
inputFolder = "./render/pnd_0.20_video_redux";              % folder containing 1.png, 2.png, etc.
outputFolder = "./render/pnd_0.20_video_redux";      % folder to save resized images
overwrite = false;                      % set true if you want to overwrite originals

% --- Create output folder if needed ---
if ~overwrite && ~exist(outputFolder, "dir")
    mkdir(outputFolder);
elseif overwrite
    outputFolder = inputFolder;  % write back to same folder
end

% --- Get list of PNG files ---
files = dir(fullfile(inputFolder, "phase*.png"));
if isempty(files)
    error("No PNG files found in the specified folder.");
end

% Sort numerically (important if names are 1.png, 2.png, ... 10.png)
[~, order] = sort(arrayfun(@(f) str2double(erase(f.name, ".png")), files));
files = files(order);

% --- Read first image to get target size ---
refImage = imread(fullfile(inputFolder, files(1).name));
targetSize = size(refImage);

fprintf("Reference image: %s (%d x %d)\n", files(1).name, targetSize(1), targetSize(2));

% --- Loop through all images ---
for k = 1:numel(files)
    filename = files(k).name;
    filepath = fullfile(inputFolder, filename);
    img = imread(filepath);

    % Resize if not already matching
    if ~isequal(size(img,1), targetSize(1)) || ~isequal(size(img,2), targetSize(2))
        imgResized = imresize(img, [targetSize(1), targetSize(2)]);
        fprintf("Resized: %s → [%d x %d]\n", filename, targetSize(1), targetSize(2));
    else
        imgResized = img;
    end

    % Save resized image
    outputPath = fullfile(outputFolder, filename);
    imwrite(imgResized, outputPath);
end

disp("✅ All images resized successfully.");
