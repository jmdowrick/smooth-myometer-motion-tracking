%% Smooth myometer displacement tracking

%% Set parameters
framerate = 25; % frames per second
dt = 1/framerate;
search_rad = 10; % pixels (radius around glitter points to average)

% Visualisation and saving toggle on/off
visualiseProgress = false; % preview tracking throughout analysis
playAnimation = false;
saveAnimation = false;
playAverageAnimation = false;
playInterpAnimation = false;

% Test features on / off
interpolatedField = false; % doesn't work well without dense coverage of glitter

% Define paths for image files and for storing outputs
folder_src = '/Users/jdow403/Desktop/AWB015_VID006';
output_file_suffix = '_tracked';

if (~isfolder(folder_src) || numel(dir(fullfile(folder_src,'*.Bmp'))) == 0)
    folder_src = uigetdir();
end

if numel(dir(fullfile(folder_src,'*.Bmp'))) == 0
    error("Selected folder does not contain an image sequence. Please select a different folder.")
end

if exist(fullfile(folder_src,'outputs'), 'dir')
    % check if processing has already been performed
    if numel(dir(fullfile(folder_src,'outputs','*.mat'))) > 0
        % ask user if they want to repeat analysis
        response = input("Motion tracking output detected, are you sure you want to proceed? (y/n) \n", "s");
        if (lower(response) ~= "y" || isempty(response))
            disp("Exiting motion tracking")
            return
        else
            disp("Proceeding with motion tracking")
        end
    end
else
    % If folder does not exist, create folder.
    mkdir(fullfile(folder_src, 'outputs'))
end

folder_outputs = fullfile(folder_src, 'outputs');

%% Load images 
dirList = dir(fullfile(folder_src,'*.Bmp'));

T = struct2table(dirList);
dirList = natsortrows(T);
dirList = table2struct(dirList);

%% Select region of interest
refImage = imread(fullfile(dirList(1).folder, dirList(1).name));
h = size(refImage,1);
w = size(refImage,2);

figure;
set(gcf, 'Name', 'Select region of interest and double click to confirm', 'NumberTitle', 'off');

% User selects region of interest (encapsulating markers)
[~, roiBox] = imcrop(refImage);
close(gcf)

%% Display region of interest
refImage_roi = insertShape(refImage, 'Rectangle', roiBox, 'Color', 'green');
figure; imshow(refImage_roi);
title("Selected region of interest");

%% Feature tracking
% Detect features to track
points = detectMinEigenFeatures(im2gray(refImage), 'ROI', roiBox);
pts = round(points.Location);
idx = sub2ind([h w], pts(:,2), pts(:,1));

% Create point tracker
pointTracker = vision.PointTracker('MaxBidirectionalError',0.025);
initialize(pointTracker, points.Location, refImage);

% Visualise tracking as it processes
w = size(refImage,1);
h = size(refImage,2);

% Initialise data storage
time = NaN(length(dirList), 1);
nPoints = size(pts, 1);

x = NaN(nPoints, length(dirList)); y = x; V = x;

if visualiseProgress
    videoPlayer = vision.VideoPlayer();
end

f = waitbar(0, 'Tracking motion...', ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)'); % progress bar
setappdata(f, 'canceling', 0);

% For loop through images in folder of interest
for imNum = 1:length(dirList)
    if getappdata(f, 'canceling')
        release(pointTracker)
        if visualiseProgress
            release(videoPlayer)
        end
        delete(f)
        disp('Tracking cancelled \n')
        return
    end

    % Update waitbar
    waitbar(imNum/length(dirList), f)

    % Load image
    currentImage = imread(fullfile(dirList(imNum).folder, dirList(imNum).name));

    [points, isFound, scores] = pointTracker(currentImage);

    x(:,imNum) = points(:,1);
    y(:,imNum) = points(:,2);
    V(:,imNum) = scores > 0.99;
    time(imNum) = (imNum - 1)*dt;

    if any(isFound)
        isFound = scores > 0.99;

        pts = round(points(isFound,:));
        idx = sub2ind([h w], pts(:,2), pts(:,1));
        currentImage = im2gray(currentImage);

        currentImage = insertMarker(currentImage, pts, '+', 'Color', 'green');
    end

    if visualiseProgress
        videoPlayer(currentImage)
    end
end
release(pointTracker)

if visualiseProgress
    release(videoPlayer)
end

% Remove entries that were ever invalid
x = x(~any(V==0,2),:);
y = y(~any(V==0,2),:);

delete(f)

%% Play animation of tracking
if playAnimation
    videoPlayer = vision.VideoPlayer();

    for imNum = 1:length(dirList)
        currentImage = imread(fullfile(dirList(imNum).folder, dirList(imNum).name));
        pts = [x(:,imNum), y(:,imNum)];
        currentImage = insertMarker(im2gray(currentImage), pts, '+', 'Color', 'green');

        videoPlayer(currentImage)
    end
    release(videoPlayer)
end

%% Save animation
if saveAnimation
    video_object = VideoWriter(fullfile(folder_outputs,output_file_suffix),'MPEG-4'); % add _tracked to end of folder name
    open(video_object);
    for imNum = 1:10:length(dirList)
        currentImage = imread(fullfile(dirList(imNum).folder, dirList(imNum).name));
        pts = [x(:,imNum), y(:,imNum)];
        currentImage = insertMarker(im2gray(currentImage), pts, '+', 'Color', 'green');
        writeVideo(video_object, currentImage);
    end
    close(video_object);
end

%% Convert from pixels to mm
fig = figure;
imshow(refImage)
title("Draw line equivalent to 10 mm and press Enter");
roi = drawline;
k = waitforbuttonpress;

length_line = sqrt((roi.Position(1,1)-roi.Position(2,1))^2 + (roi.Position(1,2)-roi.Position(2,2))^2);
ppmm = length_line/10; % pixels per mm
close(fig)

%% Manually select glitter as points to track
fig = figure;
currentImage = imread(fullfile(dirList(length(dirList)).folder, dirList(length(dirList)).name));
pts = [x(:,length(dirList)), y(:,length(dirList))];
currentImage = insertMarker(im2gray(currentImage), pts, '+', 'Color', 'green');

imshow(currentImage)
title("Select glitter locations and press Enter");
[glitter_x, glitter_y] = getpts;
close(fig);

%% Group x and y displacement per point of glitter
av_x = NaN(length(glitter_x), size(x,2));
av_y = av_x;

% I don't think this works currently - needs to consider y aspect when
% averaging per glitter point
for i = 1:length(glitter_x)
    temp_x = x(x(:,1)>=(glitter_x(i) - search_rad) & x(:,1) <= (glitter_x(i) + search_rad), :);
    av_x(i, :) = mean(temp_x);

    temp_y = y(y(:,1)>=(glitter_y(i) - search_rad) & y(:,1) <= (glitter_y(i) + search_rad), :);
    av_y(i, :) = mean(temp_y);
end

clear temp_x temp_y

%% Visualise averaged displacements (per piece of glitter)
if playAveragedAnimation
    videoPlayer = vision.VideoPlayer();

    for imNum = 1:length(dirList)
        currentImage = imread(fullfile(dirList(imNum).folder, dirList(imNum).name));
        pts = [av_x(:,imNum), av_y(:,imNum)];
        currentImage = insertMarker(im2gray(currentImage), pts, '+', 'Color', 'green');
        currentImage = insertShape(currentImage, 'circle', [av_x(:,imNum), av_y(:,imNum), repmat(search_rad,[size(av_x,1),1])], 'Color', 'green');
        videoPlayer(currentImage)
    end
    release(videoPlayer)
end

%% Interpolate velocities across rectangle
if interpolatedField
    % Calculate x and y velocities
    [b,g] = sgolay(5,25);

    dx = zeros(size(x));
    dy = zeros(size(y));

    for i = 1:size(x,1)
        dx(i,:) = conv(x(i,:), -1/dt * g(:,2), 'same');
        dy(i,:) = conv(y(i,:), -1/dt * g(:,2), 'same');
    end

    % Interpolation (in development)
    border_size = 50;
    grid_size = 25;

    [Xq, Yq] = meshgrid(round(roiBox(1)+border_size):grid_size:(round(roiBox(1))+round(roiBox(3))-border_size), ...
        round(roiBox(2)+border_size):grid_size:(round(roiBox(1))+round(roiBox(4))-border_size));

    dx_q = zeros(length(Xq(:)), size(x,2));
    dy_q = zeros(length(Yq(:)), size(y,2));

    for i = 20:(size(x, 2)-20)
        F = scatteredInterpolant(x(:,1), y(:,1), dx(:,i));
        Vq = F(Xq, Yq);
        dx_q(:,i) = Vq(:);

        F = scatteredInterpolant(x(:,1), y(:,1), dx(:,i));
        Vq = F(Xq, Yq);
        dy_q(:,i) = Vq(:);
    end

    x_q = cumtrapz(dx_q, 2)*dt + repmat(Xq(:), [1, size(x,2)]);
    y_q = cumtrapz(dy_q, 2)*dt + repmat(Yq(:), [1, size(y,2)]);

    if playInterpAnimation
        % Play animation of interpolated tracking
        videoPlayer = vision.VideoPlayer();

        for imNum = 1:length(dirList)
            currentImage = imread(fullfile(dirList(imNum).folder, dirList(imNum).name));
            pts = [x_q(:,imNum), y_q(:,imNum)];
            currentImage = insertMarker(im2gray(currentImage), pts, '+', 'Color', 'green');

            videoPlayer(currentImage)
        end
        release(videoPlayer)
    end
end
