%% Smooth myometer displacement tracking

%% Set parameters
framerate = 25; % frames per second
dt = 1/framerate;

% Visualisation and saving toggle on/off
playAnimation = true;
saveAnimation = true;
playAverageAnimation = false;

% Test features on / off
interpolatedField = false; % doesn't work well without dense coverage of glitter

% Define search information
numGlitter_X = 5;
numGlitter_Y = 3;
grid_dim = [numGlitter_X,numGlitter_Y];
search_rad = 10; % pixels (radius around glitter points to average)

% Define paths for image files and saved outputs
folder_main = '/Users/jdow403/Desktop/AWB015_VID006';
folder_src = append(folder_main, '/images/');
folder_outputs = append(folder_main, '/outputs');
output_file_prefix = '/trackedcolon';

%% Load images and select region of interest
dirList = dir(append(folder_src,'*.Bmp'));

T = struct2table(dirList);
dirList = natsortrows(T);
dirList = table2struct(dirList);

refImage = imread(append(dirList(1).folder, '/', dirList(1).name));
figure; imshow(refImage)
h = size(refImage,1);
w = size(refImage,2);

% User selects region of interest (encapsulating markers)
[~, roiBox] = imcrop(refImage);
close(gcf)

% Display region of interest
refImage_roi = insertShape(refImage, 'Rectangle', roiBox, 'Color', 'green');
figure; imshow(refImage_roi);
title("Selected Region of Interest");

%% Feature tracking
% Detect features to track
points = detectMinEigenFeatures(im2gray(refImage), 'ROI', roiBox);
figure; imshow(refImage);
hold on; plot(points);
title("Detected Features");

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
videoPlayer = vision.VideoPlayer();

% For loop through images in folder of interest
for imNum = 1:length(dirList)
    % Load image
    currentImage = imread(append(dirList(imNum).folder, '/', dirList(imNum).name));

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

    videoPlayer(currentImage)
end
release(pointTracker)

% Remove entries that were ever invalid
x = x(~any(V==0,2),:);
y = y(~any(V==0,2),:);

%% Play animation of tracking
if playAnimation
    videoPlayer = vision.VideoPlayer();

    for imNum = 1:length(dirList)
        currentImage = imread(append(dirList(imNum).folder, '/', dirList(imNum).name));
        pts = [x(:,imNum), y(:,imNum)];
        currentImage = insertMarker(im2gray(currentImage), pts, '+', 'Color', 'green');

        videoPlayer(currentImage)
    end
    release(videoPlayer)
end

%% Save animation
if saveAnimation
    video_object = VideoWriter(append(folder_outputs,output_file_prefix),'MPEG-4');
    open(video_object);
    for imNum = 1:10:length(dirList)
        currentImage = imread(append(dirList(imNum).folder, '/', dirList(imNum).name));
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
imshow(refImage)
title("Select glitter locations and press Enter");
[glitter_x, glitter_y] = getpts;
close(fig);

%% Group x and y displacement per point of glitter
av_x = NaN(length(glitter_x), size(x,2));
av_y = av_x;
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
        currentImage = imread(append(dirList(imNum).folder, '/', dirList(imNum).name));
        pts = [av_x(:,imNum), av_y(:,imNum)];
        currentImage = insertMarker(im2gray(currentImage), pts, '+', 'Color', 'green');
        currentImage = insertShape(currentImage, 'circle', [av_x(:,imNum), av_y(:,imNum), repmat(search_rad,[size(av_x,1),1])], 'Color', 'green');
        videoPlayer(currentImage)
    end
    release(videoPlayer)
end

%% Calculate average y-displacement and x-displacement in the different axes

av_x_sort = sortrows(av_x, 1);
av_y_sort = sortrows(av_y, 1);

for x_columns = 1:grid_dim(1)
    x_disp(x_columns,:) = mean(av_x_sort((grid_dim(2)*(x_columns-1) + 1):(grid_dim(2)*(x_columns)),:),1);
    x_disp_0(x_columns,:) = x_disp(x_columns,:) - x_disp(x_columns,1);
end

for y_columns = 1:grid_dim(2)
    y_disp(y_columns,:) = mean(av_y_sort((grid_dim(1)*(y_columns-1) + 1):(grid_dim(1)*(y_columns)),:),1);
    y_disp_0(y_columns,:) = y_disp(y_columns,:) - y_disp(y_columns,1);
end

y_disp_legend = y_disp(:,1)./ppmm;
x_disp_legend = x_disp(:,1)./ppmm;

%% Plot displacements - X axis
colors = get(gca,'colororder');
imNum = 1;
currentImage = imread(append(dirList(imNum).folder, '/', dirList(imNum).name));
pts = [av_x(:,imNum), av_y(:,imNum)];

[av_x_sort, idx] = sortrows(av_x, 1);
av_y_sort = av_y(idx,:);
for i = 1:numGlitter_X
    range = (numGlitter_Y*(i-1)+1):(numGlitter_Y*i);
    currentImage = insertShape(currentImage, 'filled-circle', [av_x_sort(range,imNum), ...
        av_y_sort(range,imNum), ...
        repmat(search_rad,[numGlitter_Y,1])],'Opacity', 0.3, 'ShapeColor', colors(i,:), 'LineWidth', 2);
end

figure;
subplot(1,2,1)
hold on
for i = 1:size(x_disp_0,1)
    plot(time(1:end-50),x_disp_0(i,1:(end-50))./ppmm, 'LineWidth', 1.5)
end
xlabel('Time (s)')
ylabel('x displacement (mm)')
subplot(1,2,2)
imshow(currentImage)
fontsize(18,"points")

%% Plot displacements - Y axis
colors = get(gca,'colororder');
imNum = 1;
currentImage = imread(append(dirList(imNum).folder, '/', dirList(imNum).name));
pts = [av_x(:,imNum), av_y(:,imNum)];

[av_y_sort, idx] = sortrows(av_y, 1);
av_x_sort = av_x(idx,:);

for i = 1:numGlitter_Y
    range = (numGlitter_X*(i-1)+1):(numGlitter_X*i);
    currentImage = insertShape(currentImage, 'filled-circle', [av_x_sort(range,imNum), ...
        av_y_sort(range,imNum), ...
        repmat(search_rad,[numGlitter_X,1])],'Opacity', 0.3, 'ShapeColor', colors(i,:), 'LineWidth', 2);
end

figure;
subplot(1,2,1)
hold on
for i = 1:size(y_disp_0,1)
    plot(time(1:end-50),y_disp_0(i,1:end-50)./ppmm, 'LineWidth', 1.5)
end
xlabel('Time (s)')
ylabel('y displacement (mm)')
subplot(1,2,2)
imshow(currentImage)
fontsize(18,"points")

%% Calculate x and y velocities
[b,g] = sgolay(5,25);

dx = zeros(size(x));
dy = zeros(size(y));

for i = 1:size(x,1)
    dx(i,:) = conv(x(i,:), -1/dt * g(:,2), 'same');
    dy(i,:) = conv(y(i,:), -1/dt * g(:,2), 'same');
end

%% Compare Savitzky-Golay derivative with naive differentiation
dx_alt = diff(x,1,2)/dt;
figure; hold on;  plot(dx_alt(5,20:(end-20))); plot(dx(5,20:(end-20)), 'LineWidth', 2);

%% Interpolate velocities across rectangle
if interpolatedField
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

    % Play animation of interpolated tracking
    videoPlayer = vision.VideoPlayer();

    for imNum = 1:length(dirList)
        currentImage = imread(append(dirList(imNum).folder, '/', dirList(imNum).name));
        pts = [x_q(:,imNum), y_q(:,imNum)];
        currentImage = insertMarker(im2gray(currentImage), pts, '+', 'Color', 'green');

        videoPlayer(currentImage)
    end
    release(videoPlayer)

end
