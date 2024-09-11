%% Plot motion results
%
% Run this script after motion_tracking.m on an image sequence to generate
% figures of average displacement.
%
% Author: Jarrah Dowrick
% Date: 12th Sept 2024
%

%% Set parameters
% Define search information
numGlitter_X = 5;
numGlitter_Y = 3;
grid_dim = [numGlitter_X,numGlitter_Y];

%% Load motion tracking information


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
