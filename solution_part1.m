%% Adding mmread as part of path
clc; clear all;
addpath('./MMread');

%% Part a - Getting the video data using mmread
num_frames = 3;
[video_struct, ~] = mmread('./videos/cars.avi', 1:num_frames);
fprintf("Video Structure info: %d x %d \nVideo Rate: %0.2f \n", video_struct.width, video_struct.height, video_struct.rate);
vid_data = video_struct.frames;
req_ht = 120;
req_wd = 240;
vid_gray_data = zeros(req_ht, req_wd, num_frames);
for frame = 1:num_frames
    vid_gray_data(:, :, frame) = im2double(rgb2gray(vid_data(frame).cdata(end-119:end, end-239:end , :)));
    figure;
    imshow(mat2gray(vid_gray_data(:, :, frame)));
%     title(['Frame #', num2str(frame)]);
    axis tight;
    saveas(gcf,'./sample.png');
end

%% Part b - Making the coded snapshot from the read video frames
noise_std = 2;
code_data = randi([0,1], req_ht, req_wd, num_frames);
coded_snapshot = sum(code_data.*vid_gray_data, 3);
noisy_coded_snapshot = coded_snapshot + ((noise_std/255).*randn(req_ht, req_wd));
% noisy_coded_snapshot = coded_snapshot + normrnd(0, noise_std/255, req_ht, req_wd);
figure;
imshow(mat2gray(noisy_coded_snapshot));
title(['Coded snapshot with noise for T = ' num2str(num_frames)]);
saveas(gcf, './Figures/Experiment_1/reconstructed/T_3/coded_snapshot.png');

%% Part c - Reconstruction of the video from the coded snapshot(noisy) for 3 frames
patch_size = 8;
eps = 3*(noise_std/255)*8;
tic;
[reconstructed_x] = omp_reconstruction(code_data, noisy_coded_snapshot, eps, patch_size);
toc;

for frame = 1:num_frames
    figure;
    imshow(mat2gray(reconstructed_x(:, :, frame)));
    title(['Reconstructed frame#' num2str(frame)]);
    axis tight;
    saveas(gcf, ['./Figures/Experiment_1/reconstructed/T_3/', num2str(frame), '.png']);
end

%% RMSE calculation for part c
m = mean(vid_gray_data, 'all');
RMSE = mean((reconstructed_x - vid_gray_data).^2, 'all')/mean((vid_gray_data - m).^2, 'all');
% RMSE = mean((reconstructed_x - vid_gray_data).^2, 'all')/mean((vid_gray_data).^2, 'all');
fprintf("Relative Mean Squared Error for T = %d frame reconstuction: %0.5f\n", num_frames, RMSE);

%% Part d - Reconstruction of the video from the coded snapshot(noisy) for 5 frames
num_frames = 5;
[video_struct, ~] = mmread('./videos/cars.avi', 1:num_frames);
vid_data = video_struct.frames;
vid_gray_data = zeros(req_ht, req_wd, num_frames);
for frame = 1:num_frames
    vid_gray_data(:, :, frame) = im2double(rgb2gray(vid_data(frame).cdata(end-119:end, end-239:end , :)));
    figure;
    imshow(mat2gray(vid_gray_data(:, :, frame)));
    title(['Frame #' num2str(frame)]);
    axis tight;
end

code_data = randi([0,1], req_ht, req_wd, num_frames);
coded_snapshot = sum(code_data.*vid_gray_data, 3);
noisy_coded_snapshot = coded_snapshot + ((noise_std/255).*randn(req_ht, req_wd));
% noisy_coded_snapshot = coded_snapshot + normrnd(0, noise_std/255, req_ht, req_wd);
figure;
imshow(mat2gray(noisy_coded_snapshot));
title(['Coded snapshot with noise for T = ' num2str(num_frames)]);
saveas(gcf, './Figures/Experiment_1/reconstructed/T_5/coded_snapshot.png');

patch_size = 8;
eps = 3*(noise_std/255)*8;
tic;
[reconstructed_x] = omp_reconstruction(code_data, noisy_coded_snapshot, eps, patch_size);
toc;

for frame = 1:num_frames
    figure;
    imshow(mat2gray(reconstructed_x(:, :, frame)));
    title(['Reconstructed frame#' num2str(frame)]);
    axis tight;
    saveas(gcf, ['./Figures/Experiment_1/reconstructed/T_5/', num2str(frame), '.png']);
end

%% RMSE calculation for part d
m = mean(vid_gray_data, 'all');
RMSE = mean((reconstructed_x - vid_gray_data).^2, 'all')/mean((vid_gray_data - m).^2, 'all');
% RMSE = mean((reconstructed_x - vid_gray_data).^2, 'all')/mean((vid_gray_data).^2, 'all');
fprintf("Relative Mean Squared Error for T = %d frame reconstuction: %0.5f\n", num_frames, RMSE);

%% Part e - Reconstruction of the video from the coded snapshot(noisy) for 7 frames
num_frames = 7;
[video_struct, ~] = mmread('./videos/cars.avi', 1:num_frames);
vid_data = video_struct.frames;
vid_gray_data = zeros(req_ht, req_wd, num_frames);

for frame = 1:num_frames
    vid_gray_data(:, :, frame) = im2double(rgb2gray(vid_data(frame).cdata(end-119:end, end-239:end , :)));
    figure;
    imshow(mat2gray(vid_gray_data(:, :, frame)));
    title(['Frame #' num2str(frame)]);
    axis tight;
    saveas(gcf, ['./Figures/Experiment_1/original/' num2str(frame) '.png']);
end

code_data = randi([0,1], req_ht, req_wd, num_frames);
coded_snapshot = sum(code_data.*vid_gray_data, 3);
noisy_coded_snapshot = coded_snapshot + ((noise_std/255).*randn(req_ht, req_wd));
% noisy_coded_snapshot = coded_snapshot + normrnd(0, noise_std/255, req_ht, req_wd);
figure;
imshow(mat2gray(noisy_coded_snapshot));
title(['Coded snapshot with noise for T = ' num2str(num_frames)]);
saveas(gcf, './Figures/Experiment_1/reconstructed/T_7/coded_snapshot.png');

patch_size = 8;
eps = 3*(noise_std/255)*8;
tic;
[reconstructed_x] = omp_reconstruction(code_data, noisy_coded_snapshot, eps, patch_size);
toc;

for frame = 1:num_frames
    figure;
    imshow(mat2gray(reconstructed_x(:, :, frame)));
    title(['Reconstructed frame#' num2str(frame)]);
    axis tight;
    saveas(gcf, ['./Figures/Experiment_1/reconstructed/T_7/', num2str(frame), '.png']);
end

%% RMSE calculation for part e
m = mean(vid_gray_data, 'all');
RMSE = mean((reconstructed_x - vid_gray_data).^2, 'all')/mean((vid_gray_data - m).^2, 'all');
% RMSE = mean((reconstructed_x - vid_gray_data).^2, 'all')/mean((vid_gray_data).^2, 'all');
fprintf("Relative Mean Squared Error for T = %d frame reconstuction: %0.5f\n", num_frames, RMSE);

%% Part f - Reconstruction of the flames video from coded snapshot of T = 5 frames
num_frames = 5;
[video_struct, ~] = mmread('./videos/flame.avi', 1:num_frames);
vid_data = video_struct.frames;
ht = video_struct.height;
wd = video_struct.width;
vid_gray_data = zeros(ht, wd, num_frames);
for frame = 1:num_frames
    vid_gray_data(:, :, frame) = im2double(rgb2gray(vid_data(frame).cdata));
    figure;
    imshow(mat2gray(vid_gray_data(:, :, frame)));
    title(['Frame #' num2str(frame)]);
    axis tight;
    saveas(gcf, ['./Figures/Experiment_2/original/' num2str(frame) '.png']);
end
%%
noise_std = 2;
code_data = randi([0,1], ht, wd, num_frames);
coded_snapshot = sum(code_data.*vid_gray_data, 3);
noisy_coded_snapshot = coded_snapshot + ((noise_std/255).*randn(ht, wd));
% noisy_coded_snapshot = coded_snapshot + normrnd(0, noise_std/255, req_ht, req_wd);
figure;
imshow(mat2gray(noisy_coded_snapshot));
title(['Coded snapshot with noise for T = ' num2str(num_frames)]);
saveas(gcf, './Figures/Experiment_2/reconstructed/coded_snapshot.png');
%%
patch_size = 8;
eps = 3*(noise_std/255)*8;
tic;
[reconstructed_x] = omp_reconstruction(code_data, noisy_coded_snapshot, eps, patch_size);
toc;

for frame = 1:num_frames
    figure;
    imshow(mat2gray(reconstructed_x(:, :, frame)));
    title(['Reconstructed frame #' num2str(frame)]);
    axis tight;
    saveas(gcf, ['./Figures/Experiment_2/reconstructed/', num2str(frame), '.png']);
end

%% RMSE calculation for part e
m = mean(vid_gray_data, 'all');
RMSE = mean((reconstructed_x - vid_gray_data).^2, 'all')/mean((vid_gray_data - m).^2, 'all');
% RMSE = mean((reconstructed_x - vid_gray_data).^2, 'all')/mean((vid_gray_data).^2, 'all');
fprintf("Relative Mean Squared Error for T = %d frame reconstuction: %0.5f\n", num_frames, RMSE);