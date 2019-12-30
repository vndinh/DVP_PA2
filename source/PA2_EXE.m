clc; clear; close all;

% Parameter
width = 352;
height = 288;
YUV_type = [1, 0.5, 0.5];
blk_size = 1;
num_iters = 16;
sigma = 3.3;

% Read yuv file
f_name = '..\data\Calendar_CIF30.yuv';
f_id = fopen(f_name, 'r');
fr_1 = YUV_READER(f_id, width, height, YUV_type, 10, 1);
fr_2 = YUV_READER(f_id, width, height, YUV_type, 11, 1);

psnr_10_11 = psnr(fr_1, fr_2, 255);

%% Motion Estimation using Horn-Schunck algorithm
[v1_hs, v2_hs] = ME_HS(fr_1, fr_2, blk_size);
rec_hs = MC(fr_1, v1_hs, v2_hs);
psnr_hs = psnr(rec_hs, fr_2, 255);

%%
mean_number = 8;
psnr_map = zeros(mean_number, 1);
psnr_occ = zeros(mean_number, 1);
psnr_line = zeros(mean_number, 1);
psnr_all = zeros(mean_number, 1);

for i = 1:mean_number
  %% (a) MAP without occlusion and line field
  [v1_map, v2_map] = MAP(fr_1, fr_2, v1_hs, v2_hs, num_iters, sigma);
  rec_map = MC(fr_1, v1_map, v2_map);
  psnr_map(i) = psnr(rec_map, fr_2, 255);
  
  %% (b) MAP with occlusion field
  [v1_occ, v2_occ] = MAP_OCC(fr_1, fr_2, v1_hs, v2_hs, num_iters, sigma);
  rec_occ = MC(fr_1, v1_occ, v2_occ);
  psnr_occ(i) = psnr(rec_occ, fr_2, 255);
  
  %% (c) MAP with line field
  [v1_line, v2_line] = MAP_LINE(fr_1, fr_2, v1_hs, v2_hs, num_iters, sigma);
  rec_line = MC(fr_1, v1_line, v2_line);
  psnr_line(i) = psnr(rec_line, fr_2, 255);
  
  %% (d) MAP with occlusion and line field
  [v1_all, v2_all] = MAP_ALL(fr_1, fr_2, v1_hs, v2_hs, num_iters, sigma);
  rec_all = MC(fr_1, v1_all, v2_all);
  psnr_all(i) = psnr(rec_all, fr_2, 255);
end

%% PSNR
psnr_map_mean = mean(psnr_map);
psnr_occ_mean = mean(psnr_occ);
psnr_line_mean = mean(psnr_line);
psnr_all_mean = mean(psnr_all);

%% Display
fprintf('Horn-Schunck algorithm: PSNR = %2.4f dB\n', psnr_hs);
fprintf('MAP without occlusion and line field: PSNR = %2.4f dB\n', psnr_map_mean);
fprintf('MAP with occlusion field: PSNR = %2.4f dB\n', psnr_occ_mean);
fprintf('MAP with line field: PSNR = %2.4f dB\n', psnr_line_mean);
fprintf('MAP with occlusion and line field: PSNR = %2.4f dB\n', psnr_all_mean);

figure('Name', 'Reconstruction');
subplot(2,3,1); imshow(uint8(fr_2)); title('11th Frame');
subplot(2,3,2); imshow(uint8(rec_hs)); title('Horn-Schunck');
subplot(2,3,3); imshow(uint8(rec_map)); title('MAP');
subplot(2,3,4); imshow(uint8(rec_occ)); title('MAP with occlusion field');
subplot(2,3,5); imshow(uint8(rec_line)); title('MAP with line field');
subplot(2,3,6); imshow(uint8(rec_all)); title('MAP with occlusion and line field');

figure('Name', 'Difference between 11th and reconstruction frame');
subplot(2,3,1); imshow(uint8(abs(fr_2-fr_1))); title('Difference between 10th and 11th frame');
subplot(2,3,2); imshow(uint8(abs(fr_2-rec_hs))); title('Horn-Schunck');
subplot(2,3,3); imshow(uint8(abs(fr_2-rec_map))); title('MAP');
subplot(2,3,4); imshow(uint8(abs(fr_2-rec_occ))); title('MAP with occlusion field');
subplot(2,3,5); imshow(uint8(abs(fr_2-rec_line))); title('MAP with line field');
subplot(2,3,6); imshow(uint8(abs(fr_2-rec_all))); title('MAP with occlusion and line field');


