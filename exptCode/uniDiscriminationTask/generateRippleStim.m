function stimTextureMatrix = generateRippleStim(VSinfo,ScreenInfo,windowPtr, x_loc_pixel)
%% Gaussian Atlas
% Parameters

increment = 0.2; %pi/VSinfo.numFrames;
back_sd_y = VSinfo.noise_sd;
stim_sd_y = VSinfo.stim_sd; %doesn't matter when repeating lines
height = VSinfo.height;
wBack = VSinfo.wBack;
wGauss = 1 - wBack;

% Generate the anisotropic Gaussian blob
[x, y] = meshgrid(0:ScreenInfo.xaxis-1, 0:height-1);
center_x = ScreenInfo.xmid;
center_y = height / 2;
mean_luminance = 1.0;
targetLoc = [center_y,ScreenInfo.xmid + x_loc_pixel];

back_gaussian = mean_luminance * exp(-((x - center_x).^2 / (2 * noise_sd_y^2) + (y - center_y).^2 / (2 * back_sd_y^2)));

% Perform Fourier transformation on the anisotropic Gaussian blob
background_gaussian_fft = fftshift(fft2(back_gaussian));

% Extract the amplitude and phase
amplitude = abs(background_gaussian_fft);
initial_phase = rand(size(amplitude)) * 2 * pi;

% Create random + or - 0.2 increment directions
phase_directions = randi([0, 1], size(initial_phase)) * 2 - 1;
phase_directions = phase_directions * increment;

% Initialize the 3D matrix to store all frames
all_frames = zeros(height, ScreenInfo.xaxis, VSinfo.numFrames);

% Generate all frames with updated phases
for frame = 1:VSinfo.numFrames
    initial_phase = initial_phase + phase_directions;
    new_fft = amplitude .* exp(1i * initial_phase);
    new_image_fft_shifted = ifft2(ifftshift(new_fft));
    all_frames(:, :, frame) = abs(new_image_fft_shifted);
end

% Normalize the entire 3D time-space matrix
max_val = max(all_frames(:));
normalized_frames = (all_frames / max_val) * 255;

% Renormalize the Gaussian blob
stim_gaussian = mean_luminance * exp(-((x - targetLoc(2)).^2 / (2 * stim_sd_y^2) + (y - targetLoc(1)).^2 / (2 * stim_sd_y^2)));
renomalized_gaussian_blob = (stim_gaussian / max(stim_gaussian(:))) * 255;
%%
% Pick out the center line (y = 250) of the reduced overlayed frames
center_line = renomalized_gaussian_blob(center_y, :);

% Repeat this line across the y-axis for 200 times to create new frames
repeated_center_line_frames = repmat(center_line, height, 1, VSinfo.numFrames);

% % Save the new frames as a GIF
% filename = 'reduced_center_line_repeated_animation.gif';
% for frame = 1:frames
%     im = uint8(squeeze(repeated_center_line_frames(:, :, frame)));
%     [A, map] = gray2ind(im, 256);
%     if frame == 1
%         imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 1/frame_rate);
%     else
%         imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/frame_rate);
%     end
% end
%%
% twoDStim = 0;
% if twoDStim
%     gaussian_stim = renomalized_gaussian_blob;
% else
%     gaussian_stim = repeated_center_line_frames;
% end

%%
% weigh the two
reduced_overlayed_frames = min(normalized_frames .* wBack + repeated_center_line_frames .* wGauss,255);
upperBandEdge = ScreenInfo.yaxis - ScreenInfo.liftingYaxis - center_y ;
lowerBandEdge = ScreenInfo.yaxis - ScreenInfo.liftingYaxis + center_y -1;
fullImage = zeros(ScreenInfo.yaxis,ScreenInfo.xaxis,VSinfo.numFrames);
fullImage(upperBandEdge:lowerBandEdge,:,:) = reduced_overlayed_frames;
for frame = 1:VSinfo.numFrames
    stimTextureMatrix(frame) =  Screen('MakeTexture', windowPtr, fullImage(:,:,frame));
end
% Save the overlayed frames as a GIF
% filename = 'reduced_overlayed_normalized_animation.gif';
% for frame = 1:frames
%     im = uint8(squeeze(reduced_overlayed_frames(:, :, frame)));
%     [A, map] = gray2ind(im, 256);
%     if frame == 1
%         imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 1/frame_rate);
%     else
%         imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/frame_rate);
%     end
% end


end
