%% Gaussian Atlas
% Parameters
frames = 30;
frame_rate = 30;
increment = 0.2;
image_size_x = 1024;
image_size_y = 200;
back_std_x = 60;
back_std_y = back_std_x;
stim_std_x = 60;
stim_std_y = 20; %doesn't matter when repeating lines
wBack = 0.5;
wGauss = 0.5;

% Generate the anisotropic Gaussian blob
[x, y] = meshgrid(0:image_size_x-1, 0:image_size_y-1);
center_x = image_size_x / 2;
center_y = image_size_y / 2;
mean_luminance = 1.0;

back_gaussian = mean_luminance * exp(-((x - center_x).^2 / (2 * back_std_x^2) + (y - center_y).^2 / (2 * back_std_y^2)));

% Perform Fourier transformation on the anisotropic Gaussian blob
background_gaussian_fft = fftshift(fft2(back_gaussian));

% Extract the amplitude and phase
amplitude = abs(background_gaussian_fft);
initial_phase = rand(size(amplitude)) * 2 * pi;

% Create random + or - 0.2 increment directions
phase_directions = randi([0, 1], size(initial_phase)) * 2 - 1;
phase_directions = phase_directions * increment;

% Initialize the 3D matrix to store all frames
all_frames = zeros(image_size_y, image_size_x, frames);

% Generate all frames with updated phases
for frame = 1:frames
    initial_phase = initial_phase + phase_directions;
    new_fft = amplitude .* exp(1i * initial_phase);
    new_image_fft_shifted = ifft2(ifftshift(new_fft));
    all_frames(:, :, frame) = abs(new_image_fft_shifted);
end

% Normalize the entire 3D time-space matrix
max_val = max(all_frames(:));
normalized_frames = (all_frames / max_val) * 255;

% Renormalize the Gaussian blob
stim_gaussian = mean_luminance * exp(-((x - center_x).^2 / (2 * stim_std_x^2) + (y - center_y).^2 / (2 * stim_std_y^2)));
renomalized_gaussian_blob = (stim_gaussian / max(stim_gaussian(:))) * 255;
%%
% Pick out the center line (y = 250) of the reduced overlayed frames
center_line = renomalized_gaussian_blob(center_y, :);

% Repeat this line across the y-axis for 200 times to create new frames
repeated_center_line_frames = repmat(center_line, image_size_y, 1, frames);

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
twoDStim = 0;
if twoDStim
    gaussian_stim = renomalized_gaussian_blob;
else
    gaussian_stim = repeated_center_line_frames;
end

%%
% weigh the two
reduced_overlayed_frames = min(normalized_frames .* wBack + gaussian_stim .* wGauss,255);

% Save the overlayed frames as a GIF
filename = 'reduced_overlayed_normalized_animation.gif';
for frame = 1:frames
    im = uint8(squeeze(reduced_overlayed_frames(:, :, frame)));
    [A, map] = gray2ind(im, 256);
    if frame == 1
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 1/frame_rate);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/frame_rate);
    end
end

