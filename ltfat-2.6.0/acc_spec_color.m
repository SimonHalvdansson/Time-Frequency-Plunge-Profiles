% -------------------- Initialize and Set Parameters --------------------
% Add LTFAT library to the MATLAB path (uncomment if LTFAT is installed)
% Ensure you have LTFAT installed and update the path accordingly
% addpath(genpath('ltfat-2.6.0'));
% ltfatstart;

% Parameters for frame and symbol
a = 5;            % Time-shift parameter
M = 20;           % Number of frequency channels
L = a * M;        % Frame length

% -------------------- Select Window Type --------------------
% Choose between 'Gaussian' and 'Box'
windowType = 'Gaussian'; % Change to 'Box' to use a box window

% -------------------- Generate Window --------------------
switch lower(windowType)
    case 'gaussian'
        g = pgauss(L); % Gaussian window from LTFAT
    case 'box'
        % Generate a box window with two rectangular regions
        g = zeros(L, 1);
        w = floor(L / 20); % Width of the box edges
        
        % Ensure indexing is correct
        g(1:w) = 1;
        g(end-w+1:end) = 1;
        
        % Normalize the window
        g = g / sqrt(sum(g.^2));
    otherwise
        error('Unsupported window type. Choose either "Gaussian" or "Box".');
end

% -------------------- Create Dual Frame Pair for DGT --------------------
% Create dual frames for the Discrete Gabor Transform (DGT)
[Fa, Fs] = framepair('dgt', g, 'dual', a, M);

% -------------------- Load the Binary Symbol Image --------------------
% Load a binary symbol image
% Ensure that you have a function 'load_symbol' that returns an MxM binary image
% Replace '0' with the desired symbol index or modify as needed
symbol = load_symbol(0, M);  
symbol = double(symbol > 0.5); % Convert to binary (0 and 1)

% Convert the symbol image to frame coefficients
s = framenative2coef(Fa, symbol);

% -------------------- Compute Eigenvalues --------------------
% Compute eigenvalues and eigenvectors for the frame
[V, D] = framemuleigs(Fa, Fs, s, a*M);
D = real(D); % Ensure eigenvalues are real

% -------------------- Select the Three Eigenvectors Closest to 0.5 --------------------
% Define the target eigenvalue
target_eigenvalue = 0.5;

% Compute the absolute difference from the target
diff_from_target = abs(D - target_eigenvalue);

% Find the indices of the three smallest differences
[~, sorted_indices] = sort(diff_from_target);

% Select the top three eigenvalues and their corresponding indices
num_selected = 3;
if length(D) < num_selected
    warning('Number of available eigenvalues (%d) is less than the desired number of selected eigenvectors (%d). Selecting all available eigenvectors.', length(D), num_selected);
    num_selected = length(D);
end

selected_indices = sorted_indices(1:num_selected);
selected_eigenvalues = D(selected_indices);

% -------------------- Sum Spectrograms for Selected Eigenvalues with Colors --------------------
% Initialize sum of spectrograms with three channels for RGB
sum_spectrogram = zeros(M, M, 3);

% Assign random colors to each selected eigenvalue
% Colors are scaled to avoid being too dark (values between 0.2 and 1)
colors = rand(num_selected, 3) * 0.8 + 0.2; 
color_idx = 1; % Initialize color index

% Iterate through each selected eigenvalue and accumulate spectrograms with assigned colors
for k = 1:num_selected
    idx = selected_indices(k);
    
    % Compute Discrete Gabor Transform (DGT) for the eigenvector
    spectrogram = abs(dgt(V(:, idx), g, a, M)).^2;
    
    % Normalize the spectrogram to have values between 0 and 1
    spectrogram = spectrogram / max(spectrogram(:));
    
    % Accumulate the spectrogram into each RGB channel with its unique color
    sum_spectrogram(:, :, 1) = sum_spectrogram(:, :, 1) + colors(color_idx, 1) * spectrogram;
    sum_spectrogram(:, :, 2) = sum_spectrogram(:, :, 2) + colors(color_idx, 2) * spectrogram;
    sum_spectrogram(:, :, 3) = sum_spectrogram(:, :, 3) + colors(color_idx, 3) * spectrogram;
    
    color_idx = color_idx + 1; % Move to the next color
end

% Normalize each RGB channel to the range [0, 1]
for c = 1:3
    channel = sum_spectrogram(:, :, c);
    max_val = max(channel(:));
    if max_val > 0
        sum_spectrogram(:, :, c) = channel / max_val;
    end
end

% -------------------- Display the Summed RGB Spectrogram --------------------
figure;
imshow(sum_spectrogram, 'InitialMagnification', 'fit');
axis on;
xlabel('Frequency Bins');
ylabel('Time Frames');

% Create a descriptive title
title_str = sprintf('Accumulated Spectrogram: %d Eigenvectors Closest to %.2f [%s]', num_selected, target_eigenvalue, windowType);
title(title_str);
