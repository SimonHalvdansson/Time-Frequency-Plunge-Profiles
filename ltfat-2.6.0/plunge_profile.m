% -------------------- Initialize and Set Parameters --------------------
addpath(genpath('ltfat-2.6.0'));
ltfatstart;

% Parameters for frame and symbol
a = 10;
M = 100;

% -------------------- Generate Gaussian Window --------------------
g = pgauss(a*M);

% -------------------- Create Dual Frame Pair for DGT --------------------
[Fa, Fs] = framepair('dgt', g, 'dual', a, M);

% -------------------- Load the Binary Symbol Image --------------------
symbol = load_symbol(8, M);
symbol = double(symbol > 0.5);

% -------------------- Convert Frame Coefficients to Native Coefficients --------------------
s = framenative2coef(Fa, symbol);

% -------------------- Compute Eigenvalues --------------------
D = framemuleigs(Fa, Fs, s, a*M);
D = real(D);

% -------------------- Ensure the 'symbol' Image is Binary --------------------
BW = logical(symbol);

% -------------------- Preprocess the Image --------------------
% Remove small objects to reduce noise (adjust the threshold as needed)
BW_clean = bwareaopen(BW, 50); % Removes objects with fewer than 50 pixels

% -------------------- Measure Perimeter Using regionprops --------------------
props = regionprops(BW_clean, 'Perimeter');

% Handle cases with multiple regions
if isempty(props)
    disp('No objects detected in the image.');
    totalPerimeter_regionprops = 0;
else
    totalPerimeter_regionprops = sum([props.Perimeter]);
    % Removed the original perimeter print statement
    % fprintf('Total Perimeter using regionprops: %.2f pixels\n', totalPerimeter_regionprops);
end

% -------------------- Define Parameters A and B --------------------
A = sum(s)*a/M; % this should be the area, properly normalized
B = totalPerimeter_regionprops*sqrt(a/M); % Set B to the total perimeter length

% -------------------- Compute erfc Function --------------------
% Define the range of k based on the length of D
k = 1:length(D);

% Compute the erfc function
if B ~= 0
    erfc_values = 0.5 * erfc((k - A) / (B/sqrt(2*pi)));
else
    erfc_values = zeros(size(k)); % Initialize with zeros if B is zero
    warning('Total Perimeter is zero. The erfc function cannot be computed.');
end

% -------------------- Close All Existing Figures --------------------
close all;

% -------------------- Combined Plot with Tiled Layout --------------------
% Create a new figure for all subplots with larger width and a slightly reduced height
figure('Name', 'Combined Plots', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 500]);

% Define a tiled layout with 2 rows and 3 columns to allow unequal column widths
t = tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% -------------------- Left Subplot: Symbol with Boundaries --------------------
% Span the big image across both rows in the first column
nexttile(1, [2 1]); % Starts at tile 1, spans 2 rows and 1 column

% Display the cleaned binary image
imshow(BW_clean);
hold on;

% Extract boundaries using bwboundaries
[B_boundaries, L] = bwboundaries(BW_clean); % 'noholes' excludes boundaries of holes within objects

% Define colors for different boundaries (optional)
colors = lines(length(B_boundaries));

% Plot each boundary with a unique color
for idx = 1:length(B_boundaries)
    boundary = B_boundaries{idx};
    plot(boundary(:,2), boundary(:,1), 'LineWidth', 2, 'Color', colors(idx,:));
end

hold off;
title('Symbol with Boundary', 'FontSize', 14); % Larger title font size

% -------------------- Top-Right Subplot: Eigenvalues and erfc --------------------
% Span the top-right plot across columns 2 and 3 in row 1
nexttile(2, [1 2]); % Starts at tile 2, spans 1 row and 2 columns

% Plot eigenvalues with modified label for legend
plot(k, D, '-', 'LineWidth', 1.5, 'DisplayName', '$\lambda_k^\Omega$');
hold on;

% Plot erfc function if B is not zero
if B ~= 0
    % Define the DisplayName with LaTeX formatting
    erfc_display = '$\frac{1}{2} \mathrm{erfc}\left(\sqrt{2\pi} \frac{k - |\Omega|}{|\partial \Omega|}\right)$';
    plot(k, erfc_values, 'r-', 'LineWidth', 2, 'DisplayName', erfc_display);
end

% Enhance the plot with titles, labels, and larger font sizes
title('Eigenvalues and Complementary Error Function', 'FontSize', 14); % Larger title font size
xlabel('Index (k)', 'FontSize', 12);
ylabel('Value', 'FontSize', 12);
grid on;

% -------------------- Define and Analyze the Plunge Region --------------------
if B ~= 0
    % Define the plunge region as where erfc_values are between 0.1 and 0.9
    plunge_region = (erfc_values > 0.1) & (erfc_values < 0.9);
    
    % Count the number of points in the plunge region
    num_plunge_points = sum(plunge_region);
    
    % Find the first and last k in the plunge region
    if num_plunge_points > 0
        k_low = find(plunge_region, 1, 'first');
        k_high = find(plunge_region, 1, 'last');
        
        % Plot vertical lines at k_low and k_high
        xline(k_low, '--b', 'LineWidth', 1.5, 'DisplayName', 'Plunge Start');
        xline(k_high, '--g', 'LineWidth', 1.5, 'DisplayName', 'Plunge End');
        
        % Removed the text annotations for "Start of Plunge" and "End of Plunge"
        % Previously, these were added using the `text` function
        
        % Optionally, you can also annotate the statistics near the plot
        % Position the text box in a corner to display B, num_eigenvalues, and fraction
        annotation_text = sprintf('Perimeter = %.2f \nEigenvalues in plunge region: %d\nPerimeter / Plunge Points = %.4f', ...
            B, num_eigenvalues, B_div_num_plunge_points);
        
        % Place the annotation box at a suitable location
        % Adjust 'Position' as needed: [x, y, width, height]
        annotation('textbox', [0.65, 0.75, 0.2, 0.15], 'String', annotation_text, ...
            'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
            'FontSize', 10, 'Interpreter', 'none');
    else
        warning('No points found in the plunge region.');
    end
end

% Add a legend with LaTeX interpreter, larger font, and modified label
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
hold off;

% -------------------- Bottom-Right Subplot: Difference (Error) --------------------
nexttile(5, [1 2]); % Starts at tile 5, spans 1 row and 2 columns

if B ~= 0
    % Compute the difference between eigenvalues and erfc values
    error_values = D - erfc_values.';
    
    % Optional: Apply smoothing to reduce oscillations
    % Uncomment the following line to apply a simple moving average filter
    % error_values = movmean(error_values, 5);
    
    % Plot the error with reduced line width
    plot(k, error_values, '-x', 'LineWidth', 1, 'Color', 'm', 'DisplayName', 'Error');
    
    % Enhance the plot with titles, labels, and larger font sizes
    title('Difference between Eigenvalues and erfc', 'FontSize', 14); % Larger title font size
    xlabel('Index (k)', 'FontSize', 12);
    ylabel('Difference', 'FontSize', 12);
    grid on;
else
    % Display a message if erfc cannot be computed
    % Adjust the position of the text to center it within the axes
    axesHandle = gca; % Get current axes
    xlim(axesHandle, [0, 1]);
    ylim(axesHandle, [0, 1]);
    text(0.5, 0.5, 'Error plot not available', 'HorizontalAlignment', 'center', 'FontSize', 12);
    title('Difference Not Available', 'FontSize', 14);
    xlabel('Index (k)', 'FontSize', 12);
    ylabel('Difference', 'FontSize', 12);
    grid on;
end

% -------------------- End of Plotting Section --------------------

% -------------------- Print Rescaled Perimeter B --------------------
fprintf('B (rescaled perimeter): %.2f pixels\n', B);

% -------------------- Count and Print Number of Eigenvalues in Range --------------------
num_eigenvalues = sum(D > 0.1 & D < 0.9);
fprintf('Number of eigenvalues in the range 0.1 < lambda_k < 0.9: %d\n', num_eigenvalues);

% -------------------- Define and Analyze the Plunge Region --------------------
if B ~= 0
    % Define the plunge region as where erfc_values are between 0.1 and 0.9
    plunge_region = (erfc_values > 0.1) & (erfc_values < 0.9);
    
    % Count the number of points in the plunge region
    num_plunge_points = sum(plunge_region);
    
    % Handle case where there are no points in the plunge region to avoid division by zero
    if num_plunge_points > 0
        B_div_num_plunge_points = B / num_plunge_points;
        fprintf('B divided by the number of points in the plunge region: %.4f pixels per point\n', B_div_num_plunge_points);
    else
        fprintf('No points found in the plunge region. Cannot compute B divided by the number of plunge points.\n');
        B_div_num_plunge_points = NaN; % Assign NaN to indicate undefined
    end
else
    fprintf('B is zero. Cannot define plunge region or compute B divided by the number of plunge points.\n');
    B_div_num_plunge_points = NaN; % Assign NaN to indicate undefined
end