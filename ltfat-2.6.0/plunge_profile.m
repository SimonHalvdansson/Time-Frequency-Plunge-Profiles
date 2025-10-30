addpath(genpath('ltfat-2.6.0'));
ltfatstart;

% Parameters for frame and symbol
a = 10;
M = 40;
L = a * M; % Total window length

% Choose between 'Gaussian' and 'Box'
windowType = 'Gaussian';

switch lower(windowType)
    case 'gaussian'
        g = pgauss(L);
    case 'box'
        g = zeros(L,1);
        w = floor(L/20);
        g(1:w) = 1;
        g(end-w+1:end) = 1;
    otherwise
        error('Unsupported window type. Choose either "Gaussian" or "Box".');
end

[Fa, Fs] = framepair('dgt', g, 'dual', a, M);

%load symbol
symbol = load_symbol(0, M);
symbol = double(symbol > 0.5);
BW = logical(symbol);

s = framenative2coef(Fa, symbol);

%eigenvalues
D = framemuleigs(Fa, Fs, s, a * M);
D = real(D);

props = regionprops(BW_clean, 'Perimeter');
totalPerimeter_regionprops = sum([props.Perimeter]);

%compute erfc
k = 1:L;
A = sum(s) * a / M; % Area normalized
B = totalPerimeter_regionprops * sqrt(a / M); % Perimeter normalized

erfc_values = 0.5 * erfc((k - A) / (B / sqrt(2 * pi)));

% -------------------- Count Number of Eigenvalues in Range --------------------
num_eigenvalues = sum(D > 0.1 & D < 0.9);
fprintf('Number of eigenvalues in the range 0.1 < lambda_k < 0.9: %d\n', num_eigenvalues);

% Define the plunge region as where erfc_values are between 0.1 and 0.9
plunge_region = (erfc_values > 0.1) & (erfc_values < 0.9);
num_plunge_points = sum(plunge_region);

B_div_num_plunge_points = B / num_plunge_points;
fprintf('B divided by the number of points in the plunge region: %.4f pixels per point\n', B_div_num_plunge_points);

%setup plot
close all;
figure('Name', 'Combined Plots', 'NumberTitle', 'off', 'Position', [100, 100, 1350, 400]);
t = tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

nexttile(1, [2 1]);

% left symbol subfigure
imshow(BW_clean);
hold on;

[B_boundaries, ~] = bwboundaries(BW_clean);

% Define colors for different boundaries
colors = lines(length(B_boundaries));

% Plot each boundary with a unique color
for idx = 1:length(B_boundaries)
    boundary = B_boundaries{idx};
    plot(boundary(:,2), boundary(:,1), 'LineWidth', 2, 'Color', colors(idx,:));
end

hold off;
title('Symbol and boundary', 'FontSize', 14);

%info box
infoText = sprintf('Parameters:\n a = %d\n M = %d\n Window: %s', a, M, windowType);
annotation('textbox', [0.12, 0.05, 0.9, 0.15], 'String', infoText, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontSize', 12, 'HorizontalAlignment', 'left');

%top right subplot
nexttile(2, [1 2]);

hD = plot(k, D, 'b-', 'LineWidth', 3, 'DisplayName', '$\lambda_k^\Omega$');
hold on;

erfc_display = '$\frac{1}{2} \mathrm{erfc}\left(\sqrt{2\pi} \frac{k - |\Omega|}{|\partial \Omega|}\right)$';
hE = plot(k, erfc_values, 'r--', 'LineWidth', 3, 'DisplayName', erfc_display);

title('Eigenvalues and ideal erfc', 'FontSize', 14);
xlabel('Index (k)', 'FontSize', 12);
ylabel('Value', 'FontSize', 12);
grid on;

%analyze plunge
k_low = find(plunge_region, 1, 'first');
k_high = find(plunge_region, 1, 'last');

% shaded plunge band
yl = ylim;
h = patch([k_low k_high k_high k_low], ...
          [yl(1) yl(1) yl(2) yl(2)], ...
          [0.3 0.3 0.3], ...        % darker gray
          'FaceAlpha', 0.25, ...
          'EdgeColor', 'none', ...
          'DisplayName', 'Plunge region');
uistack(h,'bottom');                 % keep curves on top

annotation_text = sprintf('Perimeter = %.2f \nEigenvalues in plunge region: %d\nPerimeter / Plunge Points = %.4f', B, num_eigenvalues, B_div_num_plunge_points);

% annotation box position set here
annotation('textbox', [0.55, 0.75, 0.2, 0.15], 'String', annotation_text, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontSize', 10, 'Interpreter', 'none');

legend([hD hE h], 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
hold off;

% bottom right subplot
nexttile(5, [1 2]);

% Compute the difference between eigenvalues and erfc values
error_values = D - erfc_values.';

% Plot the error with reduced line width
plot(k, error_values, 'k-', 'LineWidth', 2, 'DisplayName', 'Error');

title('Erfc discrepancy', 'FontSize', 14);
xlabel('Index (k)', 'FontSize', 12);
ylabel('Discrepancy', 'FontSize', 12);
grid on;

fprintf('B (rescaled perimeter): %.2f pixels\n', B);