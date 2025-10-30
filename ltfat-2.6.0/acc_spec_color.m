% Parameters
a = 8;
M = 40;
L = a*M;

win_types = {'gaussian','box'};
rgbs = cell(1,2);
titlestr = cell(1,2);

for t = 1:2
    % --- Window ---
    wt = win_types{t};
    switch wt
        case 'gaussian'
            g = pgauss(L);
        case 'box'
            g = zeros(L,1);
            w = floor(L/20);
            g(1:w) = 1;
            g(end-w+1:end) = 1;
        otherwise
            error('Unsupported window type.');
    end

    % --- Frames and symbol ---
    [Fa, Fs] = framepair('dgt', g, 'dual', a, M);
    symbol = load_symbol(0, M);
    symbol = double(symbol > 0.5);
    s = framenative2coef(Fa, symbol);

    % --- Eigenpairs of multiplier ---
    [V, D] = framemuleigs(Fa, Fs, s, a*M);
    lam = real(D); if ~isvector(lam), lam = real(diag(lam)); end

    % --- Pick 3 closest to 0.5 ---
    target = 0.5;
    [~, idx_sorted] = sort(abs(lam - target));
    sel = idx_sorted(1:3);

    % --- Build RGB image: one eigvec per channel ---
    rgb = zeros(M, M, 3);
    for k = 1:3
        v = V(:, sel(k));
        S = abs(dgt(v, g, a, M)).^2;      % M x (L/a)
        den = max(S(:)); if den > 0, S = S/den; end
        rgb(:,:,k) = S;
    end

    rgbs{t} = rgb;
    titlestr{t} = sprintf('%s window', wt);
end

% --- Plot side-by-side ---
figure();
subplot(1,2,1);
imshow(rgbs{1}, 'InitialMagnification','fit'); axis on;
xlabel('Frequency bin');
ylabel('Time frame');
title(titlestr{1});

subplot(1,2,2);
imshow(rgbs{2}, 'InitialMagnification','fit'); axis on;
xlabel('Frequency bin');
ylabel('Time frame');
title(titlestr{2});

sgtitle(sprintf('Spectrograms in RGB channels: 3 eigvectors closest to %.2f', target));