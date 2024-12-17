%addpath(genpath('ltfat-2.6.0'));
%ltfatstart;
a = 10;
M = 50;
g = pgauss(a*M);

[Fa, Fs] = framepair('dgt', g, 'dual', a, M);

symbol_small = load_symbol(0, M);
symbol_big = load_symbol(0, M);

s_small = framenative2coef(Fa, symbol_small);
s_big = framenative2coef(Fa, symbol_big);

acc_spec_small = zeros(M, M);
acc_spec_big = zeros(M, M);

acc_spec_small_all = zeros(M, M, a*M);
acc_spec_big_all = zeros(M, M, a*M);

[V_small, D_small] = framemuleigs_mod(Fa, Fs, s_small, a*M);
[V_big, D_big] = framemuleigs_mod(Fa, Fs, s_big, a*M);

for k = 1:a*M
    acc_spec_small = acc_spec_small + real(D_small(k)) * abs(dgt(V_small(:, k), g, a, M)).^2;
    acc_spec_big = acc_spec_big + real(D_big(k)) * abs(dgt(V_big(:, k), g, a, M)).^2;

    acc_spec_small_all(:, :, k) = acc_spec_small;
    acc_spec_big_all(:, :, k) = acc_spec_big;
end

close all;

fig = figure;
fig.Position = [100, 1000, 1800, 200];

for k = 1:10
    index = k*5;
    subplot(2, 10, k);
    imagesc(acc_spec_small_all(:, :, index));
    colormap(flipud(gray));
    colorbar;
    title("small " + int2str(index));

    subplot(2, 10, k+10);
    imagesc(acc_spec_big_all(:, :, index));
    colormap(flipud(gray));
    colorbar;
    title("big " + int2str(index));
end

%exportgraphics(gcf,'figures/white_noise_illustration.png','Resolution',300) 
