% By: Rocio Lopez
% 03/13/2025
clc; clear; close all;
fsamp = 20000;
fcuts = [4000 4500];
mags = [1 0]; % desired magnitude
dev = [0.1 0.05]; % ripple
[M, w, beta, ftype] = kaiserord(fcuts, mags, dev, fsamp);
M = M + rem(M,2);
hh = fir1(M, w, ftype, kaiser(M+1, beta), "noscale");
% Plot
[H, w] = freqz(hh, 1, 1024, fsamp);
plot(w, abs(H))
grid on
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Kaiser Window FIR Filter');

[M, fo, ao, w] = firpmord(fcuts, mags, dev, fsamp);
M = M + rem(M,2);
ho = firpm(M, fo, ao, w);
[Ho, w] = freqz(ho,1,1024,fsamp); 
n = 0:length(ho)-1;
stem(n, ho, 'filled');
xlabel('n');ylabel('h_0(n)');
title('Impulse Response of Equiripple FIR Filter');
freqz(ho,1,1024,fsamp); 
title('Frequency Response of Equiripple FIR Filter ');
zplane(ho, 1);
title('Zeros of H_0(z)');
deltaP = max(abs(Ho(w <= fcuts(1)))) - min(abs(Ho(w <= fcuts(1))));
deltaS = max(abs(Ho(w >= fcuts(2))));

M_half = ceil(length(ho) / 2); 
h1 = ho;
h1(M_half) = ho(M_half) + deltaS + 1e-4;
freqz(h1,1,1024,fsamp);
title('Lift H_0(e^{j\omega}) by the Stopband Ripple');
stem(n, real(h1), 'filled');
xlabel('n');ylabel('h_1(n)');
title('Impulse Response of h_1(n)');
zplane(h1, 1);
title('Zeros of H_1(z) After Increasing the Center Sample by {\delta}_s')

tol =0;
h1_zeros = roots(h1);
h1_zeros_filtered = h1_zeros(abs(h1_zeros) <= 1 + tol);
zplane(h1_zeros_filtered, 1); 
title('Removed All Zeros Outside TUC');

h2 = real(poly(h1_zeros_filtered));
n = 0:length(h2)-1;
stem(n, h2, 'filled');
xlabel('n');ylabel('h_0(n)');
title('Impulse Response of Minimum-Phase Filter.');
[H2, w] = freqz(h2,1,1024,fsamp);
freqz(h2,1,1024,fsamp);
title('Frequency Response of H_2(e^{j\omega})');
d_max = max(abs(H2(w <= fcuts(1))));
d_min = min(abs(H2(w <= fcuts(1))));
d_center = (d_max + d_min) / 2;
H2_n = (H2 /d_center);
plot(w, abs(Ho));
hold on;
plot(w, abs(H2_n), 'r--');  
legend('H(e^{j\omega})', 'H_{min}(e^{j\omega})');
title('Magnitude Response of Designed Filter and Normalized Minimum-Phase Filter.');
