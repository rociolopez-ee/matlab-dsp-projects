clc;
clear all;
close all;

N= 50;
M = N/2;
Fs = 24000;
w = linspace(0, pi, 1024);

freqs = w * Fs / (2*pi); % Hz
D_w = zeros(size(w));
D_w(freqs <= 4000) = 1;
D_w((freqs >= 7500) & (freqs <= 8500)) = linspace(0, 1, sum((freqs >= 7500) & (freqs <= 8500))); 

d = trapz(w, D_w.^2);

Q_size = zeros(M+1, M+1, length(w));
for i = 1:length(w)
    c = cos((0:M)' * w(i));
    Q_size(:,:,i) = c * c.';
end
Q = trapz(w, Q_size, 3);

p_size = zeros(M+1, length(w));
for i = 1:length(w)
    c = cos((0:M)' * w(i));
    p_size(:,i) = D_w(i) * c;
end
p = trapz(w, p_size, 2);

% Q_t
Q_t = [Q,     p;
      p.',   d];

[V,D] = eig(Q_t);
    % D: is the eigenvalues
    % V: eigenvectors
eigenvalues = diag(D);
[~, idx] = min(eigenvalues);
a_hat = V(:, idx);

scale = -1 / a_hat(end);
a = a_hat * scale;

h = zeros(M+1,1);

h(M+1) = a(1);  

for k = 1:M
    h(M+1 - k) = a(k+1)/2;
    h(M+1+k) = a(k+1)/2;
end

H = freqz(h, 1, w);

plot(freqs, D_w , 'r')
hold on;
plot(freqs, abs(H), 'b')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('D(\omega)', 'A(\omega)')
hold off;

f_notches = [4500, 8000];
E = zeros(length(f_notches), M+1);
for k = 1:length(f_notches)
    omega_k = 2 * pi * f_notches(k) / Fs;
    E(k, :) = cos((0:M) * omega_k);  % C^T(omega_k)
end
v = zeros(length(f_notches), 1); 

E_hat = [E, v];

B = null(E_hat);

Q_reduce = B' * Q_t * B;
[V2, D2] = eig(Q_reduce);
[~, idx] = min(diag(D2));
w_min = V2(:, idx);

a_hat0 = B * w_min;
a_hat0 = a_hat0 * (-1 / a_hat0(end));  % Make last value = -1
a_0 = a_hat0(1:end-1);  % Remove the -1 

h = zeros(N+1, 1);
h(M+1) = a_0(1); 
for k = 1:M
    h(M+1 - k) = a_0(k+1)/2;
    h(M+1 + k) = a_0(k+1)/2;
end

H = freqz(h, 1, w);

plot(freqs, D_w, 'r'); 
hold on;
plot(freqs, abs(H), 'b');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('D(\omega)', 'A(\omega)');
title('Eigenfilter with Notches at 4500Hz and 8000Hz');
