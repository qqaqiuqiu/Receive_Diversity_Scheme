clear all; clc;

SNRdB = 0:2:10;                     % SNR = [0, 2, 4, 6, 8, 10]
SNR = 10 .^ (0.1 .* SNRdB);         % SNR in linear scale
L = 4;                              % Branch number
N = 10^6;                           % sample signal number

x = rand(1,N) > 0.5;                % sample signal (0 or 1)

% Simulation --------------------------------------------------------------
s = [2*x - 1];

% MRC ---------------------------------------------------------------------
for i_mrc = 1 : length(SNR)

    deviate_mrc = sqrt(0.5 / SNR(i_mrc));

    for j_mrc = 1 : L
        n_mrc(:, :, j_mrc) = [randn(1, N) + j*randn(1, N)];
        h_mrc(:, :, j_mrc) = 1/sqrt(2) * [randn(1, N) + j*randn(1, N)];
    end 

    y_mrc = h_mrc .* s + deviate_mrc * n_mrc;
    r_mrc = conj(h_mrc) .* y_mrc;

    sigma = zeros(1, N);
    for j_mrc = 1 : L
        sigma = sigma + r_mrc(:, :, j_mrc);
    end
    xHat_mrc = real(sigma) > 0;
    error_mrc(i_mrc) = size(find([x - xHat_mrc]),2);
end
errorbitrate_mrc = error_mrc / N;

% theory ------------------------------------------------------------------
sigma = 0;

% teacher's method, but not practical
% for k = 0 : (L-1)
%     parta = power((1 ./ (SNR + 1)), k);
%     partb = factorial(2 * k) / (power(4, k) * power(factorial(k), 2));
%     sigma = sigma + parta * partb;
% end
% theory = 0.5 - 0.5 * (SNR ./ (SNR + 1)) .* sigma;

p = 0.5 - 0.5 * power((1 + (1 ./ SNR)), -0.5);
for k = 0: (L-1)
    parta = factorial(L-1+k) / (factorial(k) * factorial(L-1));
    partb = power(1-p, k);
    sigma = sigma + parta * partb;
end
theory = power(p, L) .* sigma;

% Upperbound --------------------------------------------------------------
upperbound = 0.5 * power(SNR, -L);

% Picture -----------------------------------------------------------------
figure(1);
semilogy(SNRdB, theory, '-o');
    hold on
semilogy(SNRdB, upperbound, '-h');
semilogy(SNRdB, errorbitrate_mrc, '-x');

axis([0 10 10^-5 10^0])
L=legend('theory', 'upperbound', 'Simulation');
set(L,'Fontsize',12);
grid on
xlabel('Eb/No');
ylabel('BitError Probability');
title('theoretical BER and upperbound of MRC for BPSK Over Rayleigh fading');


