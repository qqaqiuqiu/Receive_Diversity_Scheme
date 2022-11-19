clear all; clc;
N = 10^6;                           % sample signal number
SNRdB = 0:2:20;                     % SNR = [0, 2, 4, 6, 8, 10]
SNR = 10 .^ (0.1 .* SNRdB);         % SNR in linear scale
x = rand(1,N) > 0.5;                % sample signal (0 or 1)
L = 4;                              % Branch number

% Simulation --------------------------------------------------------------
s = [2*x - 1];

% 1Rx ---------------------------------------------------------------------
for i = 1:length(SNR)
    n = [randn(1, N) + j*randn(1, N)];                % ~CN(0, 2)
    h = 1/sqrt(2) * [randn(1, N) + j*randn(1, N)];    % ~CN(0, 1)
    deviate = sqrt(0.5 / SNR(i));                     % 處理成高斯雜訊(包含調整到預期的 SNR)

    y = h.*s + deviate*n;           % receive signal
    r = conj(h) .* y;               % post processing
                            
    xHat = real(r) > 0;                          % determine the received signal      
    error(i) = size(find([x - xHat]),2);   % recoed the number of error signal
end
errorbitrate = error/N;

% SC ----------------------------------------------------------------------
for i_sc = 1 : length(SNR)

    deviate_sc = sqrt(0.5 / SNR(i_sc));

    for j_sc = 1 : L
        n_sc(:, :, j_sc) = [randn(1, N) + j*randn(1, N)];
        h_sc(:, :, j_sc) = 1/sqrt(2) * [randn(1, N) + j*randn(1, N)];
    end

    h_sc_opt = zeros(1, N);
    n_sc_opt = zeros(1, N);
    for j_sc = 1 : N
        for k_sc = 1 : L
            if abs(h_sc(1, j_sc, k_sc)) > abs(h_sc_opt(1, j_sc))
                h_sc_opt(1, j_sc) = h_sc(1, j_sc, k_sc);
                n_sc_opt(1, j_sc) = n_sc(1, j_sc, k_sc);
            end
        end
    end

    y_sc = h_sc_opt .* s + deviate_sc * n_sc_opt;
    r_sc = conj(h_sc_opt) .* y_sc;

    xHat_sc = real(r_sc) > 0;
    error_sc(i_sc) = size(find([x - xHat_sc]),2);
end
errorbitrate_sc = error_sc / N;

% picture -----------------------------------------------------------------
figure(1);
semilogy(SNRdB, errorbitrate, '-o');
    hold on
semilogy(SNRdB, errorbitrate_sc, '-h');

axis([0 20 10^-5 10^0])
xlabel('Eb/No');
ylabel('BitError Probability');
title('Comparison of Different rho threshold for SSC');

% SSC in different rho ----------------------------------------------------
for i = 0:5:20
    ssc(i);
end

% legned ------------------------------------------------------------------
L=legend('1Rx', '4Rx SC', 'SSC with threshold 0 dB', 'SSC with threshold 5 dB', ...
    'SSC with threshold 10 dB', 'SSC with threshold 15 dB', 'SSC with threshold 20 dB');
set(L,'Fontsize',12);