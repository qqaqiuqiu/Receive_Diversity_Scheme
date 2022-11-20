clear all; clc;
N = 10^6;                           % sample signal number
SNRdB = 0:2:20;                     % SNR = [0, 2, 4, 6, 8, 10]
SNR = 10 .^ (0.1 .* SNRdB);         % SNR in linear scale
x = rand(1,N) > 0.5;                % sample signal (0 or 1)
L = 4;                              % Branch number

rho_th_dB = 6;
rho_th = 10 .^ (0.1 * rho_th_dB);

% 1Rx ---------------------------------------------------------------------
for i = 1:length(SNRdB)

    h = 1/sqrt(2) * [randn(1, N) + j*randn(1, N)];    % ~CN(0, 1)
    rho = SNR(i) * power(abs(h), 2);    % SNR in receiver

    hHat = rho < rho_th;         % detection
     
    outage(i) = size(find(hHat),2);     % record the number of error signal
end
outprob = outage / N;  % outage probability

% SC ----------------------------------------------------------------------
for i_sc = 1 : length(SNR)

    for j_sc = 1 : L
        h_sc(:, :, j_sc) = 1/sqrt(2) * [randn(1, N) + j*randn(1, N)];
    end

    h_sc_opt = zeros(1, N);
    for j_sc = 1 : N
        for k_sc = 1 : L
            if abs(h_sc(1, j_sc, k_sc)) > abs(h_sc_opt(1, j_sc))
                h_sc_opt(1, j_sc) = h_sc(1, j_sc, k_sc);
            end
        end
    end

    rho = SNR(i_sc) * power(abs(h_sc_opt), 2);
    outage_sc(i_sc) = size(find([rho < rho_th]), 2);
end
pout_sc = outage_sc / N;

% picture -----------------------------------------------------------------
figure(1);
semilogy(SNRdB, outprob, '-o');
    hold on
semilogy(SNRdB, pout_sc, '-h');

axis([0 20 10^-5 10^0])
xlabel('Eb/No');
ylabel('Outage Probability');
title('Comparison of Different rho threshold for SSC');

% SSC in different rho ----------------------------------------------------
for i = 0:5:20
    out_ssc(i);
end

% legned ------------------------------------------------------------------
L=legend('1Rx', '4Rx SC', 'SSC with threshold 0 dB', 'SSC with threshold 5 dB', ...
    'SSC with threshold 10 dB', 'SSC with threshold 15 dB', 'SSC with threshold 20 dB');
set(L,'Fontsize',12);