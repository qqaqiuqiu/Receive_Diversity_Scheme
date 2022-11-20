clear all; clc;
N = 10^6;                           % sample signal number
SNRdB = 0:2:10;                     % SNR = [0, 2, 4, 6, 8, 10]
SNR = 10 .^ (0.1 .* SNRdB);         % SNR in linear scale
x = rand(1,N) > 0.5;                % sample signal (0 or 1)
L = 4;                              % Branch number

rho_th_dB = 6;
rho_th = 10 .^ (0.1 * rho_th_dB);

% Simulation --------------------------------------------------------------
s = [2*x - 1];

% MRC ---------------------------------------------------------------------
for i_mrc = 1 : length(SNR)

    h_sigma_mrc = 0;
    for j_mrc = 1 : L
        h_sigma_mrc = h_sigma_mrc + power(abs(1/sqrt(2) * [randn(1, N) + j*randn(1, N)]), 2);
    end 

    rho = SNR(i_mrc) * h_sigma_mrc;
    outage_mrc(i_mrc) = size(find([rho < rho_th]), 2);
end
pout_mrc = outage_mrc / N;

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

% SSC ---------------------------------------------------------------------
for i_ssc = 1 : length(SNR)

    for j_ssc = 1 : L
        h_ssc(:, :, j_ssc) = 1/sqrt(2) * [randn(1, N) + j*randn(1, N)];
    end

    h_ssc_opt = zeros(1, N);
    rho_ssc_dB = 6;
    rho_ssc = 10 ^ (0.1 * rho_ssc_dB);
    tmp = 0;

    for j_ssc = 1 : N
        if  j_ssc ~= 1
            if (SNR(i_ssc) * power(abs(h_ssc_opt(1, (j_ssc - 1))), 2)) >= rho_ssc
                h_ssc_opt(1, j_ssc) = h_ssc(1, j_ssc, tmp);
                continue
            end
        end
            
        for k_ssc = 1 : L
            if abs(h_ssc(1, j_ssc, k_ssc)) > abs(h_ssc_opt(1, j_ssc))
                h_ssc_opt(1, j_ssc) = h_ssc(1, j_ssc, k_ssc);
                tmp = k_ssc;
            end
        end
    end

    rho = SNR(i_ssc) * power(abs(h_ssc_opt), 2);
    outage_ssc(i_ssc) = size(find([rho < rho_th]), 2);
end
pout_ssc = outage_ssc / N;

% EGC ---------------------------------------------------------------------
for i_egc = 1 : length(SNR)

    h_sigma_egc = 0;
    for j_egc = 1 : L
        h_sigma_egc = h_sigma_egc + abs(1/sqrt(2) * [randn(1, N) + j*randn(1, N)]);
    end 

    rho = SNR(i_egc) * power(h_sigma_egc, 2) / L;
    outage_egc(i_egc) = size(find([rho < rho_th]), 2);
end
pout_egc = outage_egc / N;

% theory MRC --------------------------------------------------------------
sigma = 0;
for k = 0 : L-1
    parta = 1 / factorial(k);
    partb = power((rho_th ./ SNR), k);
    sigma = sigma + parta * partb;
end
pout_mrc_theory = 1 - exp(-rho_th ./ SNR) .* sigma;

% theory SC ---------------------------------------------------------------
pout_sc_theory = power((1 - exp(-rho_th ./ SNR)), L);

% picture------------------------------------------------------------------
figure(1);
semilogy(SNRdB, pout_mrc, '-o');
  hold on
semilogy(SNRdB, pout_sc, '-h');
semilogy(SNRdB, pout_ssc,'-x');
semilogy(SNRdB, pout_egc,'->');
semilogy(SNRdB, pout_mrc_theory);
semilogy(SNRdB, pout_sc_theory);

axis([0 10 10^-4 10^0])
grid on
L=legend('MRC', 'SC', 'SSC with threshold 6dB', 'EGC', 'MRC Theory', ...
    'SC/SSC Theory');
set(L,'Fontsize',12);

xlabel('Eb/No');
ylabel('Outage Probability');
title('Receive diversity Comparison for BPSK Over Rayleigh fading');