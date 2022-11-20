function out_ssc(rho_ssc_dB)

N = 10^6;                           % sample signal number
SNRdB = 0:2:20;                     % SNR = [0, 2, 4, 6, 8, 10]
SNR = 10 .^ (0.1 .* SNRdB);         % SNR in linear scale
L = 4;                              % Branch number

rho_th_dB = 6;
rho_th = 10 .^ (0.1 * rho_th_dB);

for i_ssc = 1 : length(SNR)

    for j_ssc = 1 : L
        h_ssc(:, :, j_ssc) = 1/sqrt(2) * [randn(1, N) + j*randn(1, N)];
    end

    h_ssc_opt = zeros(1, N);
    rho_ssc_dB = rho_ssc_dB;
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

% picture------------------------------------------------------------------
figure(1);
semilogy(SNRdB, pout_ssc, '-o');

axis([0 20 10^-5 10^0])
grid on