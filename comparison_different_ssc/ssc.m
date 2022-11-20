function ssc(rho_threshold)

N = 10^6;                           % sample signal number
SNRdB = 0:2:20;                     % SNR = [0, 2, 4, 6, 8, 10]
SNR = 10 .^ (0.1 .* SNRdB);         % SNR in linear scale
x = rand(1,N) > 0.5;                % sample signal (0 or 1)
L = 4;                              % Branch number

% Simulation --------------------------------------------------------------
s = [2*x - 1];

% SSC ---------------------------------------------------------------------
for i_ssc = 1 : length(SNR)

    deviate_ssc = sqrt(0.5 / SNR(i_ssc));

    for j_ssc = 1 : L
        n_ssc(:, :, j_ssc) = [randn(1, N) + j*randn(1, N)];
        h_ssc(:, :, j_ssc) = 1/sqrt(2) * [randn(1, N) + j*randn(1, N)];
    end

    h_ssc_opt = zeros(1, N);
    n_ssc_opt = zeros(1, N);
    rho_th = rho_threshold;
    rho_th_linear = 10 ^ (0.1 * rho_th);
    tmp = 0;

    for j_ssc = 1 : N

        if  j_ssc ~= 1

            if (SNR(i_ssc) * power(abs(h_ssc_opt(1, (j_ssc - 1))), 2)) >= rho_th_linear
                h_ssc_opt(1, j_ssc) = h_ssc(1, j_ssc, tmp);
                n_ssc_opt(1, j_ssc) = n_ssc(1, j_ssc, tmp);
                continue
            end
        end
            
        for k_ssc = 1 : L
            if abs(h_ssc(1, j_ssc, k_ssc)) > abs(h_ssc_opt(1, j_ssc))
                h_ssc_opt(1, j_ssc) = h_ssc(1, j_ssc, k_ssc);
                n_ssc_opt(1, j_ssc) = n_ssc(1, j_ssc, k_ssc);
                tmp = k_ssc;
            end
        end
    end

    y_ssc = h_ssc_opt .* s + deviate_ssc * n_ssc_opt;
    r_ssc = conj(h_ssc_opt) .* y_ssc;

    xHat_ssc = real(r_ssc) > 0;
    error_ssc(i_ssc) = size(find([x - xHat_ssc]),2);
end
errorbitrate_ssc = error_ssc / N;

% picture------------------------------------------------------------------
figure(1);
semilogy(SNRdB, errorbitrate_ssc, '-o');

axis([0 20 10^-5 10^0])
grid on