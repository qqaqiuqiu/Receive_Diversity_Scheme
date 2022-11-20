clear all; clc;
N = 10^6;                           % sample signal number
SNRdB = 0:2:10;                     % SNR = [0, 2, 4, 6, 8, 10]
SNR = 10 .^ (0.1 .* SNRdB);         % SNR in linear scale
x = rand(1,N) > 0.5;                % sample signal (0 or 1)
L = 4;                              % Branch number

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

% SSC ---------------------------------------------------------------------
for i_ssc = 1 : length(SNR)

    deviate_ssc = sqrt(0.5 / SNR(i_ssc));

    for j_ssc = 1 : L
        n_ssc(:, :, j_ssc) = [randn(1, N) + j*randn(1, N)];
        h_ssc(:, :, j_ssc) = 1/sqrt(2) * [randn(1, N) + j*randn(1, N)];
    end

    h_ssc_opt = zeros(1, N);
    n_ssc_opt = zeros(1, N);
    rho_th = 10;
    rho_th_linear = 10 ^ (0.1 * rho_th);
    tmp = 0;

    %　　test -------------------------------------------------------------
    % h_ssc
    % ---------------------------------------------------------------------

    for j_ssc = 1 : N

        if  j_ssc ~= 1

            %　　test -----------------------------------------------------
            % rho = SNR(i_ssc) * abs(h_ssc_opt(1, (j_ssc - 1)))
            % rho_th_linear
            % -------------------------------------------------------------
            
            if (SNR(i_ssc) * power(abs(h_ssc_opt(1, (j_ssc - 1))), 2)) >= rho_th_linear
                h_ssc_opt(1, j_ssc) = h_ssc(1, j_ssc, tmp);
                n_ssc_opt(1, j_ssc) = n_ssc(1, j_ssc, tmp);
                
                % test ----------------------------------------------------
                % tmp
                % h_ssc_opt
                % ---------------------------------------------------------
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

        % test ------------------------------------------------------------
        % tmp
        % h_ssc_opt
        % -----------------------------------------------------------------
    end

    y_ssc = h_ssc_opt .* s + deviate_ssc * n_ssc_opt;
    r_ssc = conj(h_ssc_opt) .* y_ssc;

    xHat_ssc = real(r_ssc) > 0;
    error_ssc(i_ssc) = size(find([x - xHat_ssc]),2);
end
errorbitrate_ssc = error_ssc / N;

% EGC ---------------------------------------------------------------------
for i_egc = 1 : length(SNR)

    deviate_egc = sqrt(0.5 / SNR(i_egc));

    for j_egc = 1 : L
        n_egc(:, :, j_egc) = [randn(1, N) + j*randn(1, N)];
        h_egc(:, :, j_egc) = 1/sqrt(2) * [randn(1, N) + j*randn(1, N)];
    end 

    y_egc = h_egc .* s + deviate_egc * n_egc;
    post_processing = conj(h_egc) ./ abs(h_egc);
    r_egc = post_processing .* y_egc;

    sigma = zeros(1, N);
    for j_egc = 1 : L
        sigma = sigma + r_egc(:, :, j_egc);
    end
    xHat_egc = real(sigma) > 0;
    error_egc(i_egc) = size(find([x - xHat_egc]),2);
end
errorbitrate_egc = error_egc / N;

% picture -----------------------------------------------------------------
figure(1);
semilogy(SNRdB, errorbitrate_mrc, '-o');
    hold on
semilogy(SNRdB, errorbitrate_sc, '-h');
semilogy(SNRdB, errorbitrate_ssc, '-x');
semilogy(SNRdB, errorbitrate_egc, '->');

axis([0 10 10^-5 10^0])
L=legend('MRC', 'SC', 'SSC with threshold 10 dB', 'EGC');
set(L,'Fontsize',12);
grid on
xlabel('Eb/No');
ylabel('BitError Probability');
title('Receive diversity Comparison for BPSK Over Rayleigh fading');