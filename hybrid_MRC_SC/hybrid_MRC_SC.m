function hw1_4c(branch_number)

N = 10^5;                           % sample signal number
SNRdB = 0:2:10;                     % SNR = [0, 2, 4, 6, 8, 10]
SNR = 10 .^ (0.1 .* SNRdB);         % SNR in linear scale
x = rand(1,N) > 0.5;                % sample signal (0 or 1)
L = branch_number;                  % Branch number

% Simulation --------------------------------------------------------------
s = [2*x - 1];

for i_hy_mrc = 1 : length(SNR)

    deviate_hy_mrc = sqrt(0.5 / SNR(i_hy_mrc));

    for j_mrc = 1 : L
        n_mrc(:, :, j_mrc) = [randn(1, N) + j*randn(1, N)];
        h_mrc(:, :, j_mrc) = 1/sqrt(2) * [randn(1, N) + j*randn(1, N)];
    end 

    
    % compare the channel (need to be make first)
    for k_hy_mrc = 1 : N

        % choose the largest channel
        for j_h1_mrc = 1 : L
            if j_h1_mrc == 1
                h_hy_mrc(:, k_hy_mrc, 1) = h_mrc(:, k_hy_mrc, j_h1_mrc);
                n_hy_mrc(:, k_hy_mrc, 1) = n_mrc(:, k_hy_mrc, j_h1_mrc);
                h_hy_mrc(:, k_hy_mrc, 3) = j_h1_mrc;
                continue
            end

            if abs(h_mrc(:, k_hy_mrc, j_h1_mrc)) > abs(h_hy_mrc(:, k_hy_mrc, 1))
                h_hy_mrc(:, k_hy_mrc, 1) = h_mrc(:, k_hy_mrc, j_h1_mrc);
                n_hy_mrc(:, k_hy_mrc, 1) = n_mrc(:, k_hy_mrc, j_h1_mrc);
                h_hy_mrc(:, k_hy_mrc, 3) = j_h1_mrc;
            end
        end

        % choose the second largest
        for j_h2_mrc = 1 : L
            if j_h2_mrc == h_hy_mrc(:, k_hy_mrc, 3)
                continue
            end

            if j_h2_mrc == 1
                h_hy_mrc(:, k_hy_mrc, 2) = h_mrc(:, k_hy_mrc, j_h2_mrc);
                n_hy_mrc(:, k_hy_mrc, 2) = n_mrc(:, k_hy_mrc, j_h2_mrc);
                continue
            end

            if abs(h_mrc(:, k_hy_mrc, j_h1_mrc)) > abs(h_hy_mrc(:, k_hy_mrc, 2))
                h_hy_mrc(:, k_hy_mrc, 2) = h_mrc(:, k_hy_mrc, j_h1_mrc);
                n_hy_mrc(:, k_hy_mrc, 2) = n_mrc(:, k_hy_mrc, j_h1_mrc);
            end
        end
    end
    h_hy_mrc(:, :, 3) = [];

    y_hy_mrc = h_hy_mrc .* s + deviate_hy_mrc * n_hy_mrc;
    r_hy_mrc = conj(h_hy_mrc) .* y_hy_mrc;

    sigma = zeros(1, N);
    for j_hy_mrc = 1 : 2
        sigma = sigma + r_hy_mrc(:, :, j_hy_mrc);
    end
    xHat_hy_mrc = real(sigma) > 0;
    error_hy_mrc(i_hy_mrc) = size(find([x - xHat_hy_mrc]),2);
end
errorbitrate_hy_mrc = error_hy_mrc / N;

% MRC for comparison ------------------------------------------------------
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

% picture -----------------------------------------------------------------
figure(1);
semilogy(SNRdB, errorbitrate_hy_mrc, '-o');
    hold on
semilogy(SNRdB, errorbitrate_mrc, '-h');

axis([0 10 10^-6 10^0])
grid on

xlabel('Eb/No');
ylabel('BitError Probability');