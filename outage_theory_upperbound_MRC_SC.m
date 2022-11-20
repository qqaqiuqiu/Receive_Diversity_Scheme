clear all; clc;

SNRdB = 0:2:10;                     % SNR = [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30]
SNR = 10 .^ (0.1 .* SNRdB);         % SNR in linear scale
L = 4;

rho_th_dB = 6;
rho_th = 10 ^ (0.1 * rho_th_dB);

% MRC ---------------------------------------------------------------------

% Outage Probability
sigma = 0;
for k = 0 : L-1
    parta = 1 / factorial(k);
    partb = power((rho_th ./ SNR), k);
    sigma = sigma + parta * partb;
end
pout_mrc = 1 - exp(-rho_th ./ SNR) .* sigma;

% Upper bound
upper_mrc = power((rho_th ./ SNR), L) / factorial(L);

% SC ----------------------------------------------------------------------

% Outage Probability
pout_sc = power((1 - exp(-rho_th ./ SNR)), L);

% Upper bound
upper_sc = power((rho_th ./ SNR), L);

% picture------------------------------------------------------------------
figure(1);
semilogy(SNRdB, pout_mrc, '-o');
  hold on
semilogy(SNRdB, pout_sc, '-h');
semilogy(SNRdB, upper_mrc);
semilogy(SNRdB, upper_sc);

axis([0 10 10^-4 10^0])
grid on
L=legend('MRC theory', 'SC(SSC) theory', 'MRC upperbound', ...
    'SC(SSC) upperbound');
set(L,'Fontsize',12);

xlabel('Eb/No');
ylabel('Outage Probability');
title('Theory and Upper bound of Receive diversity Comparison');