clc;
clear all;
close all;

format long;
M = 16;  % 16-PSK
k = log2(M);  % Broj bita po simbolu
EbNo = (5:15)';  % Raspon Eb/No
symbol_per_frame = 100;  % Broj simbola po okviru
data_points = 144;  % Broj podataka
block_size = 8;  % Veličina bloka
cp = ceil(0.1*block_size);  % Cyclic prefix
ifft_points = block_size;
fft_points = block_size;
n = 3e4;  % Broj bita
nsamp = 1;
N = 10^5;

% Generisanje binarnih podataka
data = randi([0 1], N, 1);  % Binarni podaci

% Grupisanje binarnih podataka u 4 bita po simbolu za 16-PSK
data_symbols = reshape(data, k, length(data)/k).';  % Grupisanje u matrice
decimal_data = bi2de(data_symbols, 'left-msb');  % Pretvaranje u dekadni oblik

% 16-PSK modulacija (ručno)
% Generisanje faza za 16-PSK
phase_angles = (2*pi*(0:M-1))/M;  % Faze: 0, 2*pi/16, ..., 2*pi*(M-1)/M
psk_modulisani = exp(1j * phase_angles(decimal_data + 1));  % Dodajemo 1 jer indeksi počinju od 0

% Konstelacijski dijagram
scatterplot(psk_modulisani);
title('Konstelacijski dijagram na predajnoj strani');

% IFFT
broj_kolona = (length(psk_modulisani)/block_size); % IFFT
data_matrix = reshape(psk_modulisani, block_size, broj_kolona);
cp_start = block_size - cp;
cp_end = block_size;
for i = 1:broj_kolona
    ifft_data_matrix(:, i) = ifft(data_matrix(:, i), ifft_points);
    for j = 1:cp
       actual_cp(j, i) = ifft_data_matrix(j + cp_start, i);
    end
    ifft_data(:, i) = vertcat(actual_cp(:, i), ifft_data_matrix(:, i));
end

[rows_ifft_data, cols_ifft_data] = size(ifft_data); % Konverzija u seriju
len_ofdm_data = rows_ifft_data * cols_ifft_data;
ofdm_signal = reshape(ifft_data, 1, len_ofdm_data);
figure(3)
plot(real(ofdm_signal)); xlabel('Time'); ylabel('Amplitude');
title('OFDM Signal'); grid on;

% clipping - PAPR redukcija
avg = 0;
clipped = ofdm_signal;
for i = 1:length(clipped)
    if clipped(i) > avg
        clipped(i) = clipped(i);
    elseif clipped(i) < -avg
        clipped(i) = 0;
    end
end

figure(4)
plot(real(clipped)); xlabel('Vrijeme'); ylabel('Amplituda');
title('Clipped signal'); grid on;
axis([0 180 -1.5 1.5]);

% Kanal
channel = randn(1, block_size) + sqrt(-1) * randn(1, block_size);

% Prijemnik Rx
after_channel = filter(channel, 1, ofdm_signal);
EbNo2 = 10;
snr = EbNo2 + 10 * log10(k) - 10 * log10(nsamp);
receivedSignal = awgn(psk_modulisani, snr, 'measured');
sPlotFig = scatterplot(receivedSignal, 1, 0, 'g.');
hold on
scatterplot(psk_modulisani, 1, 0, 'k*', sPlotFig)
awgn_noise = awgn(zeros(1, length(after_channel)), 0);
primljeni_signal = awgn_noise + after_channel;
primljeni_signal_matrix = reshape(primljeni_signal, rows_ifft_data, cols_ifft_data); % Konverzija u paralelu zbog FFT
primljeni_signal_matrix(1:cp, :) = []; % Uklanjanje CP
for i = 1:cols_ifft_data % FFT
    fft_data_matrix(:, i) = fft(primljeni_signal_matrix(:, i), fft_points);
end

primljeni_clipped = reshape(fft_data_matrix, 1, (block_size * broj_kolona)); % Konverzija u seriju

% 16-PSK demodulacija (ručno)
received_angles = angle(primljeni_clipped);  % Izračunavanje faze svakog prijemnog simbola
demodulated_data = mod(round(received_angles * M / (2*pi)), M);  % Pronađite indeks najbliže faze

% Prikazivanje demodulisanih podataka
figure(7)
stem(demodulated_data, 'rx');
grid on; xlabel('Podaci');
title('Primljena sekvenca');

% BER
dataOutMatrix = de2bi(demodulated_data, k);
dataOut = dataOutMatrix(:);
[numErrors, ber] = biterr(data, dataOut);
fprintf('\nBER = %5.2e, broj gresaka %d \n', ber, numErrors);

% Prijemnik clipped signala
after_channel = filter(channel, 1, clipped);
awgn_noise = awgn(zeros(1, length(after_channel)), 0);
primljeni_signal = awgn_noise + after_channel;
primljeni_signal_matrix = reshape(primljeni_signal, rows_ifft_data, cols_ifft_data); % Konverzija u paralelu zbog FFT
primljeni_signal_matrix(1:cp, :) = []; % Uklanjanje CP
for i = 1:cols_ifft_data % FFT
    fft_data_matrix(:, i) = fft(primljeni_signal_matrix(:, i), fft_points);
end
primljeni_clipped = reshape(fft_data_matrix, 1, (block_size * broj_kolona)); % Konverzija u seriju
psk_demodulated_data = pskdemod(primljeni_clipped, M);  % 16-PSK demodulacija

figure(8)
stem(psk_demodulated_data, 'rx');

% Teorijski 16-PSK BER
EbNodB2 = 0:2:16;
EbNo2 = 10.^(EbNodB2 / 10);
x = sqrt(2 * sin(pi / M) * EbNo2);
Pb = (1 / M) * erfc(x / sqrt(2));
figure(9)
semilogy(EbNodB2, Pb, 'bs-')
title('16-PSK BER');
xlabel('EbNo');
ylabel('Pb');
