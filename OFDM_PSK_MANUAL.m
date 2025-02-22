clc;
clear all;
close all;

% Parametri
N = 64; % Broj OFDM podnosioca
M = 16; % 16-PSK
numSymbols = 1000; % Broj OFDM simbola
EbNo = 0:2:20; % Eb/No vrijednosti u dB
symbolRate = 1e6; % Brzina simbola

fs = N * symbolRate; % Frekvencija uzorkovanja
fiberLength = 10e3; % Dužina vlakna (10 km)
beta2 = -2.17e-26; % Koeficijent hromatske disperzije (s^2/m)
berSimulated = zeros(size(EbNo));
berTheoretical = zeros(size(EbNo));
numZeroCarriers = 16; % Broj nultih subnosioca

%Generisanje sluèajnih podataka
rng(10);
data = randi([0 M-1], N, numSymbols);

%16-PSK modulacija - rucno implementirana
theta = 2 * pi * data / M;
modData = exp(1j * theta);

for i = 1:length(EbNo)
     % Dodavanje nultih subnosioca
    modDataWithZeros = [modData(1:N/2, :); zeros(numZeroCarriers, numSymbols); modData(N/2+1:end, :)];

    % IFFT - zbog OFDM
    ifftData = ifft(modDataWithZeros, N + numZeroCarriers, 1);
    
    %IFFT
    %ifftData = ifft(reshape(modData, N, numSymbols), N, 1);
    
    %Dodavanje CP (Cyclic Prefix)
    cpLength = 16;
    txData = [ifftData(end-cpLength+1:end,:); ifftData];
    
    % DCO
    dcBias = 0.2;
    txDataOptical = (txData) + dcBias; 

    % Opticki kanal
    txData_after_fiber = optical_channel(txDataOptical, fiberLength, beta2, fs);
    rxData = txData_after_fiber;
    
    % Dodavanje šuma
    snr = EbNo(i) + 10*log10(log2(M)) - 10*log10(N/(N+cpLength));
    rxData = awgn(rxData, snr, 'measured');
    
    %Uklanjanje CP
    rxData = rxData(cpLength+1:end, :);
   
    % FFT
    fftData = fft(rxData - dcBias, N + numZeroCarriers, 1);

    % Uklanjanje nultih subnosioca iz sredine spektra
    fftData = [fftData(1:N/2, :); fftData(end-N/2+1:end, :)];
    
    %16-PSK demodulacija
    rxPhase = angle(fftData);
    rxPhase(rxPhase<0)=rxPhase(rxPhase<0)+2*pi;
    demodData=round(M*rxPhase / (2*pi));
    demodData(demodData==M)=0;
    
    %Osiguranje da `demodData` bude unutar opsega [0, M-1]
    demodData(demodData < 0) = 0;
    demodData(demodData >= M) = M-1;
    demodData=demodData(:);
    
    %Prebacivanje u bite
    dataBits = de2bi(data(:), log2(M), 'left-msb');
    demodDataBits = de2bi(demodData(:), log2(M), 'left-msb');
    
    %BER racunanje
    [numErrors, ber] = biterr(dataBits(:)', demodDataBits(:)');
    berSimulated(i) = ber;
    
    %Teorijski BER za 16-PSK
    berTheoretical(i) = berawgn(EbNo(i), 'psk', M, 'nondiff');
end
 
% Grafièki prikazi
figure;
 
% BER prikaz
subplot(2,2,1);
semilogy(EbNo, berSimulated, 'b-o');
hold on;
semilogy(EbNo, berTheoretical, 'r-*');
title('Grafik vjerovatnoce greske');
xlabel('Eb/No (dB)');
ylabel('BER'); ylim ([10^-4, 1]);
legend('Simulacijski BER', 'Teorijski BER - AWGN kanal');
grid on;
 
% SNR prikaz
subplot(2,2,2);
plot(EbNo, 10*log10(1./berSimulated), 'b-o');
title('Signal-to-Noise Ratio (SNR)');
xlabel('Eb/No (dB)');
ylabel('SNR (dB)');
grid on;
 
% Poslani signal
subplot(2,2,3);
plot(real(txData(:,1)));
title('Poslani Signal');
xlabel('Vrijeme');
ylabel('Amplituda');
 
% Primljeni signal
subplot(2,2,4);
plot(real(rxData(:,1)));
title('Primljeni Signal');
xlabel('Vrijeme');
ylabel('Amplituda');
 
sgtitle('O-OFDM 16PSK - rucno implementirane fje');

% Konstelacijski dijagram
scatterplot(modData(:), [],[], 'k.'); grid on;
set(gcf, 'Color', 'w'); % Postavlja bijelu pozadinu za cijelu figuru
set(gca, 'Color', 'w'); % Postavlja bijelu pozadinu za koordinatni sistem
title('\rmKonstelacijski dijagram - predaja', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

scatterplot(fftData(:),  [],[], 'b.');
title('\rmKonstelacijski dijagram - prijem', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
grid on; 
set(gcf, 'Color', 'w'); % Postavlja bijelu pozadinu za cijelu figuru
set(gca, 'Color', 'w'); % Postavlja bijelu pozadinu za koordinatni sistem
 
% PSD prikaz
[psd,f] = periodogram(txDataOptical(:), hamming(length(txDataOptical(:))), 2^14, fs, 'centered');
figure;
subplot(2,1,1);
plot(f/10^6,10*log10(psd)); title('PSD - DCO OOFDM');
xlabel('Frequency (MHz)');
ylabel('Power/Frequency (dB/Hz)'); ylim([-150 -40]);
grid on;
% hold on;
subplot(2,1,2);
[psd2,f] = periodogram(txData(:), hamming(length(txData(:))), 2^14, fs, 'centered');
plot(f/10^6,10*log10(psd2))
title('PSD - OFDM');
xlabel('Frequency (MHz)');
ylabel('Power/Frequency (dB/Hz)'); ylim([-150 -40]);
grid on;
sgtitle('SGS za rucno implementirani OFDM 16PSK');
