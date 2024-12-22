clc;
close all;

% Parametri
N = 64; % Broj OFDM podnosioca
M = 16; % 16-QAM
numSymbols = 1000; % Broj OFDM simbola
EbNo = 0:2:20; % Eb/No vrijednosti u dB
symbolRate = 1e6; % Brzina simbola

fs = N * symbolRate; % Frekvencija uzorkovanja
fiberLength = 10e3; % Dužina vlakna (10 km)
beta2 = -2.17e-26; % Koeficijent hromatske disperzije (s^2/m)
berSimulated = zeros(size(EbNo));
berTheoretical = zeros(size(EbNo));
 
% QAM modulacija i demodulacija objekti
qamMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true);
qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
 
for i = 1:length(EbNo)
    % Generisanje slucajnih podataka
    data = randi([0 1], N*log2(M), numSymbols);
    
    % 16-QAM modulacija
    modData = qamMod(data(:));
    
    % IFFT - zbog OFDM
    ifftData = ifft(reshape(modData, N, numSymbols), N, 1);
    
    % Dodavanje CP (Cyclic Prefix)
    cpLength = 16;
    txData = [ifftData(end-cpLength+1:end,:); ifftData];
    
    % DCO
    dcBias = 0.6; % optimalno izabrano 
    txDataOptical = (txData) + dcBias; % Dodavanje DC biasa

    %Optički kanal
    txData_after_fiber = optical_channel(txDataOptical, fiberLength, beta2, fs);
    rxData = txData_after_fiber;
    
    % Dodavanje šuma
    snr = EbNo(i) + 10*log10(log2(M)) - 10*log10(N/(N+cpLength));
    rxData = awgn(rxData, snr, 'measured');
    
    % Uklanjanje CP
    rxData = rxData(cpLength+1:end, :); %ono sto dodamo na predaji moramo eliminsati na prijemu
    
    % FFT
    fftData = fft(rxData-dcBias, N, 1); %ono sto dodamo na predaji moramo eliminsati na prijemu
    
    % 16-QAM demodulacija
    demodData = qamDemod(fftData(:));
    
    % BER racunanje
    [numErrors, ber] = biterr(data(:), demodData);
    berSimulated(i) = ber;
    
    % Teorijski BER za 16-QAM u kanalu sa AWGN
    berTheoretical(i) = berawgn(EbNo(i), 'qam', M);
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
 
sgtitle('O-OFDM 16QAM - koristeci gotove funkcije za modulaciju');

% Konstelacijski dijagram
scatterplot(modData(:)); grid on;
title('\rmKonstelacijski dijagram - predaja');
scatterplot(fftData(:));
title('\rmKonstelacijski dijagram - prijem');
grid on;
 
% PSD prikaz
[psd,f] = periodogram(txDataOptical(:), hamming(length(txDataOptical(:))), 2^14, fs, 'centered');
figure;
subplot(2,1,1);
plot(f/10^6,10*log10(psd)); title('PSD - DCO OOFDM');
xlabel('Frequency (MHz)');
ylabel('Power/Frequency (dB/Hz)'); ylim([-150 -40]);
grid on;
hold on;
subplot(2,1,2);
[psd2,f] = periodogram(txData(:), hamming(length(txData(:))), 2^14, fs, 'centered');
plot(f/10^6,10*log10(psd2))
title('PSD - OFDM');
xlabel('Frequency (MHz)');
ylabel('Power/Frequency (dB/Hz)'); ylim([-150 -40]);
grid on;
sgtitle('SGS za OFDM 16QAM');