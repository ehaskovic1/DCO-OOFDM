clear all
close all
clc
% % Parametri
% N = 64; % Broj OFDM podnosaèa
% M = 16; % 16-PSK
% numSymbols = 1000; % Broj OFDM simbola
% EbNo = 0:2:20; % Eb/No vrijednosti u dB
%  
% berSimulated = zeros(size(EbNo));
% berTheoretical = zeros(size(EbNo));
%  
% % PSK modulacija i demodulacija objekti
% %pskMod = comm.PSKModulator('ModulationOrder', M, 'BitInput', false);
% %pskDemod = comm.PSKDemodulator('ModulationOrder', M, 'BitOutput', false);
% rng(10);
% data = randi([0 M-1], N, numSymbols);
% % PSK modulacija
% modData = pskmod(data(:), M);
% for i = 1:length(EbNo)
%     ifftData = ifft(reshape(modData, N, numSymbols), N, 1);
%     
%     % Dodavanje CP (Cyclic Prefix)
%     cpLength = 16;
%     txData = [ifftData(end-cpLength+1:end,:); ifftData];
%     
%     % Dodavanje AWGN šuma
%     snr = EbNo(i) + 10*log10(log2(M)) - 10*log10(N/(N+cpLength));
%     rxData = awgn(txData, snr, 'measured');
%     
%     % Uklanjanje CP
%     rxData = rxData(cpLength+1:end, :);
%     
%     % FFT
%     fftData = fft(rxData, N, 1);
    
    % PSK demodulacija
 %   demodData = pskdemod(fftData(:), M);
    % Prevoðenje podataka u bitove
%     dataBits = de2bi(data(:), log2(M), 'left-msb');
%     demodDataBits = de2bi(demodData(:), log2(M), 'left-msb');
%     
%     % BER izraèunavanje
%     [numErrors, ber] = biterr(dataBits(:), demodDataBits(:));
%     berSimulated(i) = ber;
%     
%     % Teorijski BER za 16-PSK
%     berTheoretical(i) = berawgn(EbNo(i), 'psk', M, 'nondiff');
%end
%  
% % Grafièki prikazi
% figure;
%  
% % BER prikaz
% subplot(2,2,1);
% semilogy(EbNo, berSimulated, 'b-o');
% hold on;
% semilogy(EbNo, berTheoretical, 'r-*');
% title('Bit Error Rate (BER)');
% xlabel('Eb/No (dB)');
% ylabel('BER');
% legend('Simulacijski BER', 'Teorijski BER');
% grid on;
%  
% % SNR prikaz
% subplot(2,2,2);
% plot(EbNo, 10*log10(1./berSimulated), 'b-o');
% title('Signal-to-Noise Ratio (SNR)');
% xlabel('Eb/No (dB)');
% ylabel('SNR (dB)');
% grid on;
%  
% % Poslani signal
% subplot(2,2,3);
% plot(real(txData(:,1)));
% title('Poslani Signal');
% xlabel('Vrijeme');
% ylabel('Amplituda');
%  
% % Primljeni signal
% subplot(2,2,4);
% plot(real(rxData(:,1)));
% title('Primljeni Signal');
% xlabel('Vrijeme');
% ylabel('Amplituda');
%  
%  
% % Konstelacijski dijagram
% scatterplot(fftData(:));
% title('Konstelacijski Dijagram');
% grid on;
%  
% % PSD prikaz
% [psd, f] = periodogram(txData(:),[],[],1);
% figure;
% plot(f,10*log10(psd));
% title('Power Spectral Density (PSD)');
% xlabel('Frequency (Hz)');
% ylabel('Power/Frequency (dB/Hz)');
% grid on;

% %OFDM 16QAM
% 
% % Parametri
% N = 64; % Broj OFDM podnosaèa
% M = 16; % 16-QAM
% numSymbols = 1000; % Broj OFDM simbola
% EbNo = 0:2:20; % Eb/No vrijednosti u dB
%  
% berSimulated = zeros(size(EbNo));
% berTheoretical = zeros(size(EbNo));
%  
% % QAM modulacija i demodulacija objekti
% qamMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true);
% qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
%  
% for i = 1:length(EbNo)
%     % Generisanje sluèajnih podataka
%     data = randi([0 1], N*log2(M), numSymbols);
%     
%     % 16-QAM modulacija
%     modData = qamMod(data(:));
%     
%     % IFFT
%     ifftData = ifft(reshape(modData, N, numSymbols), N, 1);
%     
%     % Dodavanje CP (Cyclic Prefix)
%     cpLength = 16;
%     txData = [ifftData(end-cpLength+1:end,:); ifftData];
%     
%     % Dodavanje AWGN šuma
%     snr = EbNo(i) + 10*log10(log2(M)) - 10*log10(N/(N+cpLength));
%     rxData = awgn(txData, snr, 'measured');
%     
%     % Uklanjanje CP
%     rxData = rxData(cpLength+1:end, :);
%     
%     % FFT
%     fftData = fft(rxData, N, 1);
%     
%     % 16-QAM demodulacija
%     demodData = qamDemod(fftData(:));
%     
%     % BER izraèunavanje
%     [numErrors, ber] = biterr(data(:), demodData);
%     berSimulated(i) = ber;
%     
%     % Teorijski BER za 16-QAM
%     berTheoretical(i) = berawgn(EbNo(i), 'qam', M);
% end
%  
% % Grafièki prikazi
% figure;
%  
% % BER prikaz
% subplot(2,2,1);
% semilogy(EbNo, berSimulated, 'b-o');
% hold on;
% semilogy(EbNo, berTheoretical, 'r-*');
% title('Bit Error Rate (BER)');
% xlabel('Eb/No (dB)');
% ylabel('BER');
% legend('Simulacijski BER', 'Teorijski BER');
% grid on;
%  
% % SNR prikaz
% subplot(2,2,2);
% plot(EbNo, 10*log10(1./berSimulated), 'b-o');
% title('Signal-to-Noise Ratio (SNR)');
% xlabel('Eb/No (dB)');
% ylabel('SNR (dB)');
% grid on;
%  
% % Poslani signal
% subplot(2,2,3);
% plot(real(txData(:,1)));
% title('Poslani Signal');
% xlabel('Vrijeme');
% ylabel('Amplituda');
%  
% % Primljeni signal
% subplot(2,2,4);
% plot(real(rxData(:,1)));
% title('Primljeni Signal');
% xlabel('Vrijeme');
% ylabel('Amplituda');
%  
%  
%  
% % Konstelacijski dijagram
% scatterplot(fftData(:));
% title('Konstelacijski Dijagram');
% grid on;
%  
% % PSD prikaz
% [psd, f] = periodogram(txData(:),[],[],1);
% figure;
% plot(f,10*log10(psd));
% title('Power Spectral Density (PSD)');
% xlabel('Frequency (Hz)');
% ylabel('Power/Frequency (dB/Hz)');
% grid on;

% %OFDM 16QAM MANUELLNI
% 
% clc;
% clear all;
% close all;
%  
%  
%  
% % Parametri
% N = 64; % Broj OFDM podnosaèa
% M = 16; % 16-QAM
% numSymbols = 1000; % Broj OFDM simbola
% EbNo = 0:2:20; % Eb/No vrijednosti u dB
%  
% berSimulated = zeros(size(EbNo));
% berTheoretical = zeros(size(EbNo));
%  
% for i = 1:length(EbNo)
%     % Generisanje sluèajnih podataka
%     data = randi([0 M-1], N, numSymbols);
%     
%     % 16-QAM modulacija
%     I = 2*mod(data, sqrt(M)) - sqrt(M) + 1;
%     Q = 2*floor(data/sqrt(M)) - sqrt(M) + 1;
%     modData = I + 1j*Q;
%     
%     % IFFT
%     ifftData = ifft(reshape(modData, N, numSymbols), N, 1);
%     
%     % Dodavanje CP (Cyclic Prefix)
%     cpLength = 16;
%     txData = [ifftData(end-cpLength+1:end,:); ifftData];
%     
%     % Dodavanje AWGN šuma
%     snr = EbNo(i) + 10*log10(log2(M)) - 10*log10(N/(N+cpLength));
%     rxData = awgn(txData, snr, 'measured');
%     
%     % Uklanjanje CP
%     rxData = rxData(cpLength+1:end, :);
%     
%     % FFT
%     fftData = fft(rxData, N, 1);
%     
%     % 16-QAM demodulacija
%     I = real(fftData);
%     Q = imag(fftData);
%     demodData = floor((I + sqrt(M) - 1)/2) + sqrt(M) * floor((Q + sqrt(M) - 1)/2);
%     demodData = mod(demodData, M); % Map values to nearest symbol
%     
%     % Prevoðenje podataka u bitove
%     dataBits = de2bi(data(:), log2(M), 'left-msb');
%     demodDataBits = de2bi(demodData(:), log2(M), 'left-msb');
%     
%     % BER izraèunavanje
%     [numErrors, ber] = biterr(dataBits(:), demodDataBits(:));
%     berSimulated(i) = ber;
%     
%     % Teorijski BER za 16-QAM
%     berTheoretical(i) = berawgn(EbNo(i), 'qam', M);
% end
%  
% % Grafièki prikazi
% figure;
%  
% % BER prikaz
% subplot(2,2,1);
% semilogy(EbNo, berSimulated, 'b-o');
% hold on;
% semilogy(EbNo, berTheoretical, 'r-*');
% title('Bit Error Rate (BER)');
% xlabel('Eb/No (dB)');
% ylabel('BER');
% legend('Simulacijski BER', 'Teorijski BER');
% grid on;
%  
% % SNR prikaz
% subplot(2,2,2);
% plot(EbNo, 10*log10(1./berSimulated), 'b-o');
% title('Signal-to-Noise Ratio (SNR)');
% xlabel('Eb/No (dB)');
% ylabel('SNR (dB)');
% grid on;
%  
% % Poslani signal
% subplot(2,2,3);
% plot(real(txData(:,1)));
% title('Poslani Signal');
% xlabel('Vrijeme');
% ylabel('Amplituda');
%  
% % Primljeni signal
% subplot(2,2,4);
% plot(real(rxData(:,1)));
% title('Primljeni Signal');
% xlabel('Vrijeme');
% ylabel('Amplituda');
%  
%  
% % Konstelacijski dijagram
% scatterplot(fftData(:));
% title('Konstelacijski Dijagram');
% grid on;
%  
% % PSD prikaz
% [psd, f] = periodogram(txData(:),[],[],1);
% figure;
% plot(f,10*log10(psd));
% title('Power Spectral Density (PSD)');
% xlabel('Frequency (Hz)');
% ylabel('Power/Frequency (dB/Hz)');
% grid on;


% OFDM 16 PSK MANUELNI
% 
% Parametri
N = 64; % Broj OFDM podnosaèa
M = 16; % 16-PSK
numSymbols = 1000; % Broj OFDM simbola
EbNo = 0:2:20; % Eb/No vrijednosti u dB
 
berSimulated = zeros(size(EbNo));
berTheoretical = zeros(size(EbNo));
%Generisanje sluèajnih podataka
rng(10);
data = randi([0 M-1], N, numSymbols);
%16-PSK modulacija
theta = 2 * pi * data(:) / M;
modData = exp(1j * theta);
for i = 1:length(EbNo)
    %IFFT
    ifftData = ifft(reshape(modData, N, numSymbols), N, 1);
    
    %Dodavanje CP (Cyclic Prefix)
    cpLength = 16;
    txData = [ifftData(end-cpLength+1:end,:); ifftData];
    
    %Dodavanje AWGN šuma
    snr = EbNo(i) + 10*log10(log2(M)) - 10*log10(N/(N+cpLength));
    rxData = awgn(txData, snr, 'measured');
    
    %Uklanjanje CP
    rxData = rxData(cpLength+1:end, :);
    
    %FFT
    fftData = fft(rxData, N, 1);
    
    %16-PSK demodulacija
    rxPhase = angle(fftData);
    rxPhase(rxPhase<0)=rxPhase(rxPhase<0)+2*pi;
    demodData=round(M*rxPhase / (2*pi));
    demodData(demodData==M)=0;
    
    %Osiguranje da `demodData` bude unutar opsega [0, M-1]
    demodData(demodData < 0) = 0;
    demodData(demodData >= M) = M-1;
    demodData=demodData(:);
    
    %Prevoðenje podataka u bitove
    dataBits = de2bi(data(:), log2(M), 'left-msb');
    demodDataBits = de2bi(demodData(:), log2(M), 'left-msb');
    
    %BER izraèunavanje
    [numErrors, ber] = biterr(dataBits(:)', demodDataBits(:)');
    berSimulated(i) = ber;
    
    %Teorijski BER za 16-PSK
    berTheoretical(i) = berawgn(EbNo(i), 'psk', M, 'nondiff');
end
 
%Grafièki prikazi
figure;
 
%BER prikaz
subplot(2,2,1);
semilogy(EbNo, berSimulated, 'b-o');
hold on;
semilogy(EbNo, berTheoretical, 'r-*');
title('Bit Error Rate (BER)');
xlabel('Eb/No (dB)');
ylabel('BER');
legend('Simulacijski BER', 'Teorijski BER');
grid on;
 
%SNR prikaz
subplot(2,2,2);
plot(EbNo, 10*log10(1./berSimulated), 'b-o');
title('Signal-to-Noise Ratio (SNR)');
xlabel('Eb/No (dB)');
ylabel('SNR (dB)');
grid on;
 
%Poslani signal
subplot(2,2,3);
plot(real(txData(:,1)));
title('Poslani Signal');
xlabel('Vrijeme');
ylabel('Amplituda');
 
%Primljeni signal
subplot(2,2,4);
plot(real(rxData(:,1)));
title('Primljeni Signal');
xlabel('Vrijeme');
ylabel('Amplituda');
 
 
%Konstelacijski dijagram
scatterplot(fftData(:));
title('Konstelacijski Dijagram');
grid on;
 
%PSD prikaz
[psd, f] = periodogram(txData(:),[],[],1);
figure;
plot(f,10*log10(psd));
title('Power Spectral Density (PSD)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;
