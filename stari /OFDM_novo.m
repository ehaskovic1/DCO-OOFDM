clc
clear all
close all

format long;
M = 16;
k = log2(M);
EbNo = (5:15)';
symbol_per_frame = 100;
data_points = 144;
block_size = 8;
cp = ceil(0.1*block_size);
ifft_points = block_size;
fft_points = block_size;
n = 3e4; % broj bita
nsamp = 1;
N = 10^5;

% Predajnik Tx
data = randi([0 1],N,1); % unipolarni
figure(1)
stem(data); grid on; xlabel('Ulazna sekvenca');
title('Poslani podaci')

qam_modulisani = qammod(data,M,'InputType','bit','UnitAveragePower',true); %QAM modulacija
scatterplot(qam_modulisani);
title('Konstelacijski dijagram na predajnoj strani');

broj_kolona = (length(qam_modulisani)/block_size); % IFFT
data_matrix = reshape(qam_modulisani, block_size, broj_kolona);
cp_start = block_size-cp;
cp_end = block_size;
for i = 1:broj_kolona
    ifft_data_matrix(:,i) = ifft((data_matrix(:,i)),ifft_points);
    for j=1:cp
       actual_cp(j,i) = ifft_data_matrix(j+cp_start,i);
    end
    ifft_data(:,i) = vertcat(actual_cp(:,i),ifft_data_matrix(:,i));
end

[rows_ifft_data, cols_ifft_data] = size(ifft_data); % Konverzija u seriju
len_ofdm_data = rows_ifft_data*cols_ifft_data;
ofdm_signal = reshape(ifft_data, 1, len_ofdm_data);
figure(3)
plot(real(ofdm_signal)); xlabel('Time'); ylabel('Amplitude');
title('OFDM Signal');grid on;

% clipping - PAPR redukcija
 avg=0;
 clipped=ofdm_signal;
 for i=1:length(clipped)
 	if clipped(i) > avg
 		clipped(i) = clipped(i);
    elseif clipped(i) < -avg
 		clipped(i) = 0;
    end
 end


 figure(4)
plot(real(clipped)); xlabel('Vrijeme'); ylabel('Amplituda');
 title('Clipped signal');grid on;
axis([0 180 -1.5 1.5]);

% Kanal
channel = randn(1,block_size) + sqrt(-1)*randn(1,block_size);

% Prijemnik Rx
after_channel = filter(channel, 1, ofdm_signal);
EbNo2 = 10;
snr = EbNo2 + 10*log10(k) - 10*log10(nsamp);
receivedSignal = awgn(qam_modulisani,snr,'measured');
sPlotFig = scatterplot(receivedSignal,1,0,'g.');
hold on
scatterplot(qam_modulisani,1,0,'k*',sPlotFig)
awgn_noise = awgn(zeros(1,length(after_channel)),0);
primljeni_signal = awgn_noise+after_channel;
primljeni_signal_matrix = reshape(primljeni_signal,rows_ifft_data, cols_ifft_data); % Konverzija u paralelu zbog FFT
primljeni_signal_matrix(1:cp,:)=[]; % Uklanjanje CP
for i=1:cols_ifft_data % FFT
    fft_data_matrix(:,i) = fft(primljeni_signal_matrix(:,i),fft_points);
end

primljeni_clipped = reshape(fft_data_matrix,1,(block_size*broj_kolona)); % Konverzija u seriju
qam_demodulated_data = qamdemod(primljeni_clipped,M);
figure(7)
stem(qam_demodulated_data,'rx');
grid on;xlabel('Podaci');
title('Primljena sekvenca')   

dataOutMatrix = de2bi(qam_demodulated_data,k);
dataOut = dataOutMatrix(:);

[numErrors,ber] = biterr(data,dataOut);
fprintf('\nBER = %5.2e, broj gresaka %d \n',ber,numErrors)


% Prijemnik clipped signala
after_channel = filter(channel, 1, clipped);
awgn_noise = awgn(zeros(1,length(after_channel)),0);
primljeni_signal = awgn_noise+after_channel;
primljeni_signal_matrix = reshape(primljeni_signal,rows_ifft_data, cols_ifft_data); % Konverzija u paralelu zbog FFT
primljeni_signal_matrix(1:cp,:)=[]; % Uklanjanje CP
for i=1:cols_ifft_data %   FFT
    fft_data_matrix(:,i) = fft(primljeni_signal_matrix(:,i),fft_points);
end
primljeni_clipped = reshape(fft_data_matrix, 1,(block_size*broj_kolona)); % Konverzija u seriju
qam_demodulated_data = qamdemod(primljeni_clipped,M);

figure(8)
stem(qam_demodulated_data,'rx');

%BER
 noe=0; 
 for ii= 1:1:length(primljeni_clipped)
   if data(ii)~=primljeni_clipped(ii)
    noe=noe+1;
   else
       noe=noe;
   end 
 end 
ber=( noe/length(data)); 

% Teorijski QAM BER
EbNodB2=0:2:16;
EbNo2=10.^(EbNodB2/10);
x=sqrt(3*k*EbNo2/(M-1));
Pb=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));
figure(9) 
semilogy(EbNodB2,Pb,'bs-')
title('QAM BER');
xlabel('EbNo');
ylabel('Pb');
