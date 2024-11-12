clc
clear all
close all

format long;
M = 16;                          %   Qam signal constellation
k = log2(M);     % Bits per symbol
EbNo = (5:15)';      % Eb/No values (dB)
numSymPerFrame = 100;   % Number of QAM symbols per frame
no_of_data_points = 144;        %   have 128 data points
block_size = 8;                 %   size of each ofdm block
cp = ceil(0.1*block_size);  %   length of cyclic prefix
no_of_ifft_points = block_size;           %   128 points for the FFT/IFFT
no_of_fft_points = block_size;
n = 3e4; % Number of bits to process
nsamp = 1; % Oversampling rate
berEst = zeros(size(EbNo));

% ***** Simulation M-QAM transmission over noise ***
% with using Monte Carlo simulation
q=4;
m=2^q; % M-QAM level power of 2 
loop=20; % Monte Carlo 
N=100000;  % Frame length (x_1 x_2 ... x_N)
SNRdB=0:15; % SNR in dB
SNR=10.^(SNRdB/10);
Rate= zeros(1, length(SNRdB)); %

%   +++++   TRANSMITTER    +++++
%   1.  Generate 1 x 128 vector of random data points
data_source=randi([0 1],N,1); % 1 or 0 
figure(1)
stem(data_source); grid on; xlabel('Data Points'); ylabel('transmitted data phase representation')
title('Transmitted Data "O"')

%   2.  Perform QAM modulation
qam_modulated_data = qammod(data_source,M,'InputType','bit','UnitAveragePower',true);
ytx =qam_modulated_data;
scatterplot(qam_modulated_data);title('MODULATED TRANSMITTED DATA');

%   3.  Do IFFT on each block
num_cols=(length(qam_modulated_data)/block_size);
data_matrix = reshape(qam_modulated_data, block_size, num_cols);
cp_start = block_size-cp;
cp_end = block_size;
for i=1:num_cols
    ifft_data_matrix(:,i) = ifft((data_matrix(:,i)),no_of_ifft_points);
    for j=1:cp
       actual_cp(j,i) = ifft_data_matrix(j+cp_start,i);
    end
    ifft_data(:,i) = vertcat(actual_cp(:,i),ifft_data_matrix(:,i));
end

%   4.  Convert to serial stream for transmission

[rows_ifft_data, cols_ifft_data]=size(ifft_data);
len_ofdm_data = rows_ifft_data*cols_ifft_data;
ofdm_signal = reshape(ifft_data, 1, len_ofdm_data);
figure(3)
plot(real(ofdm_signal)); xlabel('Time'); ylabel('Amplitude');
title('OFDM Signal');grid on;

%   ---------------------------------------------------------------
%     +++++   clipping as a PAPR reduction method    +++++
%   ---------------------------------------------------------------

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
plot(real(clipped)); xlabel('Time'); ylabel('Amplitude');
 title('clipped Signal');grid on;
axis([0 180 -1.5 1.5]);

%   --------------------------------
%     +++++   CHANNEL    +++++
%   --------------------------------

channel = randn(1,block_size) + sqrt(-1)*randn(1,block_size);


%   ------------------------------------------
%     +++++   RECEIVER    +++++
%   ------------------------------------------
%   1.  Pass the ofdm signal through the channel
after_channel = filter(channel, 1, ofdm_signal);
%SNR with respect to ber
EbNo2 = 10;
snr = EbNo2 + 10*log10(k) - 10*log10(nsamp);
receivedSignal = awgn(qam_modulated_data,snr,'measured');
sPlotFig = scatterplot(receivedSignal,1,0,'g.');
hold on
scatterplot(qam_modulated_data,1,0,'k*',sPlotFig)
%   2.   Add Noise

awgn_noise = awgn(zeros(1,length(after_channel)),0);

%   3.  Add noise to signal...

recvd_signal = awgn_noise+after_channel;

%   4.  Convert Data back to "parallel" form to perform FFT

recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);

%   5.  Remove CP

recvd_signal_matrix(1:cp,:)=[];

%   6.  Perform FFT

for i=1:cols_ifft_data
    
    %   FFT
    
    fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),no_of_fft_points);

end



%   7.  Convert to serial stream

recvd_serial_data = reshape(fft_data_matrix,1,(block_size*num_cols));

%   8.  Demodulate the data

qam_demodulated_data = qamdemod(recvd_serial_data,M);

figure(7)

stem(qam_demodulated_data,'rx');

grid on;xlabel('Data Points');ylabel('received data phase representation');title('Received Data "X"')   

dataOutMatrix = de2bi(qam_demodulated_data,k);
dataOut = dataOutMatrix(:);

[numErrors,ber] = biterr(data_source,dataOut);
fprintf('\nThe binary coding bit error rate = %5.2e, based on %d errors\n',ber,numErrors)


%   ----------------------------------------------------
%   F:  %   +++++   RECEIVER of clipped signal    +++++
%   ----------------------------------------------------

%   1.  Pass the ofdm signal through the channel

after_channel = filter(channel, 1, clipped);

%   2.   Add Noise

awgn_noise = awgn(zeros(1,length(after_channel)),0);

%   3.  Add noise to signal...

recvd_signal = awgn_noise+after_channel;

%   4.  Convert Data back to "parallel" form to perform FFT

recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);

%   5.  Remove CP

recvd_signal_matrix(1:cp,:)=[];

%   6.  Perform FFT

for i=1:cols_ifft_data

    %   FFT
    
    fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),no_of_fft_points);

end

%   7.  Convert to serial stream

recvd_serial_data = reshape(fft_data_matrix, 1,(block_size*num_cols));

%   8.  Demodulate the data

qam_demodulated_data = qamdemod(recvd_serial_data,M);

figure(8)

stem(qam_demodulated_data,'rx');

%BER
 noe=0; 
 for ii= 1:1:length(recvd_serial_data)
   if data_source(ii)~=recvd_serial_data(ii)
    noe=noe+1;
   else
       noe=noe;
   end 
 end 
% ber= noe/length(ofdm_signal); 
ber=( noe/length(data_source)); 

% CALCULATE QAM BER USING FORMULA
EbNodB2=0:2:16;
EbNo2=10.^(EbNodB2/10);
x=sqrt(3*k*EbNo2/(M-1));
Pb=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));
figure(9) 
semilogy(EbNodB2,Pb,'bs-')
title('QAM bit error rate');
xlabel('EbNo');
ylabel('Pb');
