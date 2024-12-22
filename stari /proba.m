clear all;
close all;

format long;

M = 16;                 %  QAM
k = log2(M);            % broj bita po simbolu
EbNoVec = (5:15)';      % Eb/No (dB)
numSymPerFrame = 100;   % Number of QAM symbols per frame

no_of_data_points = 144;        %   have 128 data points
block_size = 8;                 %   size of each ofdm block
cp_len = ceil(0.1*block_size);  %   length of cyclic prefix
no_of_ifft_points = block_size;           %   128 points for the FFT/IFFT
no_of_fft_points = block_size;


n = 3e4; % Number of bits to process
nsamp = 1; % Oversampling rate


berEst = zeros(size(EbNoVec));

% ***** Simulation M-QAM transmission over noise ***%
% with using Monte Carlo simulation


q=4;
m=2^q; % M-QAM level power of 2 
loop=20; % Monte Carlo 
N=100000;  % Frame length (x_1 x_2 ... x_N)
SNRdB=0:15; % SNR in dB
SNR=10.^(SNRdB/10);
Rate= zeros(1, length(SNRdB)); %





%   ---------------------------------------------
%               PREDAJNIK - TX
%   ---------------------------------------------

data_source=randi([0 1],N,1); % unipolarni signal - optika
figure(1)
stem(data_source); grid on; 
xlabel('Data Points'); ylabel('transmitted data phase representation')
title('Poslana sekvenca')

qam_modulated_data = qammod(data_source,M,'InputType','bit','UnitAveragePower',true);
ytx =qam_modulated_data;
scatterplot(qam_modulated_data);
title('Modulisani poslani signal');

num_cols=(length(qam_modulated_data)/block_size);
data_matrix = reshape(qam_modulated_data, block_size, num_cols);

cp_start = block_size-cp_len;
cp_end = block_size;

for i=1:num_cols
    ifft_data_matrix(:,i) = ifft((data_matrix(:,i)),no_of_ifft_points);
    %   Compute and append Cyclic Prefix
    for j=1:cp_len
       actual_cp(j,i) = ifft_data_matrix(j+cp_start,i);
    end
    %   Append the CP to the existing block to create the actual OFDM block
    ifft_data(:,i) = vertcat(actual_cp(:,i),ifft_data_matrix(:,i));
end

[rows_ifft_data, cols_ifft_data]=size(ifft_data);
len_ofdm_data = rows_ifft_data*cols_ifft_data;

ofdm_signal = reshape(ifft_data, 1, len_ofdm_data);
figure(3)
plot(real(ofdm_signal)); xlabel('Time'); ylabel('Amplitude');
title('OFDM Signal'); grid on;


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
title('clipped Signal'); grid on;
axis([0 180 -1.5 1.5]);


%   --------------------------------
%               KANAL
%   --------------------------------


channel = randn(1,block_size) + sqrt(-1)*randn(1,block_size);


%   ------------------------------------------
%               PRIJEMNIK - RX
%   ------------------------------------------

after_channel = filter(channel, 1, ofdm_signal);
EbNo2 = 10;
snr = EbNo2 + 10*log10(k) - 10*log10(nsamp);
receivedSignal = awgn(qam_modulated_data,snr,'measured');
sPlotFig = scatterplot(receivedSignal,1,0,'g.');
hold on
scatterplot(qam_modulated_data,1,0,'k*',sPlotFig)

awgn_noise = awgn(zeros(1,length(after_channel)),0);

recvd_signal = awgn_noise+after_channel;

recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);

recvd_signal_matrix(1:cp_len,:)=[]; % uklanjanje ciklicnog prefiksa


for i=1:cols_ifft_data    
    fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),no_of_fft_points);
end


recvd_serial_data = reshape(fft_data_matrix,1,(block_size*num_cols));


qam_demodulated_data = qamdemod(recvd_serial_data,M); % demodulacija

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


after_channel = filter(channel, 1, clipped);

awgn_noise = awgn(zeros(1,length(after_channel)),0);

recvd_signal = awgn_noise+after_channel;

recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);

recvd_signal_matrix(1:cp_len,:)=[];

for i=1:cols_ifft_data
    
    fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),no_of_fft_points);

end

recvd_serial_data = reshape(fft_data_matrix, 1,(block_size*num_cols));

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
 
ber=( noe/length(data_source)); 
 
EbNodB2=0:2:16;
EbNo2=10.^(EbNodB2/10);
x=sqrt(3*k*EbNo2/(M-1));
Pb=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));
figure(9) 
semilogy(EbNodB2,Pb,'bs-')
title('QAM bit error rate');
xlabel('EbNo');
ylabel('Pb');


% ***** Simulation M-QAM transmission over noise ***%
% with using Monte Carlo simulation


% ******* Transmitter **************%
for dB= 1: length(SNRdB) % start looping by SNRdB
    for lp= 1: loop % start looping of frame data 
	
% ******* q-QAM signal generation **********%    
	x_inp=randi([0 1],N,1); % 1 or 0
    x_inp_mod=qammod(x_inp,q);

% ******* Channel **************%
	  
    y_channel=awgn(x_inp_mod,SNRdB(dB)); %  AWGN 
    % ******* Receiver ***************% 
    y=y_channel; 
    x_inp_dem=qamdemod(y,q);
    x_out=round(x_inp_dem);
    
    % ******* Bit Error Rate (BER) calulation ******%    
    
    [err, rate]= biterr(x_inp, x_out);
    Rate(dB)= Rate(dB) + rate;
    
    end % end for loop
   
    Rate(dB)= Rate(dB)/loop; % Average value over Monte Carlo simulation 
                              % loop
   
end % end Monte Carlo


% ******* Plot the simulation result *********%
    f1 = figure(12);
    set(f1,'color',[1 1 1]);
    semilogy(SNRdB,Rate, 'b-*')
    hold on;
    BER_th= (2*(sqrt(m)-1)/sqrt(m))*qfunc(sqrt((6*q/(m-1)))*sqrt(SNR)); % theoritical calculation for BER
    semilogy(SNRdB,BER_th,'r-o');
    hold on;
    axis([0 12 0.00000001  1.2]);  
    xlabel( 'Signal-to-Noise Ratio (SNR)')
    ylabel( 'Bit Error Rate (BER)')
    title('Simulation QAM transmission over noise');
    legend('BER simulation','BER calculation')
    grid on;
