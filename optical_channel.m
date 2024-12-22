function rxData = optical_channel(txData, fiberLength, beta2, fs)
    % Parametri optičkog kanala
    c = 3e8; % Brzina svjetlosti (m/s)
    lambda = 1550e-9; % Talasna dužina (m)
    f0 = c / lambda; % Centralna frekvencija (Hz)

    % Vrijeme simulacije i frekvencijski opseg
    [numSymbols, N] = size(txData);
    dt = 1 / fs; % Vremenski korak
    t = (0:N-1) * dt; % Vremenska osa
    f = (-N/2:N/2-1) * (fs / N); % Frekvencijska osa

    % Hromatska disperzija
    H_cd = exp(-1j * (2 * pi * f).^2 * beta2 * fiberLength / 2); % Disperzijski filter

    % Propagacija kroz kanal
    rxData = zeros(numSymbols, N);
    for i = 1:numSymbols
        freqData = fftshift(fft(txData(i, :))); % FFT + pomak frekvencijskog spektra
        freqData = freqData .* H_cd; % Primjena disperzije
        rxData(i, :) = ifft(ifftshift(freqData)); % Povratak u vremenski domen
    end
end
