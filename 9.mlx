% Parameters
nFFT = 64;          % FFT size
nDSC = 52;          % Number of data subcarriers
nSym = 1e4;         % Number of OFDM symbols
EbN0dB = 0:10;      % Eb/N0 range in dB
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % Es/N0

% Preallocate
nErr = zeros(1, length(EbN0dB));

for ii = 1:length(EbN0dB)
    % Transmitter
    ipBit = rand(1, nDSC * nSym) > 0.5; % Random bits
    ipMod = 2 * ipBit - 1;              % BPSK modulation
    ipMod = reshape(ipMod, nDSC, nSym).'; % Reshape to symbols
    
    % Map to subcarriers
    xF = [zeros(nSym, 6), ipMod(:, 1:nDSC/2), zeros(nSym, 1), ipMod(:, nDSC/2+1:nDSC), zeros(nSym, 5)];
    
    % IFFT and normalize power
    xt = (nFFT/sqrt(nDSC)) * ifft(fftshift(xF.')).';
    
    % Add cyclic prefix
    xt = [xt(:, 49:64), xt];
    xt = reshape(xt.', 1, nSym * 80);
    
    % Noise
    nt = (1/sqrt(2)) * (randn(1, nSym * 80) + 1j * randn(1, nSym * 80));
    yt = sqrt(80/64) * xt + 10^(-EsN0dB(ii)/20) * nt;
    
    % Receiver
    yt = reshape(yt.', 80, nSym).';
    yt = yt(:, 17:80); % Remove cyclic prefix
    
    % FFT
    yF = (sqrt(nDSC)/nFFT) * fftshift(fft(yt.')).';
    yMod = yF(:, [6+(1:nDSC/2), 7+(nDSC/2+1:nDSC)]);
    
    % Demodulate
    ipModHat = sign(real(yMod)); % BPSK demodulation
    ipBitHat = (ipModHat + 1) / 2;
    ipBitHat = reshape(ipBitHat.', nDSC * nSym, 1).';
    
    % Count errors
    nErr(ii) = sum(ipBitHat ~= ipBit);
end

% BER calculation
simBer = nErr / (nSym * nDSC);
theoryBer = 0.5 * erfc(sqrt(10.^(EbN0dB/10)));

% Plot
figure;
semilogy(EbN0dB, theoryBer, 'bs-', 'LineWidth', 2);
hold on;
semilogy(EbN0dB, simBer, 'mx-', 'LineWidth', 2);
axis([0 10 1e-5 1]);
grid on;
legend('Theory', 'Simulation');
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate');
title('BER for BPSK-OFDM');
