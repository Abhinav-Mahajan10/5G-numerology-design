close all;

numGuards = 212;

% Plot power spectral density (PSD) over all subcarriers
[sumOFDMSpec1,fOFDM1] = ofdm(1);
sumOFDMSpec1 = sumOFDMSpec1/mean(sumOFDMSpec1(1+2*numGuards:end-2*numGuards));

figure; 
plot(fOFDM1,10*log10(sumOFDMSpec1)); 
hold on
[sumOFDMSpec2,fOFDM2] = ofdm(16);
sumOFDMSpec2 = sumOFDMSpec2/mean(sumOFDMSpec2(1+2*numGuards*16:end-2*numGuards*16));
plot(fOFDM2,10*log10(sumOFDMSpec2)); 
grid on
axis([-0.5 0.5 -60 10]);
legend('60KHz','15KHz');
xlabel('Normalized frequency'); 
ylabel('PSD (dBW/Hz)')
title(['Power Spectral Density of different numerologies - OFDM']);

% Compute peak-to-average-power ratio (PAPR)
% PAPR2 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
% [~,~,paprOFDM] = step(PAPR2,ifftOut);
% disp(['Peak-to-Average-Power-Ratio (PAPR) for OFDM = ' num2str(paprOFDM) ' dB']);

function [sumOFDMSpec,fOFDM] = ofdm(N)

s = rng(211);            % Set RNG state for repeatability
% N=32;
numFFT = 512*N;           % Number of FFT points
numGuards = 212*N;         % Guard bands on both sides
K = 4;                   % Overlapping symbols, one of 2, 3, or 4
numSymbols = 300;        % Simulation length in symbols

% Prototype filter
switch K
    case 2
        HkOneSided = sqrt(2)/2;
    case 3
        HkOneSided = [0.911438 0.411438];
    case 4
        HkOneSided = [0.971960 sqrt(2)/2 0.235147];
    otherwise
        return
end
% Build symmetric filter
Hk = [fliplr(HkOneSided) 1 HkOneSided];

% QAM symbol mapper
bitsPerSubCarrier = 2;
qamMapper = comm.RectangularQAMModulator(...
    'ModulationOrder', 2^bitsPerSubCarrier, ...
    'BitInput', true, ...
    'NormalizationMethod', 'Average power');

% Transmit-end processing
%   Initialize arrays
L = numFFT-2*numGuards;  % Number of complex symbols per OFDM symbol
KF = K*numFFT;
KL = K*L;
dataSubCar = zeros(L, 1);
dataSubCarUp = zeros(KL, 1);

sumFBMCSpec = zeros(KF*2, 1);
sumOFDMSpec = zeros(numFFT*2, 1);

numBits = bitsPerSubCarrier*L/2;    % account for oversampling by 2
inpData = zeros(numBits, numSymbols);
rxBits = zeros(numBits, numSymbols);
txSigAll = complex(zeros(KF, numSymbols));
symBuf = complex(zeros(2*KF, 1));

% OFDM Modulation with Corresponding Parameters
%
% For comparison, we review the existing OFDM modulation technique, using
% the full occupied band, however, without a cyclic prefix.

    for symIdx = 1:numSymbols
    
        inpData2 = randi([0 1], bitsPerSubCarrier*L, 1);
        modData = step (qamMapper, inpData2);
        
        symOFDM = [zeros(numGuards,1); modData; zeros(numGuards,1)];
        ifftOut = sqrt(numFFT).*ifft(ifftshift(symOFDM));

        [specOFDM,fOFDM] = periodogram(ifftOut, rectwin(length(ifftOut)), ...
            numFFT*2, 1, 'centered'); 
        sumOFDMSpec = sumOFDMSpec + specOFDM;
    
    end
end