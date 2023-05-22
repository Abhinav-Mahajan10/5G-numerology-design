% Configure carrier
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 30;
carrier.CyclicPrefix = 'normal';
carrier.NSizeGrid = 273;

% Set the operating frequency and choose the phase noise model
simParameters = [];
% simParameters.Fc = 6e9; % Frequency in Hz
simParameters.PNModel = 'B'; % 'A' (TDoc R1-163984 Set A), 'B' (TDoc R1-163984 Set B), 'C' (TR 38.803)

% Get the sample rate
ofdmInfo = nrOFDMInfo(carrier);
sr = ofdmInfo.SampleRate;

% Phase noise level
foffsetLog = (2:0.2:log10(sr/2)); % Model offset from 1e4 Hz to sr/2 Hz
foffset = 10.^foffsetLog;         % Linear frequency offset
pn_PSD1 = hPhaseNoisePoleZeroModel(foffset,6e9,simParameters.PNModel); % dBc/Hz
pn_PSD2 = hPhaseNoisePoleZeroModel(foffset,28e9,simParameters.PNModel); % dBc/Hz
pn_PSD3 = hPhaseNoisePoleZeroModel(foffset,82e9,simParameters.PNModel); % dBc/Hz

% Set phase noise level
pnoise = comm.PhaseNoise('FrequencyOffset',foffset,'Level',pn_PSD1,'SampleRate',sr);
pnoise.RandomStream = "mt19937ar with seed";

% Visualize spectrum mask of phase noise
figure 
semilogx(foffset, pn_PSD1, 'LineWidth', 2) 
hold on;
semilogx(foffset, pn_PSD2, 'LineWidth', 2) 
hold on;
semilogx(foffset, pn_PSD3, 'LineWidth', 2)
hold off;

legend('F_{LO} = 6 GHz', 'F_{LO} = 28 GHz', 'F_{LO} = 82 GHz')
xlabel('Frequency offset (Hz)')
ylabel('PSD (dBc/Hz)')
title('Phase noise magnitude response')
grid on
