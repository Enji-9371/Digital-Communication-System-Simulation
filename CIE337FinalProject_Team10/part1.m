% Signal generation 
duration = 1; % 1 second
messageFreq = input('Enter message frequency (Hz): '); % message frequency in Hz
numPoints = 10000; % number of time samples for this duration

% Generate time vector
timeVec = 0:duration/numPoints:duration - duration/numPoints;

% Generate the original signal
originalSignal = input('Enter the message signal: ');

% Get Sampling Frequency from User
samplingFreq = input('Enter Sampling Frequency (Hz): '); % sampling frequency

% Apply Functions
% Apply Sampler Function 
[sampledSignal, sampledZeros, sampledTimeVec] = Sampler(originalSignal, timeVec, samplingFreq);

% Apply Quantizer Function 
numLevels = input('Enter Number of Quantization Levels: ');
peakLevel = input('Enter Peak Quantization Level: ');
mu = input('Enter μ Value: ');
[quantizedSignal, quantizationLevels] = Quantizer(sampledSignal, numLevels, peakLevel);

% Apply Encoder Function
encodingType = input('Enter Type of Encoding:\n 1-Unipolar NRZ\n 2-Polar NRZ\n 3-Manchester: ');
encodedSignal = Encoder(quantizationLevels, encodingType);

% Apply Decoder Function 
decodedSignal = Decoder(encodedSignal, encodingType, numLevels, peakLevel);

% Apply Filter Function 
reconstructedSignal = Filter(decodedSignal, numPoints, samplingFreq, messageFreq);

% Plots
close all;

% 1 - Original Signal and Sampled Signal
figure;

% Subplot for Sampled Signal
subplot(2,1,1);
stem(sampledTimeVec, sampledSignal, 'c'); % Sampled signal in cyan
title("Sampled Signal (Fs=" + samplingFreq + " Hz)");
xlabel('t (s)');
ylabel('Amplitude');
grid on;

% Subplot for Original Signal
subplot(2,1,2);
plot(timeVec, originalSignal, 'b'); % Original signal in blue
title('Original Signal');
xlabel('t (s)');
ylabel('Amplitude');
grid on;



% 2 - Sampled Signal and Quantized Signal
figure;

% Subplot for Sampled Signal
subplot(2,1,1);
stem(sampledTimeVec, sampledSignal, 'b'); % Sampled signal in blue
title('Sampled Signal');
xlabel('t (s)');
ylabel('Amplitude');
grid on;

% Subplot for Quantized Signal
subplot(2,1,2);
stairs(sampledTimeVec, quantizedSignal, 'g'); % Quantized signal in green
title('Quantized Signal');
xlabel('t (s)');
ylabel('Amplitude');
grid on;


% 3 - Encoded Signal
figure;
stairs(encodedSignal, 'magenta');
title('Encoded Signal');
xlabel('Sample');
ylabel('Amplitude');
grid on;

% 4 - Original Signal vs Reconstructed Signal
figure;

% Subplot for Original Signal
subplot(2,1,1);
plot(timeVec, originalSignal, 'b'); % Original signal in blue
title('Original Signal');
xlabel('t(s)');
ylabel('Amplitude');
grid on;

% Subplot for Reconstructed Signal
subplot(2,1,2);
plot(timeVec, reconstructedSignal, 'm'); % Reconstructed signal in magenta
title('Reconstructed Signal');
xlabel('t(s)');
ylabel('Amplitude');
grid on;



% 5 - Frequency Domain Representation
freqAxis = -numPoints/2:numPoints/2-1/numPoints;

figure;
subplot(3,1,1);
plot(freqAxis, abs(fftshift(fft(originalSignal)))/numPoints, 'r'); % FFT of Original Signal in red
title('FFT [Original Signal]');
xlabel('Freq (Hz)');
ylabel('Magnitude');
xlim([-100 100]);
grid on;

subplot(3,1,2);
plot(freqAxis, abs(fftshift(fft(sampledZeros))), 'r'); % FFT of Sampled Signal in red
title('FFT [Sampled Signal]');
xlabel('Freq (Hz)');
ylabel('Magnitude');
xlim([-100 100]);
grid on;

subplot(3,1,3);
plot(freqAxis, abs(fftshift(fft(reconstructedSignal)))/numPoints, 'r'); % FFT of Reconstructed Signal in red
title('FFT [Reconstructed Signal]');
xlabel('Freq (Hz)');
ylabel('Magnitude');
xlim([-100 100]);
grid on;



function [sampledSignal, sampledSignalWithZeros, sampledTime] = Sampler(inputSignal, inputTime, samplingFrequency)
    % Determine total number of samples and duration
    duration = inputTime(end) - inputTime(1);
    sampleInterval = ceil(length(inputSignal) / ceil(samplingFrequency * duration));
    buffer = zeros(1, length(inputSignal));
    % Populate buffer with sampled values
    for i = 1:sampleInterval:length(inputSignal)
        buffer(i) = inputSignal(i);
    end
    % Extract valid samples and time vectors
    validIndices = (buffer ~= 0);
    sampledSignal = buffer(validIndices);
    sampledSignalWithZeros = buffer;
    sampledTime = inputTime(validIndices);
end


function [quantizedSignal, quantizedLevels] = Quantizer(signal, numLevels, peakLevel, mu)
% Default mu-law parameter if not provided
    if nargin == 3
        mu = 0; 
    end
    numBits = ceil(log2(numLevels)); 
    deltaV = (2 * peakLevel) / numLevels; % Step size for quantization
    representationLevels = -peakLevel + deltaV / 2 : deltaV : peakLevel; % Levels for quantization
    decisionLevels = -peakLevel : deltaV : peakLevel; % Decision boundaries for quantization
    quantizedSignal = zeros(1, length(signal));
    quantizedLevels = [];
    % Apply mu-law compression if 'mu' is specified
    if mu ~= 0
        signal = peakLevel * log(1 + mu * (abs(signal) / peakLevel)) / log(1 + mu) .* sign(signal);
    end
    % Quantize each sample in the signal
    for i = 1:length(signal)
        for level = 1:length(decisionLevels)
            if signal(i) >= decisionLevels(level) && signal(i) < decisionLevels(level + 1)
                quantizedSignal(i) = representationLevels(level);
                quantizedLevels(i, :) = de2bi(level - 1, numBits, 'left-msb');
                break;
            elseif signal(i) >= decisionLevels(end)
                quantizedSignal(i) = representationLevels(end);
                quantizedLevels(i, :) = de2bi(length(representationLevels) - 1, numBits, 'left-msb');
                break;
            elseif signal(i) <= decisionLevels(1)
                quantizedSignal(i) = representationLevels(1);
                quantizedLevels(i, :) = de2bi(0, numBits, 'left-msb');
                break;
            end 
        end
    end
end

function decodedSignal = Decoder(encodedSignal, encodingType, numLevels, maxLevel, mu)
    if nargin == 4
        mu = 0; % Default mu-law parameter if not provided
    end
    bitsPerSymbol = ceil(log2(numLevels));
    numBits = length(encodedSignal);
    switch encodingType
        case 1
            % Decode Unipolar NRZ encoding
            decodedBuffer = reshape(encodedSignal, bitsPerSymbol, []).';
        case 2
            % Decode Polar NRZ encoding
            encodedSignal(encodedSignal == -1) = 0;
            decodedBuffer = reshape(encodedSignal, bitsPerSymbol, []).';           
        case 3
            % Decode Manchester encoding
            reshapedSignal = reshape(encodedSignal, 2, []).';
            buffer = reshapedSignal(:, 2) == 0; % Convert 0 -> 1 and 1 -> 0
            decodedBuffer = reshape(buffer, bitsPerSymbol, []).';           
        otherwise
            error('Invalid encoding type specified');
    end
    % Dequantize the values
    deltaV = (2 * maxLevel) / numLevels;
    representationLevels = linspace(-maxLevel + deltaV / 2, maxLevel - deltaV / 2, numLevels);
    decodedSignal = zeros(1, size(decodedBuffer, 1));
    for i = 1:size(decodedBuffer, 1)
        levelIndex = bi2de(decodedBuffer(i, :), 'left-msb') + 1;
        decodedSignal(i) = representationLevels(levelIndex);
    end
    % Apply mu-law expansion if mu is provided
    if mu ~= 0
        decodedSignal = maxLevel / mu * (exp(abs(decodedSignal) * log(1 + mu) / maxLevel) - 1) .* sign(decodedSignal);
    end
end
function encodedSignal = Encoder(levels, encodingType)
    % Flatten levels into a binary stream
    binaryStream = reshape(levels', 1, []);
    switch encodingType
        case 1
            % Unipolar NRZ encoding
            encodedSignal = binaryStream;
        case 2
            % Polar NRZ encoding
            encodedSignal = binaryStream;
            encodedSignal(encodedSignal == 0) = -1;
        case 3
            % Manchester encoding
            encodedSignal = zeros(1, 2 * length(binaryStream)); % Preallocate space
            for i = 1:length(binaryStream)
                if binaryStream(i) == 0
                    encodedSignal(2*i-1:2*i) = [0 1];
                else
                    encodedSignal(2*i-1:2*i) = [1 0];
                end
            end
        otherwise
            error('Invalid encoding type specified');
    end
end
function reconstructedSignal = Filter(decodedSignal, numSamples, samplingFreq, messageFreq)
    % Initialize a zero array for the signal
    signalArray = zeros(1, numSamples);
    % Calculate the interval for sampling
    samplingInterval = ceil(numSamples / ceil(samplingFreq));
    currentIndex = 1;
    % Populate the signal array with decoded values at specified intervals
    for position = 1:samplingInterval:numSamples
        signalArray(position) = decodedSignal(currentIndex);
        currentIndex = currentIndex + 1;
    end
    % Design the frequency filter
    filterSize = 2 * messageFreq + 1;
    padding = round(0.5 * (numSamples - filterSize));
    frequencyFilter = [zeros(1, padding), ones(1, filterSize), zeros(1, numSamples - filterSize - padding)];
    % Apply the frequency filter in the frequency domain
    signalSpectrum = fft(signalArray) * numSamples / samplingFreq;
    filteredSpectrum = fftshift(frequencyFilter) .* fftshift(signalSpectrum);
    % Perform the inverse FFT to reconstruct the signal
    reconstructedSignal = ifft(filteredSpectrum);
end