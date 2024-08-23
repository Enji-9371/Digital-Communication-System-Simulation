%% AM Modulation
% Load an audio signal as the message signal
[audioSignal, Fs] = audioread('C:\Users\ATD\Music\sos-signal-137144.mp3'); %%%%%%%%%%%%%% important note: Load the audio file attached with your own path
audioSignal = audioSignal(:, 1); % Use only one channel if stereo

% Resample the audio signal to match the time vector if necessary
t = linspace(0, length(audioSignal)/Fs, length(audioSignal));

% Plot the original audio signal
figure;
subplot(6, 2, 1);
plot(t, audioSignal);
title('Original Audio Signal (Time Domain)');
xlabel('time(t)');
ylabel('Amplitude');

% Frequency domain of Audio Signal
Y_m = fftshift(fft(audioSignal));
f = linspace(-Fs/2, Fs/2, length(Y_m));
subplot(6, 2, 2);
plot(f, abs(Y_m)/max(abs(Y_m)));
title('Original Audio Signal (Frequency Domain)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Carrier Signal
h = 60;
Am = max(abs(audioSignal)); % Adjust the carrier amplitude according to the audio signal
Ac = Am / h;
fc = 20000; % Carrier frequency for modulation (choose higher than the highest frequency in the audio signal)
yc = Ac * cos(2 * pi * fc * t);
subplot(6, 2, 3);
plot(t(1:10000), yc(1:10000));
title('Carrier Signal (Time Domain)');
xlabel('time(t)');
ylabel('Amplitude');

% Frequency domain of Carrier Signal
Y_c = fftshift(fft(yc));
subplot(6, 2, 4);
plot(f, abs(Y_c)/max(abs(Y_c)));
title('Carrier Signal (Frequency Domain)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Modulate the audio signal using AM
y = ammod(audioSignal, fc, Fs, 0, Ac);
subplot(6, 2, 5);
plot(t(1:10000), y(1:10000));
title('AM Modulated Signal (Time Domain)');
xlabel('time(t)');
ylabel('Amplitude');

% Frequency domain of Modulated Signal
Y_mod = fftshift(fft(y));
subplot(6, 2, 6);
plot(f, abs(Y_mod)/max(abs(Y_mod)));
title('AM Modulated Signal (Frequency Domain)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Add AWGN with different SNR values
SNRs = [3000, 200, 10]; % Array of different SNR values in dB
for i = 1:length(SNRs)
    % Add noise
    y_awgn = awgn(y, SNRs(i), 'measured');
    
    % Plot Time Domain of Noisy Signal
    subplot(6, 2, 2*i+5);
    plot(t(1:10000), y_awgn(1:10000));
    title(['AM Modulated Signal with AWGN (SNR = ' num2str(SNRs(i)) ' dB) Time Domain']);
    xlabel('time(t)');
    ylabel('Amplitude');
    
    % Plot Frequency Domain of Noisy Signal
    Y_awgn = fftshift(fft(y_awgn));
    subplot(6, 2, 2*i+6);
    plot(f, abs(Y_awgn)/max(abs(Y_awgn)));
    title(['Noisy AM Modulated Signal with AWGN (SNR = ' num2str(SNRs(i)) ' dB) Frequency Domain']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
end

% Demodulate and play back the audio signals
figure;
for i = 1:length(SNRs)
    % Demodulate
    z_awgn = amdemod(y_awgn, fc, Fs, 0, Ac);
    
    % Plot Time Domain of Demodulated Signal
    subplot(length(SNRs), 2, 2*i-1);
    plot(t(1:10000), z_awgn(1:10000));
    title(['Demodulated Audio Signal with AWGN (SNR = ' num2str(SNRs(i)) ' dB) Time Domain']);
    xlabel('time(t)');
    ylabel('Amplitude');
    
    % Plot Frequency Domain of Demodulated Signal
    Z_awgn = fftshift(fft(z_awgn));
    subplot(length(SNRs), 2, 2*i);
    plot(f, abs(Z_awgn)/max(abs(Z_awgn)));
    title(['Demodulated Audio Signal with AWGN (SNR = ' num2str(SNRs(i)) ' dB) Frequency Domain']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    
    % Play back the demodulated noisy signal
    sound(z_awgn, Fs);
    pause(5); % Pause for 5 seconds to allow playback before next sound
end
