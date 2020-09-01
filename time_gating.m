function gated_signal = time_gating(signal)

%% Step 1: creating a customized window

% Defining the envelope of the signal
y = hilbert(signal); % Hilbert Transform
env = abs(y);        % Envelope function
peak = max(env);     % Peak value of the envelope

% Initializing variables
window = zeros(length(signal),1); % Window function
Lwin = length(window);            % Window length
t1 = 0;                           % Start of the rise time
t2 = 0;                           % End of the rise time
t3 = 0;                           % Start of the fall time
t4 = 0;                           % End of the fall time
            
% Identifying t1-t4
for k=1:Lwin    
    if (t1 == 0 && env(k) >= 0.05*peak)
        t1 = k;
    elseif (t1 ~= 0) && (t2 == 0) && (env(k) >= 0.40*peak)
        t2 = k; 
    elseif(t2 ~= 0) && (t3 == 0) && (env(k) < 0.40*peak)
        t3 = k; 
    elseif (t3 ~= 0) && (t4 == 0) && (env(k) < 0.05*peak)
        t4 = k;
        break
    end    
end

% Rise and fall times
Lup = t2 - t1;             % Rise time
Ldown = t4 - t3;           % Fall time

% Window function used during rise time
fup = blackman(2*Lup)';    
fup = fup(1:Lup);

% Window function used during fall time
fdown = hamming(2*Ldown)'; 
fdown = fdown(end-Ldown+1:end);

% Creating the customized window function
window(t2+1:t3) = 1;                     % Rectangular portion
window((t2 - Lup + 1):t2,1) = fup;       % Blackman portion
window((t3 + 1):(t3 + Ldown),1) = fdown; % Hamming portion

%% Step 2: applying the window

% Gated signal
gated_signal = window.*signal;                                 % Time domain
gated_signal_spectrum = 10*log10(abs(fft(gated_signal)));      % Frequency domain

% Graphs
figure
hold on
plot(signal./max(signal))
plot(env./max(env))
plot(window)
legend('Signal','Envelope','Window')
ylabel('Normalized amplitude')
xlabel('Time')
hold off

figure
plot(gated_signal)
legend('Gated signal')
ylabel('Amplitude')
xlabel('Time')

% figure
% hold on
% plot(freqHz,10*log10(abs(s12_camara)))
% plot(freqHz,10*log10(abs(s12_lea)))
% plot(freqHz,gated_signal_spectrum(1:800))
% hold off
% legend('Chamber','Signal','Gated signal')

end