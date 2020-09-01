function [signal,noise] = data_reuse(signal,noise,copies)    

signal_oneperiod = signal;
noise_oneperiod = noise;
signal_original_length = length(signal);
noise_original_length = length(noise);

for k=1:copies
    
    signal_updated_length = length(signal);
    noise_updated_length = length(noise);
    signal = [signal, zeros(1,signal_original_length)];
    signal = signal + [zeros(1,signal_updated_length), signal_oneperiod];
    noise = [noise, zeros(1,noise_original_length)];
    noise = noise + [zeros(1,noise_updated_length), noise_oneperiod];

end

end