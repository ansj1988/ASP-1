% Antenna gain estimates from noisy measurement data
% Developed by: Artur Nogueira de São José
% Graduate Program in Electrical Engineering - UFMG, Brazil

% Technique: ASP-1
% Current status: under review by a Ph.D. jury

% If you use this code, please cite this paper:

% @article{NDESAOJOSE2020107720,
% title = "Improving antenna gain estimations in non-ideal test sites with auto-tunable filters",
% journal = "Measurement",
% volume = "159",
% pages = "107720",
% year = "2020",
% issn = "0263-2241",
% doi = "https://doi.org/10.1016/j.measurement.2020.107720",
% url = "http://www.sciencedirect.com/science/article/pii/S026322412030258X",
% author = "Artur {N. de São José} and Virginie Deniau and Úrsula {do C. Resende} and Ricardo Adriano"
% }

%% Step 1: loading variables and configuring the tool

clear
close all
clc

% Examples of input signals:
% 1. Directional patch antenna
% 2. Omnidirectional patch antenna
% 3. Horn antenna

% Details of the example files:
% S12 measurements between identical antennas
% freqHz: frequency vector
% s12_chamber: reference measurement performed on an anechoic chamber
% s12_lab and noise_floor: original lab measurements
% s12_lab_10xemi and noise_floor_10xemi: original measurements with artificially increased noise

input = 1;
if input == 1
    load('directional_patch_antenna.mat'); 
elseif input == 2
    load('omni_patch_antenna.mat');        
elseif input == 3
    load('horn_antenna') 
end

ifft_size = 5000;
signal = ifftshift(ifft(s12_lab_10xemi,ifft_size,'symmetric'))';
noise = ifftshift(ifft(noise_floor_10xemi,ifft_size,'symmetric'))';

%% Step 2: filtering the signals

% If needed, you can make signal replicas in order to increase the filter learning time
reuse = 1;                                              % Enables data reuse
if reuse == 1
    copies = 2;                                         % Number of replicas
    [signal,noise] = data_reuse(signal,noise,copies);   % Applies data reuse
end

% Obtaining the optimum parameters for the NLMS adaptive filter
[optimum_stepsize,optimum_order] = pso_main(signal,noise);

% Applying the adaptive filtering
e1 = nlms_evaluation(signal,noise,optimum_order,optimum_stepsize); % Filtered signal in time domain
e1 = e1(end-(ifft_size-1):end,1);                                  % If data reuse is enabled, only the last portion of the signal is taken

%%% FILTERED SIGNAL %%%
tg1 = time_gating(e1); % Applying time gating with automatic adjustment to the adaptive filter output

%% Graphs

af_output_freq = fft(e1);   % S12 at the adaptive filter output
tg_output_freq = fft(tg1);  % S12 at the time gating output
freq_size = length(freqHz); % Original number of frequency samples (no zero padding)

% Converting S12 into gain
d = 1;                                                                            % Distance between the antennas (in meters)
lambda = (3e8./freqHz)';                                                          % Wavelengths (in meters)
gain_lab = 10*log10( (4*pi*d./lambda).*abs(s12_lab_10xemi)' );                    % Measured signal
gain_chamber = 10*log10( (4*pi*d./lambda).*abs(s12_chamber)' );                   % Reference gain curve
gain_afoutput = 10*log10( (4*pi*d./lambda).*abs(af_output_freq(1:freq_size))' );  % Filtered signal (1st stage output)
gain_tgoutput = 10*log10( (4*pi*d./lambda).*abs(tg_output_freq(1:freq_size))' );  % Filtered signal (2nd stage output)

figure
hold on
plot(freqHz,gain_lab,'b')
plot(freqHz,gain_afoutput,'y')
plot(freqHz,gain_tgoutput,'k')
plot(freqHz,gain_chamber,'r')
ylabel('Gain (dB)')
xlabel('Frequency (Hz)')
legend('Original signal','1st stage: AF','2nd stage: AF+TG','Chamber')
hold off

% A filtering quality measure based on the correlation between the cleaned
% and reference (anechoic chamber) gain curves
correl_nofilter = corr(gain_lab',gain_chamber','Type','Pearson')
filtering_quality = corr(gain_tgoutput',gain_chamber','Type','Pearson')