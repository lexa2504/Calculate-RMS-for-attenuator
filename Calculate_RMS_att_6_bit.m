%TEST_CATT_0dBm_dB_0.5

clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Const
number_bits = 6; 
min_freq    = 4e9;
max_freq    = 10e9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read S-parameters

numfiles    = pow2(number_bits)-1;
mag_step  = 0.5;

for i = 1:(numfiles)
    num(i) = i*mag_step;
end
filename    = "PS_test__"+(num)+".s2p"; % Construct filenames
S           = sparameters(filename(1)); % Read file #1 for initial set-up
freq        = S.Frequencies; % Frequency values are the same for all files
numfreq     = numel(freq); % Number of frequency points
s21_data    = zeros(numfreq,numfiles); % Preallocate for speed
idx         = (freq >= min_freq) & (freq <= max_freq);

%% Read Touchstone files
for n   = 1:(numfiles)
    S   = sparameters(filename(n));
    s21 = rfparam(S,2,1);
    s21_data(:,n) = s21;
    s11 = rfparam(S,1,1);
    s11_data(:,n) = s11;
    s22 = rfparam(S,2,2);
    s22_data(:,n) = s22;
    s12 = rfparam(S,1,2);
    s12_data(:,n) = s12;
end
%% Work with data S21 db
s21_db = 20*log10(abs(s21_data));
s21_degree =(angle(s21_data));
s11_db = 20*log10(abs(s11_data));
s22_db = 20*log10(abs(s22_data));
s12_db = 20*log10(abs(s12_data));
%

%% Work with data S21 degree
s21_pass_db   = s21_db(idx,:);
s21_pass_degree = (unwrap(s21_degree(idx,:))*180/pi);
freq_pass_ghz   = freq(idx)/1e9;    % Normalize to GHz  
%%
% %Graph S21_db, s11_db, s22_db, S12_db
% 
t = tiledlayout(2,2);
nexttile
plot(freq/1e9,s21_db)
xlabel('Frequency (GHz)')
ylabel('S_2_1 (dB)')
grid on
nexttile
plot(freq/1e9,s11_db)
xlabel('Frequency (GHz)')
ylabel('S_1_1 (dB)')
grid on
nexttile
plot(freq/1e9,s22_db)
xlabel('Frequency (GHz)')
ylabel('S_2_2 (dB)')
grid on
nexttile
plot(freq/1e9,s12_db)
xlabel('Frequency (GHz)')
ylabel('S_1_2 (dB)')
grid on

%%

for i = 1:numfiles
            S21_db_norm(:,i) = s21_pass_db(:,1) - s21_pass_db(:,i);
end
%Graph S21_degree
figure
plot(freq_pass_ghz,S21_db_norm)
xlabel('Frequency (GHz)')
ylabel('S_2_1 (dB)')
%title('S21 в градусах для 64 состояний')
axis on
grid on

numb = length(freq_pass_ghz);
%ideal_mag_step
for j = 1:numfiles
    for k = 1:numb   %найти потом как 200 чтоб автоматом брало
    mag_step_ideal(k,j) = (j-1)*mag_step;
    end
end

%Graph S21_db_norm
figure
plot(freq_pass_ghz,mag_step_ideal)
xlabel('Frequency (GHz)')
ylabel('S_2_1 (dB)')
%title('S21 в градусах для 64 состояний')
axis on


% RMS mag CALC
mag_error=S21_db_norm - mag_step_ideal;

for i=1:numb
    m_mag_error(i,:) =mean(mag_error(i,:));
end

for i=1:(numfiles)
    m_mag_error_all(:,i) = m_mag_error(:,1);
end

mag_error_abs = mag_error - m_mag_error;

for i=1:numfiles
      mag_error_abs_mul(:,i)= mag_error_abs(:,i).*mag_error_abs(:,i);
end

for i=1:k
      mag_error_abs_summ(i,:)= sum(mag_error_abs_mul(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMS_mag_deviation = sqrt(mag_error_abs_summ/(numfiles-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graph S21_db_norm
figure
plot(freq_pass_ghz,RMS_mag_deviation)
xlabel('Frequency (GHz)')
ylabel('S_2_1 (dB)')
%title('S21 в градусах для 64 состояний')
axis on
%%

%%
%RMS phase CALC

for i=1:numb
    mean_phase_error(i,:) =mean(s21_pass_degree(i,:));
end

for i=1:(numfiles)
    mean_phase_error_all(:,i) = mean_phase_error(:,1);
end

phase_error = s21_pass_degree - mean_phase_error_all;

for i=1:numfiles
     phase_error_db_mul(:,i)= phase_error(:,i).*phase_error(:,i);
end

for i=1:k
      phase_error_db_summ(i,:)= sum(phase_error_db_mul(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMS_phase_deviation = sqrt(phase_error_db_summ/(numfiles-1));
%%%%%%%%%%%%
%RMS_phase_deviation   = RMS_phase_deviation(idx,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
figure
plot(freq_pass_ghz,RMS_phase_deviation)
xlabel('Frequency (GHz)')
ylabel('S_2_1 (dB)')
%title('S21 в градусах для 64 состояний')
axis on
