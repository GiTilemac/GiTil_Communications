%% 4-PSK
%
% Author: Tilemachos S. Doganis
%
%   This system transmits a binary message using the Modulation type 4-PSK
% Phase Shift Keying using 4 Symbols.
%
% 	Initially, the message is produced as a random text containing the 
% first four lower case letters of the English alphabet. The text is then
% converted into bits using Gray coding. 
%
% The transmission system consists of four components:
% 1) Modulator: It produces the transmitted waveform by shifting its phase
% component.
%
% Ps_m(t) = g_T(t) * cos(2*pi*f_c*t + 2*pi*m/M)
% <=> g_T(t) * A_mc * cos(2*pi*f_c*t) - g_T(t) * A_ms * sin(2*pi*f_c*t)
% where A_mc = cos(2*pi*m/M) & A_ms = sin(2*pi*m/M) represent the phase
% shift
%
% The signal is therefore expressed as the sum of two carrier signals, 
% orthogonal to each other, modulated by two baseband signals.
%
% 2) Channel: During transmission, AWGN is added to the transmitted signal.
% 5 different instances of transmission are simulated, with increasing
% SNR. 
% 
% 3) Demodulator: Here, the correlation between the received noisy signal
% and the two carrier signals is calculated.
% 
% 4) Detector: The modulation the smallest norm-2 distance from the
% corresponding correlation pair is selected.
%
% Using the Gray coding dictionary, the estimated bit pair is converted
% back into a symbol.

%% Source Signal Initialization
close all
clearvars -except Fs_m
clc

L_b = 10000;  				% Length in bits
symb = {'a' 'b' 'c' 'd'};   % Available symbols

% Gray coding
cod = {'00' '01' '11' '10'};
dict = {symb cod};

% Create random initial text
inp = char(randi([97 100],1,L_b/2));

% Encode text using Gray Coding dictionary
for i = 1:4
    inp = strrep(inp,symb{i},cod{i});
end

%% Transmission parameters
E_s = 1;        % Energy per symbol
E_b = E_s/2;    % Energy per bit
T_sym = 40;     % Symbol period
f_c = 250000;   % Carrier frequency
T_c = 4;        % Carrier period
g_T = sqrt(2/T_sym);      % Baseband (square) pulse coefficient
t = linspace(0,T_sym,40); % 0 <= t <= T_s

% Store the symbol phase coefficients
s = zeros(2,4);
for m = 1:4
    s(:,m) = [ cos(2*pi*(m-1)/4) ; sin(2*pi*(m-1)/4) ]; 
end

siz = T_sym*L_b/2;
Ps_m = zeros(1,siz);
r_t = zeros(1,siz);

SNR = 0:2:8;
BER = zeros(5,1);
P_b = zeros(5,1);

for j = 1:5
    %% Specific SNR
    % Produce AWGN
    var = E_b / (2 * 10^(SNR(j)/10));
    noise = sqrt(var)*randn(1,L_b/2*40); 

    % output buffer
    out = blanks(1);

    % bit pair to send
    buf = blanks(2);

    %% Transmission loop
    for i = 1:2:L_b-1
       
       sindex = 1+(ceil(i/2)-1)*40:ceil(i/2)*40; % Current symbol indices
        
       % Mapper (Bit pair -> Symbol coefficient)
       buf = strcat(inp(i),inp(i+1));
       m = find(strcmp(cod,buf));

       % Modulator (Phase shift)
       Ps_m(sindex) = s(1,m).*g_T.*cos(2*pi*f_c*t) + s(2,m).*g_T.*sin(2*pi*f_c*t);
       
    end
    for i = 1:2:L_b-1
        
       sindex = 1+(ceil(i/2)-1)*40:ceil(i/2)*40; % Current symbol indices
        
       % Channel (Adds AWGN)
       r_t(sindex) = Ps_m(sindex) + noise(sindex);

       % Demodulator (Correlation)
       r1 = sum(r_t(sindex).*g_T.*cos(2*pi*f_c*t));
       r2 = sum(r_t(sindex).*g_T.*sin(2*pi*f_c*t));
       r = [r1 ; r2];

       % Detector
       D = zeros(4,1);
       for m = 1:4
          D(m) = norm(r-s(:,m),2);
       end
       [~,f] = min(D);

       % Demapper (Received symbol coefficient -> bit pair)
       out = strcat(out,cod{f});  
    end
     %% Results
    BER(j) = 1-sum(out == inp)/L_b; 
    P_b(j) = 1/2 * erfc(sqrt(E_b/(2*var)));
    fprintf('4-PSK using %d samples (SNR = %d)\n',L_b,SNR(j));
    fprintf('BER = %f, P_b = %f\n',BER(j),P_b(j));
end