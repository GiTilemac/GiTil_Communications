%% 4-FSK
%
% Author: Tilemachos S. Doganis
%
%   This system transmits a binary message using the Modulation type 4-FSK
% Frequency Shift Keying using 4 Symbols.
%
% 	Initially, the message is produced as a random text containing the 
% first four lower case letters of the English alphabet. The text is then
% converted into bits using Gray coding. 
%
%   The transmission system consists of four components:
% 1) Modulator: It produces the transmitted waveform by shifting its
% carrier frequency 'f_c'
%
% Fs_m(t) = sqrt(2*E_s/T_s) * cos[2*pi*(f_c + m*Delta_f)*t]
% where Delta_f is the distance between consequent symbol frequencies.
% The minimal Delta_f for base signal orthogonality is 1/(2*T_s). To ensure
% phase continuity, a greater distance of at least 1/T needs to be
% selected.
%
% 2) Channel: During transmission, AWGN is added to the transmitted signal.
% 5 different instances of transmission are simulated, with increasing
% SNR. 
% 
% 3) Demodulator: Here, the correlation between the received noisy signal
% and all possible modulations is calculated.
% 
% 4) Detector: The modulation with the highest correlation to the received
% signal is selected as the most probable bit pair.
%
% Using the Gray coding dictionary, the estimated bit pair is converted
% back into a symbol.

%% Source Signal Initialization
close all
clearvars -except Ps_m
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

siz = T_sym*L_b/2;
Fs_m = zeros(1,siz);
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
        
       % Mapper (Bit pair -> Carrier frequency)
       buf = strcat(inp(i),inp(i+1));
       m = find(strcmp(cod,buf));

       % Modulator (Frequency shift)
       Fs_m(sindex) = g_T.*cos(2*pi*(f_c + (m-1)/T_sym).*t);

     end
     for i = 1:2:L_b-1
       
       sindex = 1+(ceil(i/2)-1)*40:ceil(i/2)*40; % Current symbol indices
         
       % Channel (Adds AWGN)
       r_t(sindex) = Fs_m(sindex) + noise(sindex);

       % Demodulator (Correlation)
       r = zeros(4,1);
       for m = 1:4
            r(m) = sum(r_t(sindex).*g_T.*cos(2*pi*(f_c + (m-1)/T_sym).*t));
       end

       % Detector
       [~,f] = max(r);

       % Demapper (Received frequency -> bit pair)
       out = strcat(out,cod{f});  
       
     end
     %% Results
     BER(j) = 1-sum(out == inp)/L_b; 
     fprintf('4-FSK using %d samples (SNR = %d)\n',L_b,SNR(j));
     fprintf('BER = %f\n',BER(j));
end
