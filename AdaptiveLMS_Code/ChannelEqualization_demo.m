%% Channel Equalization using an Adaptive LMS Equalizer
%
% Author: Tilemachos S. Doganis
%
% 1. Packets of random symbols are produced and coded using 4-PAM
% 2. The packets are transmitted through an emulated channel with ISI and
% AWGN.
% 3. After receiving the noisy symbol sequence, it is equalized using an
% Adaptive LMS Equalizer, and then decoded.
% 4. Lastly, statistics on the Bit Error Rate (BER) and the LMS learning
% curves for the Equalizer are taken.
%
% This script contains 5 versions of the above procedure for comparison
% purposes:
% Channel with / without ISI
% Equalization / no Equalization after the Maximum Likelihood Detector
% Time-Variant Channel with ISI, AWGN, and Equalization

close all
clear
clc

%% Initialization

% M-PAM characteristics
M = 4;
Am = 2*(1:M)'-1-M; % M-PAM dictionary
k = log2(M);       % Bits used
Eg = 1;            % Waveform energy
Es = Eg*(M^2-1)/3; % Mean Symbol energy
Eb = Es/k;         % Mean bit energy

% Definitions
chnl = cell(2,1);      % Channels
chnl{1} = 1; % Channel 1 (no ISI)
chnl{2} = [0.04 -0.05 0.07 -0.21 -0.5...
           0.72 ...
           0.36   0   0.21  0.03  0.07]'; %Channel 2 (with ISI)
SNR = 0:10;     % SNRs (dB)
Q = 500;        % First Q symbols of packet are training sequence
P = 10500;      % Whole packet contains P symbols
N_pack = 50;   % Number of packets to transmit

% Memory allocation
eql = cell(11,3);         % Equalizers per SNR per channel
P_M = zeros(11,1);        % Theoretical BERs per SNR
BER_ML = zeros(11,2);     % Experimental ML BERs per SNR
BER_LMS = zeros(11,2);    % Experimental LMS BERs per SNR
BER_tv_LMS = zeros(11,1); % Experimental LMS BERs  (time variant channel)
MSE_1 = zeros(P,2);       % MSE curves for the average packet
MSE_2 = zeros(P,2); 
MSE_3 = zeros(P,2); 
MSE_tv = zeros(P,3);      % MSE curves (time variant channel)

% Initialize Equalizers
[eql{:, 1}] = deal(1);
[eql{:, 2:3}] = deal([ zeros(5,1) ; 1 ; zeros(5,1)]);

%% Transmission and Equalization loops

for j = 1:N_pack

    if mod(j,10) == 1
        fprintf(['Transmitting packets ' num2str(j) ' to ' num2str(j+9) '...\n']);
    end

    % Generate learning sequence
    lrn = randi(4,Q,1);
    lrn(lrn==2) = -1;
    lrn(lrn==4) = -3;

    % Generate packet bit sequence
    bit_packet = randi(2,(P-Q)*k,1)-1;

    % Gray encoding for bit error minimization
    sym_packet = enc_4pam(bit_packet);

    % Final packet
    a = [lrn; sym_packet];
    
    for c_i = 1:2

        h = chnl{c_i};   % Channel response
        L = numel(h);    % Channel length
        N = floor(L/2);  % Delay

        %% Transmit and equalize sequence packets
        for i = 1:11
            %% Packet transmission
            % Noise parameters
            N0 = Eb/(10^(SNR(i)/10));
            sigma_w = Eg*N0/2;
            
            % Theoretical BER
            P_M(i) = (2*(M-1)/M) * qfunc(sqrt( 6*log2(M)/(M^2-1)* Eb/N0  ));
            
            % Generate AWGN
            w_n = sqrt(sigma_w)*randn(P,1);
            
            % Apply ISI and AWGN
            y = my_isi(a,h) + w_n;  

            %% Task 1 - ML Detector
            z = zeros(P-Q,1);
            for l = Q+1:P
                % Max likelihood = Min Euclidean distance
                [~,dcs] = min(abs(Am - y(l))); 
                z(l-Q) = Am(dcs);
            end
            
            % Calculate BER
            BER_temp1 = sum(dec_4pam(z) ~= bit_packet)/(2*(P-Q));

            %% Task 3 - Adaptive LMS Equalizer
            [~, rcv_sym, eql{i,c_i}] = lms_eq(y, a, Am, Q, P, L, eql{i,c_i}, 1);
            
            % Calculate BER
            BER_temp2 = sum(dec_4pam(rcv_sym) ~= bit_packet)/(2*(P-Q));
            
            % SNR-wide mean BER
            BER_ML(i,c_i) = BER_ML(i,c_i) + BER_temp1;
            BER_LMS(i,c_i) = BER_LMS(i,c_i) + BER_temp2;  
            
            %% Time Variant Channel (Task 6)
            if c_i == 2
                % Apply Time-Variant ISI and AWGN
                y = my_tv_isi(a,h,(P+Q)/2) + w_n;

                % Equalize using LMS
                [~, rcv_sym, eql{i,3}] = lms_eq(y, a, Am, Q, P, L, eql{i,3}, 1);

                % Calculate BER
                BER_temp3 = sum(dec_4pam(rcv_sym) ~= bit_packet)/(2*(P-Q));
                BER_tv_LMS(i) = BER_tv_LMS(i) + BER_temp3;
            end
 
        end
        
        %% Task 4 - Packet transmission (SNR = 20 dB)
        % Noise parameters
        N0 = Eb/(10^2); %SNR = 20dB
        sigma_w = N0/2;

        % Generate AWGN
        w_n = sqrt(sigma_w)*randn(P,1);

        % Apply ISI and AWGN
        y = my_isi(a,h) + w_n;
        
        % Initialize Equalizers
        eql_init = [ zeros(N,1) ; 1 ; zeros(N,1)];
        
        % Adaptive LMS Equalizer
        [MSE_temp1, ~, ~] = lms_eq(y, a, Am, Q, P, L, eql_init, 1);
        [MSE_temp2, ~, ~] = lms_eq(y, a, Am, Q, P, L, eql_init, 0.01);
        [MSE_temp3, ~, ~] = lms_eq(y, a, Am, Q, P, L, eql_init, 2.5);

        % Cumulative packet sum
        MSE_1(:,c_i) = MSE_1(:,c_i) + MSE_temp1';
        MSE_2(:,c_i) = MSE_2(:,c_i) + MSE_temp2';
        MSE_3(:,c_i) = MSE_3(:,c_i) + MSE_temp3';
        
        %% Time Variant Channel (Task 6)
        if c_i == 2
            % Apply Time-Variant ISI and AWGN
            y_tv = my_tv_isi(a,h,(P+Q)/2) + w_n;
            
            % Adaptive LMS Equalizer
            [MSE_temp11, ~, ~] = lms_eq(y_tv, a, Am, Q, P, L, eql_init, 1);
            [MSE_temp22, ~, ~] = lms_eq(y_tv, a, Am, Q, P, L, eql_init, 0.01);
            [MSE_temp33, ~, ~] = lms_eq(y_tv, a, Am, Q, P, L, eql_init, 2.5);

            % Cumulative packet sum
            MSE_tv(:,1) = MSE_tv(:,1) + MSE_temp11';
            MSE_tv(:,2) = MSE_tv(:,2) + MSE_temp22';
            MSE_tv(:,3) = MSE_tv(:,3) + MSE_temp33';
        end
    end
end

% Packet-wide mean BER
BER_ML = BER_ML/N_pack;
BER_LMS = BER_LMS/N_pack;
BER_tv_LMS = BER_tv_LMS/N_pack;
        
% Packet-wide mean MSE
MSE_1 = MSE_1/N_pack;
MSE_2 = MSE_2/N_pack;
MSE_3 = MSE_3/N_pack;
MSE_tv = MSE_tv/N_pack;

for c_i = 1:2
    %% Plot BER
    %P_M - ML
    figure(c_i),
    subplot(2,1,1),
    hold on,
    semilogy(SNR,BER_ML(:,c_i),'b'),
    semilogy(SNR,P_M,'r'),
    hold off,
    title(['Channel ' num2str(c_i) ]),
    grid on,
    xlabel('SNR|_{dB}'),
    ylabel('log_{10}(BER)'),
    legend('ML', 'Theory');

    %P_M - LMS
    subplot(2,1,2),
    hold on,
    semilogy(SNR,BER_LMS(:,c_i),'b'),
    semilogy(SNR,P_M,'r'),
    hold off,
    title(['Channel ' num2str(c_i)]),
    grid on,
    xlabel('SNR|_{dB}'),
    ylabel('log_{10}(BER)'),
    legend('LMS', 'Theory');
    saveas(gcf, ['ch' num2str(c_i) '_BER.png']);
    
    %% Plot MSE     
    figure(2+c_i)
    hold on,
    plot(MSE_1(:,c_i),'b'),
    plot(MSE_2(:,c_i),'r'),
    plot(MSE_3(:,c_i),'g'),
    title(['Equalizer Learning Curves (Channel ' num2str(c_i) ')']),
    xlabel('Symbols'),
    ylabel('MSE'),
    xlim([0 Inf]),
    grid on,
    legend('\Delta', '0.01\Delta', '2.5\Delta');
    saveas(gcf, ['ch' num2str(c_i) '_train.png']);
    
    %% Task 5 - Plot Equalized Channel
    cnv = conv(eql{11,c_i},chnl{c_i},'full');
    L = numel(chnl{c_i});
    N = floor(L/2);
    cnv = cnv(N+1:3*N+1);
    figure(4+c_i),

    if c_i == 1
        hold on,
        stem(chnl{c_i}),
        stem(eql{11,c_i},'g'),
        stem(cnv,'r'),
        set(gca,'XTIck',1,'XTickLabel','n'),
        hold off
    else
        hold on,
        plot(chnl{c_i}),
        plot(eql{11,c_i},'g'),
        plot(cnv,'r'),
        nlabel = {'n-5', 'n-4', 'n-3', 'n-2', 'n-1' ,'n' ,'n+1', 'n+2', 'n+3', 'n+4', 'n+5'};
        set(gca,'XTick',1:L,'XTickLabel',nlabel),
        hold off,
    end
    title(['Channel ' num2str(c_i) ' Equalization']),
    grid on,
    legend('Original Channel','Equalizer','Equalized Channel');
    saveas(gcf, ['ch' num2str(c_i) '_EQ.png']);
end

%% Plot LMS BER (Time Variant Channel - Task 6)
figure(7)
hold on
semilogy(SNR,BER_LMS(:,2),'b'),
semilogy(SNR,BER_tv_LMS,'r'),
title('LMS (Channel 2)'),
grid on,
xlabel('SNR|_{dB}'),
ylabel('log_{10}(BER)'),
legend('Time Invariant', 'Time Variant');
saveas(gcf, 'ch_tv_LMS_BER.png');

%% Plot MSE (Time Variant Channel - Task 6)     
figure(8)
hold on,
plot(MSE_tv(:,1),'b'),
plot(MSE_tv(:,2),'r'),
plot(MSE_tv(:,3),'g'),
title('Equalizer Learning Curves (Time Variant Channel)'),
xlabel('Symbols'),
ylabel('MSE'),
xlim([0 Inf]),
grid on,
legend('\Delta', '0.01\Delta', '2.5\Delta');
saveas(gcf, 'ch_tv_train.png');
