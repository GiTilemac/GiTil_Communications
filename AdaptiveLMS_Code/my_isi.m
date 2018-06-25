%% InterSymbol Interference Emulation Function
%
% Author: Tilemachos S. Doganis
%
% This function emulates the addition of ISI to a sequence 'a' passing
% through a channel with impulse response 'h'
%

function [ y ] = my_isi( a, h )
    L = numel(h);
    N = floor(L/2);
    M = numel(a);
    y = zeros(M,1);
    
    for k = 1:M
        
        % Future symbols
        L1 = min(M-k,N);
        ISI_L1 = 0;
        for l = 1:L1
            ISI_L1 = ISI_L1 + a(k+l)*h(N-l+1);
        end
        
        % Past symbols
        L2 = min(k-1,N);
        ISI_L2 = 0;
        for l = 1:L2
            ISI_L2 = ISI_L2 + a(k-l)*h(N+l+1);
        end
        
        y(k) = a(k)*h(N+1) + ISI_L1 + ISI_L2;
        
    end
end

