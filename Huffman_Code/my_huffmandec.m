%% Huffman Symbol Decoding
% Author: Tilemachos S. Doganis
%
%% Huffman Symbol Decoding
% Author: Tilemachos S. Doganis
%
function [ dec ] = my_huffmandec( sig, dict )
%
% Input: Source Sequence, Dictionary
% Output: Decoded symbols (cell)
%
% 'minv' contains the minimum word length in the dictonary
% 1. For each iteration, the first 'minv' bits are stored in a temporary
% search variable 'comp' which is then compared to the dictionary entries.
% If no valid key is found, the next bit is added to 'comp' and the search
% is repeated until a valid dictionary key has been found.
% 2. Once a match has been found, the corresponding dictionary entry is
% added in the decoded string 'dec'.
%
% Potential Improvement: A key-wise sorting of the huffman dictionary so
% that only keys with the same length as 'comp' are checked.
% 
    %% Initialization
    N = size(sig,2);            % Encoded sequence length
    M = size(dict,1);           % Number of dictionary entries
    t = cellfun('length',dict); % Dictionary key lengths
    minv = min(t(:,2));         % Minimum key length
    maxv = max(t(:,2));         % Maximum key length           
    
    %% Decoding
    % Convert input to character string
    sig = char(sig); 
    
    % Initialize decoding string
    dec = '';
    i = 1;
    while i <= N
        flag = 0;
        k = minv;
        while (flag == 0) && (k <= maxv) 
            % Substring to decode
            comp = sig(i:(i+k-1));
            
            % Check Dictionary (!Brute force version!Very time consuming!)
            for j = 1:M
                if (strcmp(comp,dict(j,2)))
                    dec = strcat(dec,dict(j,1));
                    flag = 1;
                end
            end
            % Check for longer string
            k = k+1;
        end
        i = i+k-1;    
    end
end

