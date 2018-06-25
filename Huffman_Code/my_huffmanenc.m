%% Huffman Symbol Encoding
% Author: Tilemachos S. Doganis
%
function [ enc ] = my_huffmanenc( sig, dict, varargin )
%
% Input: Source Sequence, Dictionary, Data Expansion Rank
% Output: Cell containing encoded bit sequence
%
% The symbol to be encoded is stored in 'comp'. The variable is then
% compared with every symbol in the dictionary, and the corresponding
% coding is saved in variable 'enc'. The process is repeated until every
% symbol has been encoded.
%
%% Initialization
    if nargin > 2
        L = varargin{1};
    else
        L = 1;
    end
    N = size(sig,2)-L+1;
    M = size(dict,1);
    enc = '';
    
%% Word encoding
    for i=1:L:N
        % Word to encode
        comp = sig(i:(i+L-1));

        % Search in dictionary
        for j=1:M
            if (strcmp(dict(j,1),comp))
               enc = strcat(enc,dict(j,2)); % Encode symbol
            end
        end
    end
end
