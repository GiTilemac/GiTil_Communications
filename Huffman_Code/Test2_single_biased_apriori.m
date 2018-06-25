%% Huffman Coding Test 2 (Biased Source, A Priori Symbol Probabilities)
%
% Author: Tilemachos S. Doganis
%
% This test uses a set of words beginning with the letter 'k' as a source,
% and the a priori letter probabilities for the English alphabet.
%
% Due to the fact that the input contained only English words beginning
% with a specific letter, the a priori probabilities are not valid and
% therefore the average coding length is increased by 0.5 bits/symbol,
% compared to source A.
%
clc
clearvars -except txtA
fprintf('Running Huffman Encoding of Source B with Given Probabilities... \n');

% Load kwords.txt in matrix T
fid = fopen('kwords.txt','r');
T = textscan(fid,'%s');
fclose(fid);
T = T{1};

% Store list T as string txtB
txtB = T{1};
tsize = numel(T);
for i=2:tsize
   txtB = strcat(txtB,T(i)); 
end
txtB = char(txtB);


% Remove special characters
for i = 33:96
    remPos = strfind(txtB,char(i));
    remLen = numel(remPos);
    for j = remLen:-1:1
        txtB(remPos(j))='';
    end
end
total_chars = numel(txtB);

%Ensure that length is even
if ( mod( numel(txtB), 2 ) > 0)
   txtB = strcat(txtB, 'k'); 
end
    
% Probabilities for lower case english letters
probB = [ 8.167 1.492 2.782 4.253 12.702 2.228 2.015 6.094 6.966 0.153 ...
          0.772 4.025 2.406 6.749 7.507 1.929 0.095 5.987 6.327 9.056 ...
          2.758 0.978 2.361 0.150 1.974 0.074 ];
% Normalization
probB = probB/100;

% Alphabet is lower case English alphabet
alphB = 'a':'z';

% Construct Huffman Dictionary
dictB = my_huffmandict(probB,alphB);

% Calculate Theoretical Average Code Length
c = size(dictB,1);
L_t = 0;
for i = 1:c
    L_t = L_t + (numel(dictB{i,2})*probB(i));
end

% Encode Text using Huffman Encoding
encB = my_huffmanenc(txtB,dictB);
encB = encB{1};
encLenB = numel(encB);

% Decode Encoding
decB = my_huffmandec(encB,dictB);
decB = decB{1};

% Check if encoding / decoding are successful
if (strcmp(txtB,decB) == 1)
    disp('Encoding & Decoding Successful');
else
    disp('Encoding / Decoding Error');
end

L = encLenB/total_chars;
fprintf('Text length: %i \n',total_chars);
fprintf('Encoding length: %i \n',encLenB);
fprintf('Average Code Length: L = %.4f bits/symbol (Theory: %.4f)\n',L,L_t);

H = -sum(probB.*log2(probB));
fprintf('Source B Entropy: H(B) = %f\n',H);
fprintf('Source B Efficiency: n = %.2f%%\n',H/L*100);