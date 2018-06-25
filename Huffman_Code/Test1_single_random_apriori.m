%% Huffman Coding Test 1 (Random Source, A Priori Symbol Probabilities)
% Author: Tilemachos S. Doganis
%
% This test uses 10.000 random characters containing lower case English
% letters, as a source, and the a priori letter probabilities for the
% English alphabet.
%
% The average coding length is shown to be close to the lower bound for
% lossless coding, represented by the source entropy.
%
clc
clear
fprintf('Running Huffman Encoding of Source A with A Priori Probabilities... \n');

% Probabilities for lower case english letters
% (https://en.wikipedia.org/wiki/Letter_frequency)
probA = [ 8.167 1.492 2.782 4.253 12.702 2.228 2.015 6.094 6.966 0.153 ...
          0.772 4.025 2.406 6.749 7.507 1.929 0.095 5.987 6.327 9.056 ...
          2.758 0.978 2.361 0.150 1.974 0.074 ];
% Normalization
probA = probA/100;

% Alphabet is lower case English alphabet
alphA = 'a':'z';

% Produce text of 10.000 characters from random source A
inp = [double(alphA) ; probA];
txtA = char(randsrc(1,10000,inp));

% Construct Huffman Dictionary
dictA = my_huffmandict(probA,alphA);

% Calculate Theoretical Average Code Length
c = size(dictA,1);
L_t = 0;
for i = 1:c
    L_t = L_t + (numel(dictA{i,2})*probA(i));
end

% Encode Text using Huffman Coding
encA = my_huffmanenc(txtA,dictA);
encA = encA{1};

% Decode Code
decA = my_huffmandec(encA,dictA);
decA = decA{1};

% Check if encoding / decoding are successful
if (strcmp(txtA,decA) == 1)
    disp('Encoding & Decoding Successful');
else
    disp('Encoding / Decoding Error');
end

L = numel(encA)/numel(txtA);
fprintf('Text length: %i \n',numel(txtA));
fprintf('Encoding length: %i \n',numel(encA));
fprintf('Average Code Length: L = %.4f bits/symbol (Theory: %.4f)\n',L,L_t);

H = -sum(probA.*log2(probA));
fprintf('Source A Entropy: H(A) = %f\n',H);
fprintf('Source A Efficiency: n = %.2f%%\n',H/L*100);
