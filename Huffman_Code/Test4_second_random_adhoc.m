%% Huffman Coding Test 4 (Random 2nd-Extension Source, A Priori Symbol Probabilities)
% Author: Tilemachos S. Doganis
%
% This test is similar to Test 1, but the source is treated as a
% 2nd-extensions source, which means that symbols are coded in pairs
% instead of single characters. Pair probabilities are estimated as the
% product of single character probabilities, although this assumes they are
% independant of each other, which is not true in any natural language.
%
% The average coding length as well as the source entropy is cut by
% approximately half, leading to an increase of about 0.28% in efficiency.
%
clc
clearvars -except txtA
fprintf('Running Huffman Encoding of 2nd-Extension Source A with Given Probabilities... \n');

% Probabilities for lower case english letters
% (https://en.wikipedia.org/wiki/Letter_frequency)
probA = [ 8.167 1.492 2.782 4.253 12.702 2.228 2.015 6.094 6.966 0.153 ...
          0.772 4.025 2.406 6.749 7.507 1.929 0.095 5.987 6.327 9.056 ...
          2.758 0.978 2.361 0.150 1.974 0.074 ];
% Normalization
probA = probA/100;

% Create alphabet of lower case English letter pairs
asize = 26;
t = cell(asize^2,1);
for i=1:asize
    t(1+asize*(i-1):asize*i) = cellstr(strcat(char(96+i),('a':'z')'));
end

% Calculate Probabilities for pair values
prob = zeros(asize);
for i=1:asize
   for j=1:asize
       prob(i,j) = probA(i)*probA(j);
   end
end

% Convert Probability Matrix to Single Column
ext_prob = reshape(prob,[asize^2 1]);

% Create Dictionary
ext_prob = ext_prob';
alph = char(t); %Cell Array -> Column of strings
alph = cellstr(alph)';
dict = my_huffmandict( ext_prob, alph);

% Calculate Theoretical Average Code Length
c = size(dict,1);
L_t = 0;
for i = 1:c
    L_t = L_t + (numel(dict{i,2})*ext_prob(i));
end

if exist('txtA','var')
    % Use text from Part 2
    ext_txt = txtA;
else
    % Produce Text
    inp = [double('a':'z') ; probA];
    ext_txt = char(randsrc(1,10000,inp));
end

% Encode Text
enc = my_huffmanenc( ext_txt, dict, 2);
enc = enc{1};

% Decode Text
dec = my_huffmandec ( enc, dict);
dec = dec{1};

% Check if encoding / decoding are successful
if (strcmp(ext_txt,char(dec)) == 1)
    disp('Encoding & Decoding Successful');
else
    disp('Encoding / Decoding Error');
end
L = numel(enc)/5000;
fprintf('Text length: %i \n',10000);
fprintf('Encoding length: %i \n',numel(enc));
fprintf('Average Code Length: L = %.4f bits/symbol (Theory: %.4f)\n',L,L_t);

H = -sum(ext_prob.*log2(ext_prob));
fprintf('Source A Entropy: H(A) = %f (Theory: 8.1892)\n',H);
fprintf('Source A Efficiency: n = %.2f%%\n',H/L*100);
