%% Huffman Coding Test 3 (Biased Source, Source-specific Probabilities)
%
% Author: Tilemachos S. Doganis
%
% This test uses a set of words beginning with the letter 'k' as a source,
% and the a priori letter probabilities for the English alphabet.
%
% As the theoretical a priori probabilities for the whole English language
% were proven to be ineffective in the 2nd Test, this time the
% probabilities are calculated using a frequency histogram of the input.
% Dividing each frequency with the number of characters results in the
% probabilities.
%
% The results show much improvement compared to the a priori probabilities,
% such as the reduction of the coding length by 13974 bits, and of the
% average coding length by approximately 0.5 bits/symbol. Furthermore, due
% to the biased form of the source input, the source entropy (uncertainty)
% is reduced.
%
clc
clearvars -except txtA
fprintf('Running Huffman Encoding of Source B with Source-Specific Probabilities... \n');

% Load kwords.txt in matrix T
fid = fopen('kwords.txt','r');
T = textscan(fid,'%s');
fclose(fid);
T=T{1};

% Store ASCII table in matrix asc
freq = cell(1,2);
freq{1} = transpose('!':'~');
freq{2} = zeros(94,1);

% For every word in T
total_chars = 0;
tsize = numel(T);
for i=1:tsize
    cur_word = T{i};
    m = numel(cur_word);
    total_chars = total_chars + m;
   
    %Count each character
    for j=1:m
        pos = strfind(transpose(freq{1}),cur_word(j));
        freq{2}(pos) = freq{2}(pos) + 1;
    end
    
end

% Estimate probability
freq{2} = freq{2}/total_chars;

% Remove matrix rows where characters have 0 probability
n = numel(freq{1});
i=1;
while (i<=n)
   if (freq{2}(i)==0)
      freq{1}(i) = [];
      freq{2}(i) = [];
      n = n-1;
   else
      i = i+1;
   end
end

% Alphabet and Probabilities for Source B
alph = freq{1}';
prob = freq{2}';

%Create Huffman Dictionary
dict = my_huffmandict(prob,alph);

% Calculate Theoretical Average Code Length
c = size(dict,1);
L_t = 0;
for i = 1:c
    L_t = L_t + (numel(dict{i,2})*prob(i));
end

% Store list T as string txtB
txtB = T{1};
tsize = numel(T);
for i=2:tsize
   txtB = strcat(txtB,T(i)); 
end
txtB = char(txtB);
total_chars = numel(txtB);

% Encode Text using Huffman Coding
enc = my_huffmanenc(txtB,dict);
enc = enc{1};
encLenB = numel(enc);

% Decode Code
dec = my_huffmandec(enc,dict);
dec = dec{1};

%Check if encoding / decoding are successful
if (strcmp(txtB,dec) == 1)
    disp('Encoding & Decoding Successful');
else
    disp('Encoding / Decoding Error');
end

L = encLenB/total_chars;
fprintf('Text length: %i \n',total_chars);
fprintf('Encoding length: %i \n',encLenB);
fprintf('Average Code Length: L = %.4f bits/symbol (Theory: %.4f)\n',L,L_t);

H = -sum(prob.*log2(prob));
fprintf('Source B Entropy: H(B) = %f\n',H);
fprintf('Source B Efficiency: n = %.2f%%\n',H/L*100);
