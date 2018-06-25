%% Huffman Coding Test 6 (Biased 2nd-Extension Source, Source-specific Probabilities)
% Author: Tilemachos S. Doganis
%
% This test is similar to Test 6, with Source Specific probabilities
% instead of the a priori ones. This version has increased the efficiency
% for the biased source by 9.96% compared to the a priori dictionary, and
% by 0.58% compared to the Single Source version.
clc
clear
fprintf('Running Huffman Encoding of 2nd-Extension Source B with Source-Specific Probabilities... \n');

% Load kwords.txt in string ext_txt_2
fid = fopen('kwords.txt','r');
T = textscan(fid,'%s');
fclose(fid);
T=T{1};

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

% Ensure that length is even
if ( mod( numel(txtB), 2 ) > 0)
   txtB = strcat(txtB, 'k'); 
end
      
% Create alphabet of lower case English letter pairs
asize = 26;
t = cell(asize^2,1);
for i=1:asize
    t(1+asize*(i-1):asize*i) = cellstr(strcat(char(96+i),transpose('a':'z')));
end

% Initialize Pair Frequency Cell
freq = cell(1,2);
freq{1} = t;
freq{2} = zeros(numel(t),1);

% Count words using dictionary
total_chars = numel(txtB);
for j=1:2:total_chars-1
   k = strcat(txtB(j),txtB(j+1));
   pos1 = strfind('a':'z',k(1));
   pos2 = strfind('a':'z',k(2));
   fin_pos = (pos1-1)*26 + pos2;
   freq{2}(fin_pos) = freq{2}(fin_pos)+1;
end

% Calculate probability
freq{2} = freq{2}/(total_chars/2);

% Remove matrix rows where characters have 0 probability
n = asize^2;
i = 1;
while (i<=n)
   if(freq{2}(i)==0)
      freq{1}(i)=[];
      freq{2}(i)=[];
      i = i-1;
      n = n-1;
   else
     i = i+1;
   end
end

% Alphabet and Probabilities for Source B
alphB = transpose(freq{1});
probB = transpose(freq{2});

% Create Huffman Dictionary
dictB = my_huffmandict(probB,alphB);

% Calculate Theoretical Average Code Length
c = size(dictB,1);
L_t = 0;
for i = 1:c
    L_t = L_t + (numel(dictB{i,2})*probB(i));
end

% Encode Text
encB = my_huffmanenc( txtB, dictB, 2);
encB = cell2mat(encB);
encLenB = numel(encB);

% Decode Text
decB = my_huffmandec ( encB, dictB );
decB = cell2mat(decB);

% Check if same
if (strcmp(txtB,char(decB)) == 1)
    disp('Encoding & Decoding Successful');
else
    disp('Encoding / Decoding Error');
end

L = encLenB/(total_chars/2);
fprintf('Text length: %i \n',total_chars);
fprintf('Encoding length: %i \n',encLenB);
fprintf('Average Code Length: L = %.4f bits/symbol (Theory: %.4f)\n',L,L_t);

H = -sum(probB.*log2(probB));
fprintf('Source B Entropy: H(B) = %f\n',H);
fprintf('Source B Efficiency: n = %.2f%%\n',H/L*100);
