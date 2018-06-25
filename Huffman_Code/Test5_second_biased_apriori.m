%% Huffman Coding Test 5 (Biased 2nd-Extension Source, A Priori Symbol Probabilities)
% Author: Tilemachos S. Doganis
%
% The biased source is now treated as a 2nd-Extension Source, leading to an
% increase of about 0.70% in efficiency, compared to the Single Extension
% version. The inaccurate probabilities still prevent the construction of
% an optimized dictionary though.
%
clc
clear
fprintf('Running Huffman Encoding of 2nd-Extension Source B with Given Probabilities... \n');

% Load kwords.txt in string ext_txt_2
fid = fopen('kwords.txt','r');
T = textscan(fid,'%s');
fclose(fid);
T=T{1};

% Store list T as string txtA
txtA = T{1};
tsize = numel(T);
for i=2:tsize
   txtA = strcat(txtA,T(i)); 
end
txtA = char(txtA);

% Remove special characters
for i = 33:96
    remPos = strfind(txtA,char(i));
    remLen = numel(remPos);
    for j = remLen:-1:1
        txtA(remPos(j))='';
    end
end

%Ensure that length is even
if ( mod( numel(txtA), 2 ) > 0)
   txtA = strcat(txtA, 'k'); 
end
total_chars = numel(txtA);
      
% Probabilities for lower case english letters
probA = [ 8.167 1.492 2.782 4.253 12.702 2.228 2.015 6.094 6.966 0.153 ...
          0.772 4.025 2.406 6.749 7.507 1.929 0.095 5.987 6.327 9.056 ...
          2.758 0.978 2.361 0.150 1.974 0.074 ];
% Normalization
probA = probA/100;

% Create alphabet of lower case English letter pairs
asize = 26;
t = cell(asize^2,1);
for i=1:asize
    t(1+asize*(i-1):asize*i) = cellstr(strcat(char(96+i),transpose('a':'z')));
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

%Create Dictionary
ext_prob = transpose(ext_prob);
alphA = char(t); %Cell Array -> Column of strings
alphA = transpose(cellstr(alphA));
dictA = my_huffmandict( ext_prob, alphA);

% Calculate Theoretical Average Code Length
c = size(dictA,1);
Lt = 0;
for i = 1:c
    Lt = Lt + (numel(dictA{i,2})*ext_prob(i));
end

% Encode Text
encA = my_huffmanenc( txtA, dictA, 2);
encA = encA{1};
encLenA = numel(encA);

% Decode Text
decA = my_huffmandec ( encA, dictA );
decA = decA{1};

% Check if same
if (strcmp(txtA,char(decA)) == 1)
    disp('Encoding & Decoding Successful');
else
    disp('Encoding / Decoding Error');
end

L = encLenA/(total_chars/2);
fprintf('Text length: %i \n',total_chars);
fprintf('Encoding length: %i \n',encLenA);
fprintf('Average Code Length: L = %.4f bits/symbol (Theory: %.4f)\n',L,Lt);

H = -sum(ext_prob.*log2(ext_prob));
fprintf('Source B Entropy: H(B) = %f\n',H);
fprintf('Source B Efficiency: n = %.2f%%\n',H/L*100);
