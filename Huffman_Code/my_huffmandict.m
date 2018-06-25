%% Huffman Dictionary Construction
% Author: Tilemachos S. Doganis

function [ HuffDict ] = my_huffmandict( prob, alph )
%
% Input: Symbol alphabet, Symbol probabilities
% Output: Cell containing alphabet and coding
%
% Other variables:
% P contains the probabilities of symbols and their combinations
% T is a record of transitions and merges during Huffman tree construction
% sortCol(a,b) sorts matrix a in descending order, applies the same
% transitions in b and stores them in I.
%
%% Constructing the Huffman tree
% Initially the symbol probabilities are sorted and stored in the 1st
% column of P, with matrix A recording the link between each symbol and
% its probability. The first column of T is initialized as 1:N.
% 1. The last two probabilities are added and the resulting sum P' is stored
% in the (N-1)th position of the next column. The rest of the column is 
% transferred to the next one as is.
% 2. The whole column is sorted in descending order.
% 3. If the summed probability P' is equal to an existing one, Bubble Sort
% is used to transfer P' to the top of other equal probabilities, in order
% to minimize the variance.
%
alph = cellstr(transpose(alph));

N = numel(prob);
P = zeros(N); % Probability matrix
T = zeros(N); % Transition matrix

% Initialize matrices
P(:,1) = prob;
T(:,1) = (1:N)';
A = (1:N)';
[P(:,1),A] = sortCol(P(:,1),A);

% Starting at the first column, add the last 2 elements
for m = N:-1:2
   % Add the last 2 elements' probabilities
   P(m-1,N-m+2) = P(m,N-m+1)+P(m-1,N-m+1);
   
   % Combine the last 2 elements into the minimum
   T(m-1,N-m+2) = m-1; 
   
   % Transfer the rest of the column as is
   P(1:(m-2),N-m+2) = P(1:(m-2),N-m+1);
   T(1:(m-2),N-m+2) = 1:(m-2);

   % Sort the column using the Probability matrix
   [P(1:(m-1),N-m+2),T(1:(m-1),N-m+2),I] = sortCol(P(1:(m-1),N-m+2),T(1:(m-1),N-m+2));

   % Bubble sort the combined element to the top of others with equal
   % probability
   val = find(T(:,N-m+2) == max(T(:,N-m+2)) );
   if (val>1)
       while (P(val,N-m+2)==P(val-1,N-m+2))
           temp = I(val);
           I(val)=I(val-1);
           I(val-1)=temp;
           temp = P(val);
           P(val)=P(val-1);
           P(val-1)=temp;
           val=val-1;
       end
   end
   T(1:m-1,N-m+2) = I;
   
end
%% Produce Huffman Coding
%
% Matrix Enc will contain the coding for each symbol.
% Variable 'pos' keeps track of the current symbol's transitions between
% columns
% 1. Beginning at the 1st column, each symbol's path during the tree
% construction is backtracked between columns of Enc.
% 2. At every merge point a 0 or a 1 is added to the least significant bit
% of that symbol's coding.
%

% Initialize Coding matrix with 9's
Enc = repmat('9',N,N);

% Calculate Huffman Coding for each symbol (element of first column)
for i=1:N
    % Check each element of first column
    pos = i;
    for j=2:N
        if (pos == N-j+1)
            % If it's second-to-last, add '0' to coding
            Enc(i,N-j+2) = '0';
            
        elseif (pos == N-j+2)
            % If symbol is at bottom of column, it will be combined with
            % the above symbol in the next step, therefore the position
            % index is reduced by 1 before adding '1' to the coding.
            Enc(i,N-j+2) = '1';
            pos = pos-1;
        end
        % Find symbol's position in next column of T
        pos = find( T(:,j) == pos);
    end
end
Enc = cellstr(Enc); % Convert to array of strings
Enc = regexprep(Enc, '9', ''); %Replace 9's with blanks

% Matrix A contains the initial symbol probability sorting information, and
% is used to link each probability to the corresponding symbol
Enc2(A(:)) = Enc;
Enc2 = Enc2';
HuffDict = [ alph  Enc2 ];
