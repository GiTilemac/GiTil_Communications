%
% Author: Tilemachos S. Doganis
%
function [ c,d, I ] = sortCol( a,b )
	% Sort a and b according to a (descending order)
	% Store transmutations in I
	[a,I] = sort(a,'descend');
	temp = b;
	N = size(a,1);
	for i = 1:N
		b(i) = temp(I(i));
	end
	c = a;
	d = b;
end