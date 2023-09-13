function c = GraphColor(A)
% function c = GraphColor(A)
%
% Greedy graph (vertex) coloring.
% Input: A, a sparse matrix.
% Output: c, an array of integers corresponding to a vertex coloring.
%
% Colors start at 1. (0 means a vertex has not been colored yet.)

% Check for square matrix
[m,n] = size(A);
if (m~=n)
  error('Input matrix must be square!');
end

% Symmetrize
A = abs(A)+abs(A)';

% Greedy coloring. Use natural order.
% TODO: Add other orderings that may reduce #colors.
c = zeros(n,1);
c(1) = 1;
for j = 2:n
  avail = ones(1,max(c)+1);
  colored_nbors = find(A(1:j-1,j));
  avail(c(colored_nbors)) = 0;
  c(j) = find(avail,1);
end

