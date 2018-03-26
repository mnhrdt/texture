function B = incidence(A)             
	[i,j] = find(triu(A));        % find the (i,j) positions
	n = size(A,1);                  % number of vertices
	m = size(i,1);                  % total number of edges
	B1 = sparse(1:m, i, 1, m, n); % matrix for destination vertices
	B2 = sparse(1:m, j, 1, m, n); % matrix for source vertices
	B = B1 - B2;                  % signed incidence matrix
end