function P = optimizing_permutation(A, m) 
	m = A^3 * m;                      % dilate the mask by 3 pixels
	p = [find(m); find(~m)];          % permutation indices
	P = sparse(1:size(m), p, 1);      % permutation matrix
end