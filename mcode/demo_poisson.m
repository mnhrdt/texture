function miniodemo(fname_in1, fname_in0, fname_m, x, y, eq_id, fu_id,
			fname_out, epsilon="0.001")
	g = double(1+imread(fname_in1));    # read background image
	f = double(1+imread(fname_in0));    # read foreground image
	m = imread(fname_m);                # read mask
	[mx, my] = find(m);                 # extract barycenter of mask
	x = x - round(mean(mx));            # x-shift displacement by barycenter
	y = y - round(mean(my));            # y-shift displacement by barycenter
	[w,h,pd] = size(g);                 # get reference (background) size
	[fw,fh,fpd] = size(f);              # get foreground size
	tmp_m = zeros(w, h, 1);             # create mask of the reference size
	tmp_m(1:fw,1:fh) = m;               # copy mask values
	m = 0 < circshift(tmp_m, [x, y]);   # shift mask at requested position
	tmp_f = zeros(w, h, pd);            # create foreground of the ref. size
	tmp_f(1:fw,1:fh,:) = f;             # copy foreground values
	f = circshift(tmp_f, [x, y]);       # shift foreground
	out = zeros(w, h, pd);              # create output image
	for l = 1:pd                        # call the solver for each band
		out(:,:,l) = gray_cloning(m, f(:,:,l), g(:,:,l), eq_id, fu_id,
								epsilon);
	end
	imwrite(uint8(out-1), fname_out);   # write the output image
end

function out = gray_cloning(m, f, g, method_id, fusion_id, epsilon)
	[w,h] = size(m);                     # extract problem dimensions
	f = double(f(:));                    # extract interior data vector
	g = double(g(:));                    # extract boundary data vector
	M = spdiags(m(:) > 0, 0, w*h, w*h);  # build mask operator
	I = speye(w*h);                      # build identity operator
	A = grid_graph(w, h);                # build adjacency matrix
	grad = incidence(A);                 # build gradient operator
	div = -grad';                        # build divergence operator
	lap = div * grad;                    # build laplacian operator
	C = abs(grad)/2;                     # build centering operator
	P = optimizing_permutation(A, m(:)); # build permutation operator
	switch fusion_id                     # select a fusion criterion
		case "sum"  fusion = @(x,y) x + y;
		case "avg"  fusion = @(x,y) (x + y)/2;;
		case "max"  fusion = @(x,y) x + (y - x) .* (abs(y) > abs(x));
		otherwise   fusion = @(x,y) x;
	end
	switch method_id                     # build the linear system "Ax=b"
		case "poisson"
			d = fusion(grad*f, grad*g);
			A =  I - M    - M*lap  ;
			b = (I - M)*g - M*div*d;
		case "logpoisson"
			d = fusion(grad*log(f), grad*log(g));
			A =  I - M    - M*lap  ;
			b = (I - M)*log(g) - M*div*d;
		case "osmosis"
			df = (grad * f) ./ (C * f); # drift vector field
			dg = (grad * g) ./ (C * g);
			D = diag(fusion(df, dg));
			A =  I - M + M*div*(grad - D*C);
			b = (I - M)*g;
		case "osmovar"
			df = (grad * f) ./ (C * f);
			dg = (grad * g) ./ (C * g);
			D = diag(fusion(df, dg));
			A =  I - M + M*(grad - D*C)'*(grad - D*C);
			b = (I - M)*g;
		case "osmovareps"
			df = (grad * f) ./ (C * f);
			dg = (grad * g) ./ (C * g);
			d = fusion(df, dg);
			D = diag(d);
			#epsilon = 0.0001;
			X = grad - D*C;
			A =  I - M + M*((1-epsilon)*X'*X - epsilon*lap);
			b = (I - M)*g;
		otherwise
			A = I; # trivial (copy-paste) operator
			b = (I - M)*g + M*fusion(f,g);
	end
	x = P' * ( (P*A*P') \ (P*b) );       # solve the system "x=A\b"
	if strcmp("logpoisson", method_id)
		x = exp(x);
	end
	out = reshape(x, w, h);              # return a rectangular image
end

function A = grid_graph(w, h)                      # build a grid graph WxH
	px = sparse(1:w-1, 2:w, 1, w, w);          # path graph of length W
	py = sparse(1:h-1, 2:h, 1, h, h);          # path graph of length H
	A = kron(py,speye(w)) + kron(speye(h),px); # kronecker sum
	A = A + A';                                # symmetrization
end

function B = incidence(A)             # build incidence from adjacency matrix
	[i,j] = find(triu(A));        # find the (i,j) positions
	n = rows(A);                  # number of vertices
	m = rows(i);                  # total number of edges
	B1 = sparse(1:m, i, 1, m, n); # matrix for destination vertices
	B2 = sparse(1:m, j, 1, m, n); # matrix for source vertices
	B = B1 - B2;                  # signed incidence matrix
end

function P = optimizing_permutation(A, m) # build a permutating preconditioner
	m = A^3 * m;                      # dilate the mask by 3 pixels
	p = [find(m); find(!m)];          # permutation indices
	P = sparse(1:size(m), p, 1);      # permutation matrix
end
