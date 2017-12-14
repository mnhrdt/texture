function miniodemo(fname_in, fname_e, fname_s, fname_t, x, y, eq_id, fu_id, fname_out)
	%morsi disk2 gradient mask.png edge.tif
    %miniodemo('k.tif', 'edge.tif', 'mask.png', 'terminateur.png', 0, 0, 'shadow', 'shd', 'out.tif')
    epsilon = 0.001;
    g = double(1+imread(fname_in));     % read background image
	f = imread(fname_e);                % read edge mask 
	m = imread(fname_s);                % read shadow mask
	m = double(m>0);
    t = imread(fname_t);
    t = double(t(:,:,1)>0);
    [w,h,pd] = size(g);                 % get reference (background) size
	out = zeros(w, h, pd);              % create output image
	for l = 1:pd                        % call the solver for each band
		out(:,:,l) = gray_cloning(m, t, f(:,:,l), g(:,:,l), eq_id, fu_id);
	end
	imwrite_with_tiff(out-1, fname_out);   % write the output image
	%dlmwrite(fname_out, out-1);   % write the output image
	%imwrite((out-1), fname_out);   % write the output image
end

function out = gray_cloning(m, t, f, g, method_id, fusion_id, epsilon)
	[w,h] = size(m);                     % extract problem dimensions
	f = double(f(:));                    % extract interior data vector
	g = double(g(:));                    % extract boundary data vector
	M = spdiags((double(m(:))+f-1000*t(:)) > 0, 0, w*h, w*h);  % build mask operator
	I = speye(w*h);                      % build identity operator
	A = grid_graph_with_t(w, h, t(:));                % build adjacency matrix
	grad = incidence(A);                 % build gradient operator
	div = -grad';                        % build divergence operator
	lap = div * grad;                    % build laplacian operator
	C = abs(grad)/2;                     % build centering operator
	P = optimizing_permutation(A, m(:)); % build permutation operator
	switch fusion_id                     % select a fusion criterion
		case "sum"  
            fusion = @(x,y) x + y;
		case "avg"  
            fusion = @(x,y) (x + y)/2;
		case "max"  
            fusion = @(x,y) x + (y - x) .* (abs(y) > abs(x));
        case "shd"  
            fusion = @(x,y) x .*(y==0);
        otherwise
            fusion = @(x,y) x;
	end
	switch method_id                     % build the linear system "Ax=b"
        case "shadow"
			dg = (grad * g) ./ (C * g);
            ne = length(dg);
            D = sparse(1:ne, 1:ne, fusion(dg, (C*f)>0));
			A =  I - M + M*div*(grad - D*C);
			b = (I - M)*g;
		case "poisson"
			d = fusion(grad*f, grad*g);
			A =  I - M    - M*lap  ;
			b = (I - M)*g - M*div*d;
		case "logpoisson"
			d = fusion(grad*log(f), grad*log(g));
			A =  I - M    - M*lap  ;
			b = (I - M)*log(g) - M*div*d;
		case "osmosis"
			df = (grad * f) ./ (C * f); % drift vector field
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
			%epsilon = 0.0001;
			X = grad - D*C;
			A =  I - M + M*((1-epsilon)*X'*X - epsilon*lap);
			b = (I - M)*g;
		otherwise
			A = I; % trivial (copy-paste) operator
			b = (I - M)*g + M*fusion(f,g);
	end
	x = P' * ( (P*A*P') \ (P*b) );       % solve the system "x=A\b"
	if strcmp("logpoisson", method_id)
		x = exp(x);
	end
	out = reshape(x, w, h);              % return a rectangular image
end

function A = grid_graph(w, h)                      % build a grid graph WxH
	px = sparse(1:w-1, 2:w, 1, w, w);          % path graph of length W
	py = sparse(1:h-1, 2:h, 1, h, h);          % path graph of length H
	A = kron(py,speye(w)) + kron(speye(h),px); % kronecker sum
	A = A + A';                                % symmetrization
end

function A = grid_graph_with_t(w, h, t)                      % build a grid graph WxH
	px = sparse(1:w-1, 2:w, 1, w, w);          % path graph of length W
	py = sparse(1:h-1, 2:h, 1, h, h);          % path graph of length H
	A = kron(py,speye(w)) + kron(speye(h),px); % kronecker sum
	A = A + A';                                % symmetrization
    A(t>0,:) = 0;
    A(:,t>0) = 0;
end

function B = incidence(A)             % build incidence from adjacency matrix
	[i,j] = find(triu(A));        % find the (i,j) positions
	n = size(A,1);                  % number of vertices
	m = size(i,1);                  % total number of edges
	B1 = sparse(1:m, i, 1, m, n); % matrix for destination vertices
	B2 = sparse(1:m, j, 1, m, n); % matrix for source vertices
	B = B1 - B2;                  % signed incidence matrix
end

function P = optimizing_permutation(A, m) % build a permutating preconditioner
	m = A^3 * m;                      % dilate the mask by 3 pixels
	p = [find(m); find(~m)];          % permutation indices
	P = sparse(1:size(m), p, 1);      % permutation matrix
end
