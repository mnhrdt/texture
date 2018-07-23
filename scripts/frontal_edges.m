function frontal_edges(RGB, SCALARS, B, C, dirichlet, fname_out)
% takes for each vertex the colour from the most frontal image
% INPUTS: RGB image of size nv x nb_im x 3 
%         



% get number of edges and images
[ne, nb_im] = size(SCALARS);

% get number of vertices
nv = size(B,2);

% initialise output to 0
out = zeros(1, nv, 3);

% index of the most frontal image for each edge
[~, idx] = max(SCALARS, [], 2);

% loop over the rgb channels
for l = 1 : 3
    tic
    % initialise the drift-fields to zeros
    DG = zeros(ne,nb_im);
    
    % loop over the images to compute their drift-field
    for i = 1 : nb_im
        im = RGB(:,i,l);
        % getting rid of the abherent colours
        im(im <= 0) = 0;
        % computing drift-field while keeping NaNs
        DG(:,i) = B * im ./ (C * im);
    end

    % choosing for each edge the drift-field corresponding to 
    % the most frontal image
    k = sub2ind(size(DG), (1:ne)', idx);
    dg = DG(k);
    
    % replacing NaNs by 0 (pure diffusion when no information)
    dg(isnan(dg)) = 0;
    D = sparse(1:ne, 1:ne, dg);

    % dirichlet boundary conditions
    cdir = zeros(size(im)); cdir(dirichlet,1) = 1;
    Mdir = sparse(1:nv, 1:nv, cdir);
    im_ref = RGB(:,1,l); 
    im_ref(isnan(im_ref)) = 1; im_ref(im_ref >= 0) = 1;
    b = Mdir * im_ref; 

    % (laplacian - div(d )) on the roi
    % identity for boundary conditions
    A = (speye(nv) - Mdir) * (- B' * (B - D * C)) + Mdir;

    % solving the osmosis equation
    out(:,:,l) = A \ b;
    toc
end

% write result
imwrite_with_tiff(out, fname_out);

end
















