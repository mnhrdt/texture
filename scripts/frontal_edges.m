function frontal_edges(RGB_contrast, SCALARS, B, C, dirichlet, fname_out)

% get number of vertices and images
[ne, nb_im] = size(SCALARS);

% load edges
nv = size(B,2);

% initialise output to 0
out = zeros(1, nv, 3);

% index of the most frontal image for each vertex
[~, idx] = max(SCALARS, [], 2);

% loop over the images to fill the output with the colour of the most frontal
% image for each vertex
for l = 1 : 3
    tic
    DG = zeros(ne,nb_im);
    for i = 1 : nb_im
        im = RGB_contrast(:,i,l);
        im(im <= 0) = 0;
        DG(:,i) = B * im ./ (C * im);
    end

    k = sub2ind(size(DG), (1:ne)', idx);
    dg = DG(k);
    dg(isnan(dg)) = 0;

    D = sparse(1:ne, 1:ne, dg);

    A = - B' * (B - D * C);

    cdir = zeros(size(im)); cdir(dirichlet,1) = 1;
    size(cdir)
    nv
    size(dirichlet)
    size(sparse(1:nv, 1:nv,ones(nv,1)))
    Mdir = sparse(1:nv, 1:nv, cdir);
    im_ref = RGB_contrast(:,1,l); 
    im_ref(isnan(im_ref)) = 1; im_ref(im_ref >= 0) = 1;

    A = (speye(nv) - Mdir) * (- B' * (B - D * C)) + Mdir;
    b = Mdir * im_ref; 

    out(:,:,l) = A \ b;
    toc
end

imwrite_with_tiff(out, fname_out);

end
















