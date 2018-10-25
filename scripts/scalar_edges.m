function scalar_edges(RGB_contrast, SCALARS, B, C, dirichlet, fname_out)

% get number of vertices and images
[ne, nb_im] = size(SCALARS);

% load edges
nv = size(B,2);

% initialise output to 0
out = zeros(1, nv, 3);

% power for the weight
p = 200;

% loop over the images to fill the output with the mean pondered by the scalar
% of the drifts fields for each edge
for l = 1 : 3
    tic
    % initialise nominator and denominator
    den = zeros(ne,1);
    sumD = zeros(ne,1);

    % loop over the images
    for i = 1 : nb_im
        im = RGB_contrast(:,i,l);
        im(im <= 0) = 0;
        sumD = sumD + (SCALARS(:,i) .^ p) .* (B * im ./ (C * im));
        den = den + SCALARS(:,i) .^ p;
        %dd = (B * im ./ (C * im));
        %dd(isnan(dd)) = 0 ;
        %sumD = sumD + (SCALARS(:,i) .^ p) .* dd;
        %den = den + (dd~=0).*SCALARS(:,i) .^ p;
    end

    d = sumD ./ den;
    d(isnan(d)) = 0;

    D = sparse(1:ne, 1:ne, d);

    A = - B' * (B - D * C);

    cdir = zeros(size(im)); cdir(dirichlet,1) = 1;
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
















