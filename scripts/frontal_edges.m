function frontal_vertices(RGB, S, f_edges, dirichlet, fname_out)

% get number of vertices and images
[nv nb_im] = size(S);

% load edges
edges = load(f_edges);
ne = length(edges);

% graph matrices
A = sparse(edges(:,1)+1, edges(:,2)+1, 1, nv, nv);
A = A + A':
B = incidence(A);
C = abs(B)/2;

% initialise output to 0
out = zeros(1, nv, 3);

% index of the most frontal image for each vertex
[~, idx] = max(S, [], 2);

% loop over the images to fill the output with the colour of the most frontal
% image for each vertex
d = zeros(ne,1); 
for l = 1 : 3
    for i = 1 : nb_im
        Mv = sparse(1:nv, 1:nv, idx==i);
        Me = sparse(1:ne, 1:ne, C*Mv*ones(nv,1));
        im = RGB(:,i,l);
        im(isnan(im)) = 1;
        d = d + Me*(B*im./(C*im));
    end
    D = sparse(1:ne, 1:ne, d);
    linA = - B' * (B - D * C);
    out(:,:,l) = out(:,:,l) + M*RGB(:,i,l);

    cdir = zeros(size(im)); cdir(dirichlet,1) = 1;
    Mdir = sparse(1:nv, 1:nv, cdir);

    im_ref = RGB(:,1,l); im_ref(isnan(im_ref)) = 1;
    linA = (speye(nv) - Mdir) * linA + Mdir;
    b = Mdir * im_ref; 

    out(:,:,i) = linA \ b;
end

imwrite_with_tiff(out, fname_out);

end
















