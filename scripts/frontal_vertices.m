function frontal_vertices(RGB, S, fname_out)

% get number of vertices and images
[nv nb_im] = size(S);

% initialise output to 0
out = zeros(1, nv, 3);

% index of the most frontal image for each vertex
[~, idx] = max(S, [], 2);

% loop over the images to fill the output with the colour of the most frontal
% image for each vertex
for i = 1 : nb_im
    M = sparse(1:nv, 1:nv, idx==i);
    for l = 1 : 3 
        out(:,:,l) = out(:,:,l) + (M*RGB(:,i,l))';
    end
end
fname_out
imwrite_with_tiff(out, fname_out);
end
















