function scalar_vertices(RGB, S, fname_out)

% get number of vertices and images
[nv nb_im] = size(S);

% initialise output to 0
out = zeros(1, nv, 3);

% power for weights and initialise dominator
p = 3;

% loop over the images to fill the output with the mean pondered by the scalar
% of the colour for each vertex
% image for each vertex
for l = 1 : 3
    den = zeros(nv,1);
    sumC = zeros(nv,1);
    for i = 1 : nb_im 
        sumC = sumC + (S .^ p) .* RGB(:,i,l);
        den = den + (S .^ p);
    end
    out(1,:,l) = sumC ./ den;
end
imwrite_with_tiff(out, fname_out);
end
















