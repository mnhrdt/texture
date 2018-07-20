function matlab_wrapper(d, contrast, fusion, f_edges, fname_out)

% create structures containing all required images
rgb = dir(strcat(d, '/rgb/rgb_*.tif'));
scalars = dir(strcat(d, '/scalars/scalars_*.tif'));

% get number of images and number of vertices
nb_im = length(rgb);
nv = length(imread(strcat(rgb(1).folder, '/', rgb(1).name)));

% super matrices for all colour images and scalar product for visible vertices or edges
RGB = zeros(nv, nb_im, 3);

if (strcmp(fusion, "frontal_vertices") || strcmp(fusion, "scalar_vertices"))
    SCALARSv = zeros(nv, nb_im);
end

if (strcmp(fusion, "frontal_edges") || strcmp(fusion, "scalar_edges"))
    edges = load(f_edges);
    ne = length(edges);
    SCALARSf = zeros(ne, nb_im);
    A = sparse(edges(:,1)+1, edges(:,2)+1, 1, nv, nv);
    A = A + A';
    B = incidence(A);
    C = abs(B)/2;
end


% loop over the images to fill the above matrices
for i = 1 : nb_im
    % read the images
    im = double(1 + iio_read(strcat(rgb(i).folder,'/', rgb(i).name)));
    scalar = double(iio_read(strcat(scalars(i).folder,'/', scalars(i).name)));

    % fill the colour matrix
    RGB(:,i,:) = permute(im, [2 1 3]);

    % keep scalar product for visible vertices and put the other vertices 
    % scalar product to -10 
    if (strcmp(fusion, "frontal_vertices") || strcmp(fusion, "scalar_vertices"))
        Mv = sparse(1:nv, 1:nv, ~isnan(im(:,:,1)));
        SCALARSv(:, i) = Mv*scalar(:,:,1)' - ...
            10 * (speye(nv) - Mv) * ones(size(scalar(:,:,1)'));
    end
    if (strcmp(fusion, "frontal_edges") || strcmp(fusion, "scalar_edges"))
        s = C * scalar(:,:,1)';
        me = isnan(C * RGB(:,i,1));
        SCALARSf(:, i) = s .* (me == 0) - 10 * me;
    end
end

all_nonan = ~isnan(sum(RGB,2));

switch (contrast)
    case "mean_var"
        nb_nonan = sum(all_nonan);

        im_ref = RGB(:,1,:);
        im_ref(isnan(im_ref)) = 0;
        im_ref = im_ref.*all_nonan;

        mean_ref = sum(im_ref)./nb_nonan;  
        std_ref = sqrt(sum(im_ref.^2)./nb_nonan - mean_ref.^2); 

        rgb_nonan = RGB;
        rgb_nonan(isnan(rgb_nonan)) = 0;
        rgb_nonan = rgb_nonan.*all_nonan;

        mean_rgb = sum(rgb_nonan)./nb_nonan;
        std_rgb = sqrt(sum(rgb_nonan.^2)./nb_nonan - mean_rgb.^2);

        RGB_contrast = ((RGB-mean_rgb)./std_rgb).*std_ref + mean_ref;
end

switch (fusion)
    case "frontal_vertices"
        frontal_vertices(RGB_contrast, SCALARSv, fname_out);

    case "scalar_vertices"
        scalar_vertices(RGB_contrast, SCALARSv, fname_out);

    case "frontal_edges"
        im_ref = RGB_contrast(:,1,:); 
        im_ref(isnan(im_ref)) = 0;
        im_ref = im_ref.*all_nonan;
        F = find(im_ref(:,:,1));
        dirichlet = F(floor(end/2):floor(end/2)+10);
        frontal_edges(RGB_contrast, SCALARSf, B, C, dirichlet, fname_out);

    case "scalar_edges"
        im_ref = RGB_contrast(:,1,:); 
        im_ref(isnan(im_ref)) = 0;
        im_ref = im_ref.*all_nonan;
        F = find(im_ref(:,:,1));
        dirichlet = F(floor(end/2):floor(end/2)+10);
        scalar_edges(RGB_contrast, SCALARSf, B, C, dirichlet, fname_out);
end


end


