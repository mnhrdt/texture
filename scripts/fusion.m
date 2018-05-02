a = dir('../output/rgb_*.tif');
a = a([2:19 21:32 35:46]);
b = dir('../output/scalars_*.tif');
b = b([2:19 21:32 35:46]);
c = dir('../output/real_sun*.tif');
c = c([2:19 21:32 35:46]);

edges = load('../edges.txt');

sz = length(a);
M = double(imread(strcat(a(1).folder, '/', a(1).name)));
nv = size(M,2);
ne = size(edges,1);
all_rgb = zeros(sz, nv, 3);
all_scalars = zeros(sz, nv, 3);
all_sun = zeros(sz, nv, 3);
all_drift = zeros(sz, ne,3);
all_scalars_edges = zeros(sz, ne,3);
all_sun_edges = zeros(sz, ne,3);

for i = 1:sz
    strcat(a(i).folder, '/', a(i).name)
    all_rgb(i, :, :) = ...
        double(imread(strcat(a(i).folder, '/', a(i).name)));
    all_scalars(i, :, :) = ...
        double(imread(strcat(b(i).folder, '/', b(i).name)));
    all_sun(i, :, :) = ...
        double(imread(strcat(c(i).folder, '/', c(i).name)));
    [drift, scalars, sun] = ...
        compute_drift_field(all_rgb(i, :, :),all_scalars(i, :, :),...
        all_sun(i, :, :),edges);
    all_drift(i, :, :) = drift;
    all_scalars_edges(i, :, :) = scalars;
    all_sun_edges(i, :, :) = sun;
end

imwrite_with_tiff(all_rgb, '../output/all_rgb.tif');
imwrite_with_tiff(all_scalars, '../output/all_scalars.tif');
imwrite_with_tiff(all_sun, '../output/all_sun.tif');
imwrite_with_tiff(all_drift, '../output/all_drift.tif');
imwrite_with_tiff...
    (all_scalars_edges, '../output/all_scalars_edges.tif');
imwrite_with_tiff(all_sun_edges, '../output/all_sun_edges.tif');


%% weighted mean and median for values
for p = 0:4:30
    p
    mean = (all_scalars.^p).*all_rgb;
    denominateur = sum((1-isnan(all_rgb)).*all_scalars.^p);
    result = nansum(mean)./denominateur;
    str = strcat(a(1).folder, '/mean_', num2str(p), '.tif');
    imwrite_with_tiff(result, str)

    m = repmat(result,42,1,1);

    for i = 1:5
        mean = (all_scalars.^p).*all_rgb./abs(all_rgb-m);
        denominateur = all_scalars.^p./abs(all_rgb-m);
        result = nansum(mean)./nansum(denominateur);
        m = repmat(result,42,1,1);
    end
    str = strcat(a(1).folder, '/median_', num2str(p), '.tif');
    imwrite_with_tiff(result, str)
end

%% weighted mean and median with sun for values
for p = 24:4:28
    p
    mean = (all_scalars.^p).*all_rgb.*all_sun;
    denominateur = sum((1-isnan(all_rgb)).*all_scalars.^p.*all_sun);
    result = nansum(mean)./denominateur;
    str = strcat(a(1).folder, '/mean_sun_', num2str(p), '.tif');
    imwrite_with_tiff(result, str)

    m = repmat(result,42,1,1);

    for i = 1:5
        mean = (all_scalars.^p).*all_rgb./abs(all_rgb-m);
        denominateur = all_scalars.^p./abs(all_rgb-m);
        result = nansum(mean)./nansum(denominateur);
        m = repmat(result,42,1,1);
    end
    str = strcat(a(1).folder, '/median_sun_', num2str(p), '.tif');
    imwrite_with_tiff(result, str)

end










%% weighted mean for drift
for pw = 0:4:30
    pw
    mean = (all_scalars_edges.^pw).*all_drift;
    denominateur = sum((1-isnan(all_drift)).*all_scalars_edges.^pw);
    result = nansum(mean)./denominateur;

    str = strcat(a(1).folder, '/mean_', num2str(pw), '.tif');
    imwrite_with_tiff(result, str)
    
    drift = result;
    m = isnan(all_rgb(1,:,1));
    m(1:end-100) = 1;
    M = sparse(1:nv, 1:nv, m);
    drift(isnan(result))=0;
    ne = size(edges,1);
    nv = size(all_rgb,2);
    I = speye(nv);
    A = sparse(edges(:,1)+1, edges(:,2)+1, 1, nv, nv);
    B1 = sparse(1:ne, edges(:,1)+1, 1, ne, nv);
    B2 = sparse(1:ne, edges(:,2)+1, 1, ne, nv);
    grad = B1 - B2;  
    A = A + A';
    m2 = A^3 * m'; 
    p = [find(m2); find(~m2)];
    P = sparse(1:size(m2), p, 1);
    %grad = incidence(A);
    div = -grad';
    C = abs(grad)/2;
    out = zeros(nv,3);
    for i = 1:3
        tic
        D = sparse(1:ne, 1:ne, drift(:,:,i));
        linA =  I - M + M*div*(grad - D*C);
        b = (I - M)*all_rgb(1,:,i)';
        %vertex = sparse(nv,nv);
        %vertex(1,1) = 1;
        %linA = vertex + (I-vertex)*div*(grad-D*C);
        %b = zeros(nv,1);
        %b(1) = all_rgb(1,1,i);
        x = linA \ b;
        out(:,i) = x;
        toc
    end
    str = strcat(a(1).folder, '/osmosis_', num2str(pw), '.tif')
    imwrite_with_tiff(out, str)
end

%mean2 = (all_scalars.^p2).*all_rgb;
%denominateur2 = sum((1-isnan(all_rgb)).*all_scalars.^p2);
%result2 = nansum(mean2)./denominateur2;


%% weighted mean with sun for drift
p = 10;
mean = (all_scalars_edges.^p).*all_drift.*all_sun_edges;
denominateur = sum((1-isnan(all_drift))...
    .*all_scalars_edges.^p.*all_sun_edges);
result = nansum(mean)./denominateur;

imwrite_with_tiff(result, '../output/drift_mean_10_sun.tif')


%% weighted median with sun for drift
p = 10;
m = zeros(size(all_drift));

for i = 1:5
    mean = (all_scalars_edges.^p)...
        .*all_drift.*all_sun_edges./(all_drift-m);
    denominateur = all_scalars_edges.^p...
        .*all_sun_edges./(all_drift-m);
    m = nansum(mean)./nansum(denominateur);
end

imwrite_with_tiff(m, '../output/drift_median_10_sun.tif')
%%
drift = result;
drift(isnan(result))=0;
ne = size(edges,1);
nv = size(all_rgb,2);
I = speye(nv);
A = sparse(edges(:,1)+1, edges(:,2)+1, 1, nv, nv);
A = A + A';
grad = incidence(A);
div = -grad';
C = abs(grad)/2;
out = zeros(nv,3);
for i = 1:3
    tic
    D = sparse(1:ne, 1:ne, drift(:,:,i));
    vertex = sparse(nv,nv);
    vertex(1,1) = 1;
    linA = vertex + (I-vertex)*div*(grad-D*C);
    b = zeros(nv,1);
    b(1) = all_rgb(1,1,i);
    x = linA \ b;
    out(:,i) = x;
    toc
end
 


