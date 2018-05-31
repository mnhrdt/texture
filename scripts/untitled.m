%% essai de récupération d'image par osmose en gardant un pixel.

image = imread('../output/rgb_02.tif');
out = zeros(size(image));
for i = 1:3
im = image(:,:,i)';
%im(isnan(im)) = 1;
edges = load('../edges.txt');

nv = length(im);
ne = length(edges);
A = sparse(edges(:,1)+1, edges(:,2)+1, 1, nv, nv);
A = A + A';
ind = find(isnan(image(:,:,i)));
A(:, ind) = 0;
A(ind, :) = 0;
A = A + sparse
B = incidence(A);
C = abs(B)/2;
div = -B';
d = B*im./(C*im);
D = sparse(1:ne, 1:ne, d);

linA = div*(B-D*C);
linA(1,:) = 0; linA(1,1) = 1;
b = zeros(size(im));
b(1,1) = im(1,1);

out(:,:,i) = linA \ b;
end
imwrite_with_tiff(out, '../output/simple_osmosis.tif')

%% osmose avec deux images

image1 = double(1+imread('../output/rgb_02.tif'));
image2 = double(1+imread('../output/rgb_08.tif'));
out = zeros(size(image1));
for i = 1:3
im1 = image1(:,:,i)';
im1(isnan(im1)) = 1;
im2 = image2(:,:,i)';
im2(isnan(im2)) = 1;
edges = load('../edges.txt');

nv = length(im1);
%ne = length(edges);
A = sparse(edges(:,1)+1, edges(:,2)+1, 1, nv, nv);
A = A + A';
liste = (isnan(image1(:,:,i)')+...
    isnan(image2(:,:,i)') == 2);
ind = find(liste);
A(ind,:) = 0;
A(:,ind) = 0;
%A = A + sparse(1:nv, 1:nv, liste);
B = incidence(A);
C = abs(B)/2;
div = -B';
ne = size(B,1);
d1 = B*im1./(C*im1);
d2 = B*im2./(C*im2);
d = max(d1.*(~isnan(C*image1(:,:,i)')),d2.*(~isnan(C*image2(:,:,i)')));
%  d = (d1.*(~isnan(C*image1(:,:,i)'))+...
%      d2.*(~isnan(C*image2(:,:,i)')))./...
%      max(1,(~isnan(C*image1(:,:,i)'))+...
%     (~isnan(C*image2(:,:,i)')));
D = sparse(1:ne, 1:ne, d);
tic
linA = div*(B-D*C);
linA(1,:) = 0; linA(1,1) = 1;
b = zeros(size(im1));
b(1,1) = im1(1,1);

out(:,:,i) = linA \ b;
toc
end
imwrite_with_tiff(out, '../output/complex_osmosis.tif')
























