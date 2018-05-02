function [I, Sc, Sn] = compute_drift_field(img, sc, sn, edges)

    [~, nv, pd] = size(img);
    image = zeros(nv,pd);
    scalars = zeros(nv,pd);
    sun = zeros(nv,pd);
    image(:,:) = img(1,:,:);
    image(isnan(image))=0;
    scalars(:,:) = sc(1,:,:);
    sun(:,:) = sn(1,:,:);
    
    size(image)
    [ne, ~] = size(edges);

    A = sparse(edges(:,1)+1, edges(:,2)+1, 1, nv, nv);
    A = A + A';

    B = incidence(A);
    C = abs(B)/2;

    drift_field = zeros(ne, 3);
    scalars_edges = zeros(ne, 3);
    sun_edges = zeros(ne, 3);

    for i=1:3
        drift_field(:,i) = (B*image(:,i))./(C*image(:,i));
        scalars_edges(:,i) = C*scalars(:,i);
        sun_edges(:,i) = C*sun(:,i);
    end
    
    I = drift_field;
    Sc = scalars_edges;
    Sn = sun_edges;
end

