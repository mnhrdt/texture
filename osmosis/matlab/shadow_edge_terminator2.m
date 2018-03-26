function out = shadow_edge_terminator2(fname_in1, fname_in2, fname_m, fname_e)
% entrées = background image
%           foreground image
%           masque (gris à remplacer, rouge conditions neumann)
%           bords d'ombre où annuler le drift-field
tic
    G = double(1+(imread(fname_in1)));
    %G = G(2800:3400,700:2500,:);
    F = double(1+(imread(fname_in2)));
    %F = F(2800:3400,700:2500,:);
    mask = double(imread(fname_m));
    %mask = mask(2800:3400,700:2500,:);
    e = double(imread(fname_e));
    %e = e(2800:3400,700:2500,:);
    e = double(e(:,:,1)>10);
    m = double(mask(:,:,2) > 10);
    t = double(mask(:,:,1) > 250);
    %m = double(mask(:,:,3) > 10);
    %t = double(mask(:,:,2) > 250);
    me = (m+e)>0;
    ind = find(t(:));
    [w,h,pd] = size(G);
    I = speye(w*h); 
    A = grid_graph(w, h);
    %for k = 1:length(ind)
    %    for l = 1:length(ind)
    %        A(ind(k),ind(l)) = 0;
    %    end
    %end
    A(ind,:) = 0;
    A(:,ind) = 0;
    for k = 1:length(ind)
        A(ind(k),ind(k)) = 1;
    end
    grad = incidence(A);
    div = -grad';
    lap = div * grad;
    C = abs(grad)/2;
    M = spdiags((m(:)+e(:)) > 0, 0, w*h, w*h);
    P = optimizing_permutation(double(A), double(m(:)));
    %P = optimizing_permutation(double(A), double(e(:)));
    [mx, my] = find(m);
    out = zeros(w, h, pd);
    for l = 1:pd
        g = G(:,:,l);
        f = F(:,:,l);
        f = double(f(:));
        g = double(g(:));
        dg = (grad * g) ./ (C * g);
        df = (grad * f) ./ (C * f);
        df = df.*(1-C*e(:));
        %g = f.*me + g.*(1-me);
        %dg = dg.*(1-C*e(:));
        ne = length(df);
        D = sparse(1:ne, 1:ne, df);
        %corner = sparse(w*h,w*h);
        %corner(1,1) = 1;
        %A =  corner + (I-corner)*div*(grad - D*C);
        A =  I - M + M*div*(grad - D*C);
        %b = zeros(w*h,1);
        %b(1) = g(1);
        b = (I - M)*g;
        x = P' * ( (P*A*P') \ (P*b) );
        out(:,:,l) = reshape(x, w, h);
    end
    %out = hsv2rgb(out);
    toc
end
