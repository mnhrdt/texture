function thresh = get_shadow_threshold(filename_proj, filename_crop)

    a = imread(filename_proj);
    b = imread(filename_crop);

    t = 150;
    c = (100 * (b <= t) + 255 * (b > t)) .* (a > 0);
    l1_old = sum(sum((a==c).*(a>0)));
    l1 = l1_old;
    T = zeros(1, 850);
    L = zeros(1, 850);
    i = 0;
    
    while (l1_old >= l1 && t < 1000)
        t = t + 1;
        i = i + 1;
        c = (100 * (b <= t) + 255 * (b > t)) .* (a > 0);
        l1_new = sum(sum((a==c).*(a>0)));
        diff = l1_new - l1_old;
        l1_old = l1_new;
        T(i) = t;
        L(i) = l1_new;
    end

%    printf("%f", t);
    [m, idx] = max(L);
    thresh = 150 + idx;
end


    
