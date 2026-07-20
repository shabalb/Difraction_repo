function y = d2cyl_dr2(kind, n, r, q)
% Second derivative with respect to r of cyl(kind,n,q*r).
    z = q*r;
    Z = cyl(kind, n, z);
    dZdr = q * dcyl(kind, n, z);
    y = (n^2/r^2 - q^2)*Z - dZdr/r;
end