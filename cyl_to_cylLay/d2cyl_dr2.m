function y = d2cyl_dr2(kind, n, r, q)

    z = q*r;
    Z = cyl(kind, n, z);
    dZdr = q * dcyl(kind, n, z);
    y = (n^2/r^2 - q^2)*Z - dZdr/r;
end