function c = colT(kind, n, r, q, mu)

    Z = cyl(kind, n, q*r);
    dZ = q * dcyl(kind, n, q*r);
    ddZ = d2cyl_dr2(kind, n, r, q);

    c = [
        1i*n/r * Z;
        -dZ;
        2i*n*mu/r * (dZ - Z/r);
        mu * (dZ/r - ddZ - n^2/r^2*Z)
    ];
end