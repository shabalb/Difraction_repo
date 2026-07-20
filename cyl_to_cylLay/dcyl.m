function y = dcyl(kind, n, z)
% Derivative with respect to the argument z.
    y = 0.5 * (cyl(kind, n-1, z) - cyl(kind, n+1, z));
end