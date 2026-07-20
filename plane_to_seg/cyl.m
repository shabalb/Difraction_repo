function y = cyl(kind, n, z)
    switch kind
        case 'J'
            y = besselj(n, z);
        case 'Y'
            y = bessely(n, z);
        case 'H'
            y = besselh(n, 1, z);
        otherwise
            error('Unknown cylindrical function kind: %s', kind);
    end
end