function W = compute_W_j(r, k_l, k_t, lambda, mu, n)

    xl = k_l * r;
    xt = k_t * r;
    Jl = besselj(n, xl);
    Nl = bessely(n, xl);
    Jt = besselj(n, xt);
    Nt = bessely(n, xt);
    
    dJl = k_l * 0.5 * (besselj(n-1, xl) - besselj(n+1, xl));
    dNl = k_l * 0.5 * (bessely(n-1, xl) - bessely(n+1, xl));
    dJt = k_t * 0.5 * (besselj(n-1, xt) - besselj(n+1, xt));
    dNt = k_t * 0.5 * (bessely(n-1, xt) - bessely(n+1, xt));

    d2Jl = ( (n^2/r^2) - k_l^2 ) * Jl;
    d2Nl = ( (n^2/r^2) - k_l^2 ) * Nl;

    d2Jt = ( (n^2/r^2) - k_t^2 ) * Jt;
    d2Nt = ( (n^2/r^2) - k_t^2 ) * Nt;
    
    W = zeros(4,4);

    W(1,1) = dJl;
    W(1,2) = dNl;
    W(1,3) = 1i * n / r * Jt;
    W(1,4) = 1i * n / r * Nt;
    

    W(2,1) = 1i * n / r * Jl;
    W(2,2) = 1i * n / r * Nl;
    W(2,3) = -dJt;
    W(2,4) = -dNt;
    

    W(3,1) = lambda * d2Jl + (1/r^2)*(lambda + 2*mu) * (-n^2 * Jl) + (1/r)*(lambda + 2*mu) * dJl;
    W(3,2) = lambda * d2Nl + (1/r^2)*(lambda + 2*mu) * (-n^2 * Nl) + (1/r)*(lambda + 2*mu) * dNl;

    W(3,3) = (2*mu / r^2) * (1i * n * Jt) - (2*mu / r) * (1i * n * dJt);
    W(3,4) = (2*mu / r^2) * (1i * n * Nt) - (2*mu / r) * (1i * n * dNt);

    W(4,1) = 2*mu * (1i*n/r) * (dJl - Jl/r);
    W(4,2) = 2*mu * (1i*n/r) * (dNl - Nl/r);
    W(4,3) = mu * ( (1/r)*dJt - d2Jt - (n^2/r^2)*Jt );
    W(4,4) = mu * ( (1/r)*dNt - d2Nt - (n^2/r^2)*Nt );
end