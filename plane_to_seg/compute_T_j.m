
function T = compute_T_j(r, k_l, k_t, lambda, mu, n)

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

    T = zeros(4,4);
    T(1,1) = dJl;          
    T(1,2) = dNl;         
    T(1,3) = 1i * n / r * Jt;
    T(1,4) = 1i * n / r * Nt;
    
    T(2,1) = 1i * n / r * Jl;
    T(2,2) = 1i * n / r * Nl;
    T(2,3) = -dJt;
    T(2,4) = -dNt;

    T(3,1) = lambda * (d2Jl + (1/r)*dJl - (n^2/r^2)*Jl) + 2*mu * d2Jl;
    T(3,2) = lambda * (d2Nl + (1/r)*dNl - (n^2/r^2)*Nl) + 2*mu * d2Nl;

    T(3,3) = 2*mu * (1i*n/r) * (dJt - Jt/r);
    T(3,4) = 2*mu * (1i*n/r) * (dNt - Nt/r);

    T(4,1) = 2*mu * (1i*n/r) * (dJl - Jl/r);
    T(4,2) = 2*mu * (1i*n/r) * (dNl - Nl/r);

    T(4,3) = mu * ( (1/r)*dJt - d2Jt - (n^2/r^2)*Jt );
    T(4,4) = mu * ( (1/r)*dNt - d2Nt - (n^2/r^2)*Nt );
end