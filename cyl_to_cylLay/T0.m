function T = T0(n, r, kl, kt, lambda, mu)
%T0 Regular core matrix with columns J_n(kl*r) and J_n(kt*r).

    L1 = colL('J', n, r, kl, lambda, mu);
    
    S1 = colT('J', n, r, kt, mu);

    T = [L1, S1];
end
