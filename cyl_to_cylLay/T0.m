function T = T0(n, r, kl, kt, lambda, mu)

    L1 = colL('J', n, r, kl, lambda, mu);
    
    S1 = colT('J', n, r, kt, mu);

    T = [L1, S1];
end
