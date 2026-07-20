function T = Tmat(n, r, kl, kt, lambda, mu)
    L1 = colL('J', n, r, kl, lambda, mu);
    L2 = colL('Y', n, r, kl, lambda, mu);

    S1 = colT('J', n, r, kt, mu);
    S2 = colT('Y', n, r, kt, mu);

    T = [L1, L2, S1, S2];
end
