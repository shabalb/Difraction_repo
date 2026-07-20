function T = Tmat(n, r, kl, kt, lambda, mu)
%TMAT Matrix relating layer coefficients to displacement-stress components.
% Columns correspond to J_n(kl*r), Y_n(kl*r), J_n(kt*r), Y_n(kt*r).

    L1 = colL('J', n, r, kl, lambda, mu);
    L2 = colL('Y', n, r, kl, lambda, mu);

    S1 = colT('J', n, r, kt, mu);
    S2 = colT('Y', n, r, kt, mu);

    T = [L1, L2, S1, S2];
end
