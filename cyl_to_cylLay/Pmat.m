function P = Pmat(n, r, kl, kt, lambda, mu)
%PMAT Transfer matrix through all coating layers for one harmonic n.
% r contains layer boundaries: r(1)=r0, ..., r(N+1)=rN.
% kl, kt, lambda, mu contain parameters of coating layers 1..N.

    Nlayers = numel(r) - 1;
    P = eye(4);

    for j = 1:Nlayers
        Tleft = Tmat(n, r(j),   kl(j), kt(j), lambda(j), mu(j));
        Tright = Tmat(n, r(j+1), kl(j), kt(j), lambda(j), mu(j));
        Pj = Tright / Tleft;
        P = Pj * P;
    end
end
