function D = compute_D_j(r, k_l, k_t, lambda, mu, n)
   
    T = Tmat(n, r, k_l, k_t, lambda, mu);
    D = T([1,3,4], :);  
end