
function T0 = compute_T_0(r, k_l, k_t, lambda, mu, n)
    T_full = compute_T_j(r, k_l, k_t, lambda, mu, n);
    T0 = T_full(:, [1, 3]);   
end