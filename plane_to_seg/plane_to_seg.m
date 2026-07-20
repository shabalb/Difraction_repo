clear;
clc;
% ==============================================================
%  Параметры задачи 
% ==============================================================
r0 = 1.0;              % радиус цилиндра
r1 = 1.5;              % внешний радиус покрытия
Nseg = 2;              
phi = linspace(0, 2*pi, Nseg+1);   % угловые границы
phi_s = 0.0;          

% Материалы сегментов 
rho_j = [ 11, 2.7, 7, 7 ];       
lambda_j = [ 0.5e11, 5.3e11, 10e11, 10e11 ];    
mu_j = [ 0.14e11, 2.6e11, 8.2e11, 8.2e11 ];

%lambda_j = [ 0.5e11, 1, 1, 1 ];    
%mu_j = [ 0.14e11, 1, 1, 1 ];

%rho_j = [ 1 ];
%lambda_j = [ 1e10 ];    
%mu_j = [ 1e10 ];


rho0 = 2.7;
lambda0 = 5.3e11;
mu0 = 2.6e11;

rho_f = 5;
c_f = 148500.0;

x = 3.0;
omega = x/r1 * c_f;
k = omega / c_f;

Nmax = 5;              
n_vals = -Nmax : Nmax;
M = length(n_vals);   

% ==============================================================
%  Волновые числа
% ==============================================================
% цилиндр
c_l0 = sqrt((lambda0 + 2*mu0) / rho0);
c_t0 = sqrt(mu0 / rho0);
k_l0 = omega / c_l0;
k_t0 = omega / c_t0;

% сегмент
c_l_j = sqrt((lambda_j + 2*mu_j) ./ rho_j);
c_t_j = sqrt(mu_j ./ rho_j);
k_l_j = omega ./ c_l_j;
k_t_j = omega ./ c_t_j;

k_f = omega / c_f;

I_mn = zeros(M, M, Nseg);  
for j = 1:Nseg
    phi_j = phi(j);
    phi_jp1 = phi(j+1);
    delta = phi_jp1 - phi_j;
    for idx_m = 1:M
        m = n_vals(idx_m);
        for idx_n = 1:M
            n = n_vals(idx_n);
            if n == m
                I_mn(idx_m, idx_n, j) = delta;
            else
                I_mn(idx_m, idx_n, j) = ( exp(1i*(n-m)*(phi_jp1 - phi_s)) - exp(1i*(n-m)*(phi_j - phi_s)) ) / (1i*(n-m));
            end
        end
    end
end

% Общее число неизвестных:
n_unknowns = Nseg * M * 4 + M * 2 + M;

nrows = 4 * Nseg * M;
A_inner = sparse(nrows, n_unknowns);

row = 0;
%A_inner = zeros((Nseg-1)*M*4 + (M-1)*4 + 4,(Nseg-1)*M*4 + (M-1)*4 + 4);
for j = 1:Nseg
    % параметры сегмента
    lambda = lambda_j(j);
    mu = mu_j(j);
    k_l = k_l_j(j);
    k_t = k_t_j(j);
    
    for idx_m = 1:M
        m = n_vals(idx_m);
        for comp = 1:4   % компонента вектора S = (u_r, u_phi, sigma_rr, sigma_rphi)
            row = (j-1)*M*4 + (idx_m-1)*4 + comp;
            
            for idx_n = 1:M
                n = n_vals(idx_n);
                Tj = Tmat(n, r0, k_l, k_t, lambda, mu);
                T_0 = T0(n, r0, k_l0, k_t0, lambda0, mu0);
                factor = I_mn(idx_m, idx_n, j);
                % вклад от сегмента
                for comp_k = 1:4
                    col = (j-1)*M*4 + (idx_n-1)*4 + comp_k;
                    A_inner(row, col) = A_inner(row, col) + Tj(comp, comp_k) * factor;
                end
                % вклад от цилиндра
                for comp_k = 1:2
                    col = Nseg*M*4 + (idx_n-1)*2 + comp_k;
                    A_inner(row, col) = A_inner(row, col) - T_0(comp, comp_k) * factor;
                end
            end
        end
    end
end

b_inner = zeros(nrows, 1);

% ---- Внешняя граница ----
nrows_outer = 3 * Nseg * M;
A_outer = sparse(nrows_outer, n_unknowns);
b_outer = zeros(nrows_outer, 1);

row = 0;
for j = 1:Nseg
    lambda = lambda_j(j);
    mu = mu_j(j);
    k_l = k_l_j(j);
    k_t = k_t_j(j);
    
    for idx_m = 1:M
        m = n_vals(idx_m);
        for comp = 1:3   % u_r sigma_rr sigma_rphi
            row = row + 1;
            for idx_n = 1:M
                n = n_vals(idx_n);
                factor = I_mn(idx_m, idx_n, j);
                
                Dj = compute_D_j(r1, k_l, k_t, lambda, mu, n);
                for comp_k = 1:4
                    col = (j-1)*M*4 + (idx_n-1)*4 + comp_k;
                    A_outer(row, col) = A_outer(row, col) + Dj(comp, comp_k) * factor;
                end
                
                if comp == 1
                    val_A = dcyl('H', n,  k_f * r1);
                elseif comp == 2
                    val_A = -1i * rho_f * omega * besselh(n, 1, k_f * r1);
                else 
                    val_A = 0;
                end
                col_A = Nseg*M*4 + M*2 + idx_n;  
                A_outer(row, col_A) = A_outer(row, col_A) - val_A * factor;
            end
            
            for idx_n = 1:M
                n = n_vals(idx_n);
                if comp == 1
                    val_inc = (1i)^n * k_f * 0.5*(besselj(n-1, k_f*r1) - besselj(n+1, k_f*r1)); 
                elseif comp == 2
                    val_inc = (1i)^n * (-1i * rho_f * omega * besselj(n, k_f*r1));
                else
                    val_inc = 0;
                end
                b_outer(row) = b_outer(row) + val_inc * I_mn(idx_m, idx_n, j);
            end
        end
    end
end

% Число точек коллокации
Nr = 2*M; 
% узлы и веса на [-1,1]
[x_gauss, w_gauss] = lgwt(Nr,r0,r1); 

r_colloc = 0.5*(r1 - r0) * x_gauss + 0.5*(r1 + r0);

nrows_colloc = 4 * Nseg * Nr;  
A_colloc = sparse(nrows_colloc, n_unknowns);
b_colloc = zeros(nrows_colloc, 1);

row = 0;
for j = 1:Nseg
    
    j_next = mod(j, Nseg) + 1;  
    
    
    Lambda_j = lambda_j(j); Mu_j = mu_j(j); K_l_j = k_l_j(j); K_t_j = k_t_j(j);
    Lambda_jn = lambda_j(j_next);
    Mu_jn = mu_j(j_next);
    K_l_jn = k_l_j(j_next);
    K_t_jn = k_t_j(j_next);
    
    phi_j = phi(j);  
    
    for p = 1:Nr
        rp = r_colloc(p);
        for comp = 1:4 
            row = row + 1;
            for idx_n = 1:M
                n = n_vals(idx_n);
                phase = exp(1i * n * (phi_j - phi_s));
                
               
                Wj = compute_W_j(rp, K_l_j, K_t_j, Lambda_j, Mu_j, n);
               
                Wjn = compute_W_j(rp, K_l_jn, K_t_jn, Lambda_jn, Mu_jn, n);
                
                
                for comp_k = 1:4
                    col = (j-1)*M*4 + (idx_n-1)*4 + comp_k;
                    A_colloc(row, col) = A_colloc(row, col) + Wj(comp, comp_k) * phase;
                end
                
                for comp_k = 1:4
                    col = (j_next-1)*M*4 + (idx_n-1)*4 + comp_k;
                    A_colloc(row, col) = A_colloc(row, col) - Wjn(comp, comp_k) * phase;
                end
            end
            
        end
    end
end

% общая матрица и правая часть
A_total = [A_inner; A_outer; A_colloc];
b_total = [b_inner; b_outer; b_colloc];


X = A_total \ b_total;

K_seg = reshape(X(1 : Nseg*M*4), [4, M, Nseg]);
K0 = reshape(X(Nseg*M*4 + 1 : Nseg*M*4 + M*2), [2, M]);
A_n = X(Nseg*M*4 + M*2 + 1 : end);

rN = r1;
r_obs = 100*rN;
phi = linspace(0, 2*pi, 721);
Psi_s = zeros(size(phi));

for idx = 1:numel(n_vals)
    n = n_vals(idx);
    Psi_s = Psi_s + A_n(idx) * cyl('H', n, k*r_obs) ...
        .* exp(1i*n*(phi - phi_s));
end

Ffar = zeros(size(phi));

for idx = 1:numel(n_vals)
    n = n_vals(idx);
    Ffar = Ffar + A_n(idx) * exp(1i*n*pi/2) ...
        .* exp(1i*n*(phi - phi_s));
end

D = abs(Ffar);
D = D / max(D);
X = cos(phi).*D;
Y = sin(phi).*D;
%axis equal;

%,'-','LineWidth', 1.1,'Color','k'
polarplot(X,Y,'-','LineWidth', 1.1,'Color','k');
%grid on;
%hold on;
rlim([0 0.7])
%xlim([0 2])
%ylim([-2 2])
%xlabel('x');
%ylabel('y');