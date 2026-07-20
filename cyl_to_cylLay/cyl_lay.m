
x = 3.0;              
R0 = 2.0;
R1 = 1.0;%  радиус центра
Nmaple = 3;          

rho_ext = 5;
c_ext = 148500.0;

rho_core = 2.7;
lambda_core = 5.3e11;
mu_core = 2.6e11;

rho_base = 2.7;
lambda_base = 5.3e11;
mu_base = 2.6e11;

CRO = [1.0, 0.0, 0.0];
CLA = [1e12, 0.0, 0.0];
CMU = [1e12, 0.0, 0.0];

Nlayers = Nmaple - 1;
r0 = R1;
rN = R0;
r = linspace(r0, rN, Nlayers + 1);

omega = x/R0 * c_ext;
k = omega / c_ext;

cl_core = sqrt((lambda_core + 2*mu_core) / rho_core);
ct_core = sqrt(mu_core / rho_core);
kl_core = omega / cl_core;
kt_core = omega / ct_core;

rho = zeros(1, Nlayers);
lambda = zeros(1, Nlayers);
mu = zeros(1, Nlayers);
kl = zeros(1, Nlayers);
kt = zeros(1, Nlayers);

for j = 1:Nlayers
    r_mid = 0.5*(r(j) + r(j+1));

    rho(j) = rho_base * polyval(fliplr(CRO), r_mid);
    lambda(j) = lambda_base * polyval(fliplr(CLA), r_mid);
    mu(j) = mu_base * polyval(fliplr(CMU), r_mid);

    cl = sqrt((lambda(j) + 2*mu(j)) / rho(j));
    ct = sqrt(mu(j) / rho(j));

    kl(j) = omega / cl;
    kt(j) = omega / ct;
end

A = 1;
rs = 2*rN;
phi_s = 0.0;


nmax = ceil(2*x + 2);
nvals = -nmax:nmax;
An = zeros(size(nvals));
Bcore = zeros(size(nvals));
Ccore = zeros(size(nvals));

for idx = 1:numel(nvals)
    n = nvals(idx);
    [An(idx), Bcore(idx), Ccore(idx)] = solveMode( ...
        n, r, r0, rN, kl, kt, lambda, mu, ...
        kl_core, kt_core, lambda_core, mu_core, ...
        omega, rho_ext, k, A, rs);
end


r_obs = 100*rN;
phi = linspace(0, 2*pi, 721);
Psi_s = zeros(size(phi));

for idx = 1:numel(nvals)
    n = nvals(idx);
    Psi_s = Psi_s + An(idx) * cyl('H', n, k*r_obs) ...
        .* exp(1i*n*(phi - phi_s));
end

%в дальней зоне
Ffar = zeros(size(phi));

for idx = 1:numel(nvals)
    n = nvals(idx);
    Ffar = Ffar + An(idx) * exp(1i*n*pi/2) ...
        .* exp(1i*n*(phi - phi_s));
end

D = abs(Ffar);
D = D / max(D);
X = cos(phi).*D;
Y = sin(phi).*D;
axis equal;
grid on;

plot(X,Y,'-.','LineWidth', 1.1,'Color','k');
hold on;
xlim([-0.7 1.1])
ylim([-0.5 0.5])
xlabel('x');
ylabel('y');
legend('rs = 50*rN','rs = 2*rN');
%figure;
%if exist('polarplot', 'file')
    %polarplot(phi, D, 'LineWidth', 1.5,'Color','k');
%else
    %polar(phi, D);
%end

%
%figure;
%plot(phi, abs(Psi_s), 'LineWidth', 1.5);
%xlim([0, 2*pi]);
%grid on;
%xlabel('phi');
%ylabel('|Psi_s|');
%

%figure;
%stem(nvals, abs(An), 'filled');
%grid on;
%xlabel('n');
%ylabel('|A_n|');
%title('Scattering coefficients');

function [A_n, B_core, C_core] = solveMode( ...
    n, r, r0, rN, kl, kt, lambda, mu, ...
    kl_core, kt_core, lambda_core, mu_core, ...
    omega, rho_ext, k, A, rs)

    P = Pmat(n, r, kl, kt, lambda, mu);
    T_0 = T0(n, r0, kl_core, kt_core, lambda_core, mu_core);
    M = P*T_0;

    Jn = cyl('J', n, k*rN);
    dJn_dr = k * dcyl('J', n, k*rN);
    Hn = cyl('H', n, k*rN);
    dHn_dr = k * dcyl('H', n, k*rN);
    Iinc = A * cyl('H', n, k*rs);

    L = [
        dHn_dr,            1i*omega*M(1,1), 1i*omega*M(1,2);
        1i*rho_ext*omega*Hn, M(3,1),        M(3,2);
        0,                 M(4,1),          M(4,2)
    ];

    rhs = [
        -Iinc * dJn_dr;
        -1i*rho_ext*omega * Iinc * Jn;
        0
    ];

    x = L \ rhs;

    A_n = x(1);
    B_core = x(2);
    C_core = x(3);
end
