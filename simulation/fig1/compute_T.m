function T = compute_T(wsim, geo)

xs = geo.xs;
ys = geo.ys;
ng = geo.ng;
hx = geo.hx;
chi_s = geo.chi_s;

wl = 2*pi*1/wsim;
k = 2*pi/wl;

% *** build GF on scatter grid
R = zeros(ng,ng);
G = zeros(ng,ng);
for ir = 1:ng
    xx = xs - xs(ir);
    yy = ys - ys(ir);
    rr = sqrt(xx.^2+yy.^2);
    R(ir,:) = rr;
    G(ir,:) = 1j.*k.^2 ./ 4 .* besselh(0,k*rr);
    G(ir,ir) = 0;
end

quad_order = 6;
% *** D0 and D1
[D0_t,D1_t]= planck_win_get_D(k,hx);
D0=D0_t*hx^2; D1=D1_t*hx^4;
% *** Duan-Rokhlin correction
Tao = zeros(ng,ng);
switch quad_order
    case 0 % no correction, do nothing
    case 4 % forth order correction
        Tao(1:ng+1:end) =  D0.*ones(1,ng);
    case 6 % sixth order correction
        for ir1 = 1:ng
            for ir2 = 1:ng
                xx = xs(ir1) - xs(ir2);
                yy = ys(ir1) - ys(ir2);
                rr = sqrt(xx^2+yy^2);
                if abs(rr)<1e-14
                    Tao(ir1,ir2) = D0 - 2*D1/hx^2;
                elseif abs(rr-hx)<1e-14
                    Tao(ir1,ir2) = 1/2 * D1 / hx^2;
                end
            end
        end
end

G = G + 1j*k^2/4/hx^2 * Tao;

% *** material matrix
Xi_vec = -1./chi_s;
Xi = zeros(ng,ng);
Xi(1:ng+1:end) = Xi_vec;

% *** T matrix
T = -inv(G*hx^2 + Xi);

end