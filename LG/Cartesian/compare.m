%==========================================================================
DATA.geometry = 'Cartesian';         % 'Carteisan' or 'Polar'
DATA.geom1    = 60000;                   % x1 or r1
DATA.geom2    = 150000;              % x2 or r2    
DATA.h0       = 3;                   % Bathymetric coefficent
DATA.n        = 0;                   % Bathymetric power
DATA.amp      = 0.3;                 % Tidal amplitude of forcing
DATA.freq     = 2*pi*1/(12.4*3600);  % Tidal frequency of forcing (Note: this frequency corresponds to an M2 tide)
DATA.phase    = 0;                   % Tidal Phase of forcing
DATA.tau      = 0.005;                   % Linear bottom friction factor
DATA.g        = 9.81;

T = 5*86400; % seconds
dt = 1;
noutge = 360;
frames = T / (dt*noutge);
hs = [1875 3750 7500 15000];

%% Loop through different meshes
err_zeta = [];
err_eta = [];
ms = [1 2 3 4];
offset = 0;
for k = ms
    dataset = ['Mesh' num2str(k)];
    m = msh([dataset '/fort.14']);
    coords = m.p;
    N = size(m.p, 1);
    solution = [];
    ts = linspace(0,T,frames);
    
    for i = 1:frames
        [ze, u, v] = LG2D_Solutions(coords, i*dt*noutge, DATA);
        solution = [solution u];
    end
    zeta = ncread([dataset '/fort.64.dg.nc'], "u-vel");
    eta = ncread([dataset '/fort.64.adc.nc'], "u-vel");
    
    
    %% L2 norm at final snapshot
    err_zeta = [err_zeta sqrt((1/N) * sum( (zeta(:,end-offset)-solution(:,end-offset)).^2))];
    err_eta = [err_eta sqrt((1/N) * sum( (eta(:,end-offset)-solution(:,end-offset)).^2))];

end
err_zeta = fliplr(err_zeta);
loglog(hs(ms), err_zeta, '-s', 'DisplayName', 'DG-CG');
hold on
err_eta = fliplr(err_eta);
loglog(hs(ms), err_eta, '-o', 'DisplayName', 'CG');

loglog(hs(ms), hs(ms).^1 / 1e6, '--', 'DisplayName', 'h');
ylabel('L^2 error')
xlabel('h (m)')
hold off
legend
set(gca, 'XTick', hs)
%xtickformat('%.1e')


p1 = polyfit(log(hs(ms)), log(err_eta), 1);
p2 = polyfit(log(hs(ms)), log(err_zeta), 1);


%%
pt = round(N-15);
plot(ts, solution(pt,:), 'DisplayName','Analytical');
hold on
plot(ts ,eta(pt,:), 'DisplayName', 'ADCIRC');
plot(ts, zeta(pt,:), 'DisplayName','DG-CG');
hold off
legend