%==========================================================================
DATA.geometry = 'Cartesian';         % 'Carteisan' or 'Polar'
DATA.geom1    = 60000;                   % x1 or r1
DATA.geom2    = 150000;              % x2 or r2    
DATA.h0       = 3;                   % Bathymetric coefficent
DATA.n        = 0;                   % Bathymetric power
DATA.amp      = 0.3;                 % Tidal amplitude of forcing
DATA.freq     = 2*pi*1/(12.4*3600);  % Tidal frequency of forcing (Note: this frequency corresponds to an M2 tide)
DATA.phase    = 0;                   % Tidal Phase of forcing
DATA.tau      = 0.0025;                   % Linear bottom friction factor
DATA.g        = 9.81;

T = 5*86400; % seconds
dt = 2;
noutge = 360;
frames = T / (dt*noutge);

%% 
coords = m.p;

solution = [];
ts = linspace(0,T,frames);

for i = 1:frames
    solution = [solution LG2D_Solutions(coords, i*dt*noutge, DATA)];
end
zeta = ncread("Mesh1/fort.63.dg.nc", "zeta");
eta = ncread("Mesh1/fort.63.adc.nc", "zeta");

zeta_end = zeta(:,end);

%%
pt = 28;
plot(ts, solution(pt,:), 'DisplayName','Analytical');
hold on
plot(ts ,eta(pt,:), 'DisplayName', 'ADCIRC');
plot(ts, zeta(pt,:), 'DisplayName','DG-CG');
hold off
legend