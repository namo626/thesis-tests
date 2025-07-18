m = msh('fort.14');
%m = msh('../LG/Cartesian/Mesh3/fort.14');

%% Make a hump in the middle, in the x-direction
xs = m.p(:,1);
ys = m.p(:,2);
center = 0.5 * (max(xs) + min(xs));
width = 1e-4;
%bathy = -exp(-width*(m.p(:,1)-center) .^ 2) + 5;


bathy = -exp(-(width*(m.p(:,1)-center)) .^ 2) - 1;

%% Remove boundary conditions
m.op = [];
m.bd = [];
m = make_bc(m, 'outer', 1);

%% set
m.b = bathy;
write(m, 'fv4', 'f14');

%% Test
c = deg2km(center) * 1000;
x1 = deg2km(min(xs)) * 1000;
x2 = deg2km(max(xs)) * 1000;

fun = @(x) exp(-(width*(x-c)) .^ 2) + 1;
fun2 = @(x) exp(-(width*(x-c)) .^ 2) + 1.001;

q = integral(fun, x1, x2);
q2 = integral(fun2, x1, x2);

y1 = deg2km(min(ys)) * 1000;
y2 = deg2km(max(ys)) * 1000;

area = (y2-y1) * (x2-x1); % m^2
hump_volume = q*(y2-y1); % m^3
time = 86400*2;
rain_volume = 7.0556e-6 * time * area; % m^3

xss = linspace(x1, x2, 1000);
yss = fun2(xss);
real_area = arclength(xss, yss) * (y2 - y1);
real_rain_volume = 7.0556e-6 * time * real_area; 

sim_volume = 2.4173 * real_area - hump_volume;
err = abs(sim_volume - real_rain_volume) / real_rain_volume;