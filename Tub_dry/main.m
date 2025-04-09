%m = msh('fort.14');
m = msh('mesh3.14');

%% Make a hump in the middle, in the x-direction
xs = m.p(:,1);
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
write(m, 'fv3', 'f14');

%% Test
c = deg2km(center) * 1000;
x1 = deg2km(min(xs)) * 1000;
x2 = deg2km(max(xs)) * 1000;

fun = @(x) exp(-(width*(x-c)) .^ 2) + 1;
fun2 = @(x) exp(-(width*(x-c)) .^ 2) + 1.001;

q = integral(fun, x1, x2);
q2 = integral(fun, x1, x2);