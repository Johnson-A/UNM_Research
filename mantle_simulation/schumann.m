function schumann
%SCHUMANN reproduce schumann (1929) results computationally
%   Find magnitude of relative motion
%   Find kf st fluid temperature decreses by less than 10% over interval
%   Consider different ks (solid tracks fluid at what ks)
n = 400;
k = 1.0;
T_0 = 1.0;

nstreams = 10;

interval = 1.0;
x = linspace(0,interval,n);
region_1 = x < (interval / 2);
region_2 = x >= (interval / 2);

porosity = [0.1, 0.2];
porosityA(region_1) = porosity(1);
porosityA(region_2) = porosity(2);
k_solid = k ./ (1.5 * (1 - porosityA));
k_fluid = k ./ (1.0 * porosityA);

dx = interval / (n-1);

% Conserve fluid mass flux by changing velocity across the boundary
% Velocity change at boundary
v = interval / 1.0;

dt = dx / v / 500; % Stability Criterion

T_solid = zeros(1, n);
T_fluid = zeros(1, n);
T_solid_steps = [];
T_fluid_steps = [];
t_steps = [];

figure(1); hold on;
solid_plot = plot(x, T_solid);
fluid_plot = plot(x, T_fluid);
ylim([0,T_0]);
legend('Solid', 'Fluid');
ylabel('Temperature'); xlabel('x');

tic;

step = 0;
for t = 0:dt:10.0/2, step = step + 1;
    in = x <= v * t;
    
    % To eliminate numerical instabilities, impose a boundary condition on
    % all new points which have been reached since the last time step.
    % These points will be at a temperature according to equation (8).
    % Note that dt must be small enough that at most 1 new partition is
    % added every several steps
    new_points = in & x > v * (t - dt);
    if any(new_points & region_2)
        n_mid = find(region_2, 1,'first');
        t_half = x(n_mid) / v;
        T_fluid(new_points) = T_0 * exp(-k_fluid(1) * t_half + ...
                                        -k_fluid(n_mid) * (t - t_half));
    else
        T_fluid(new_points) = T_0 * exp(-k_fluid(1) * t);
    end
    
    dT_fluid_dx = gradient(T_fluid(in), dx);
    delta_T = T_fluid(in) - T_solid(in);
    
    dT_fluid_dt = -k_fluid(in) .* delta_T - v * dT_fluid_dx;
    dT_solid_dt = k_solid(in) .* delta_T;
    
    
    T_fluid(in) = T_fluid(in) + dT_fluid_dt * dt;
    T_solid(in) = T_solid(in) + dT_solid_dt * dt;
    T_fluid(1) = T_0;
    T_fluid(end/2+1) = T_fluid(end/2);
    T_solid(end/2+1) = T_solid(end/2);
    
    if toc > 1/30
        set(solid_plot, 'YDATA', T_solid);
        set(fluid_plot, 'YDATA', T_fluid);
        drawnow
        tic;
    end
    
    if mod(step, floor(1 / dt / nstreams)) == 0
        t
        t_steps = [t_steps; t];
        T_solid_steps = [T_solid_steps; T_solid(1:n/nstreams:end)];
        T_fluid_steps = [T_fluid_steps; T_fluid(1:n/nstreams:end)];
    end
end

x_points = x(1:n/nstreams:end);
time_since_contact = x_points / v;
y = k_fluid * time_since_contact;
% z = k_solid * (t_steps - time_since_contact);

figure(); hold on; title('Solid');
xlabel('Z = k_{solid} (t - x / v)'); ylabel('Ts / To');
ind = 1;
for stream = T_solid_steps
    z = k_solid * (t_steps - time_since_contact(ind));
    plot(z, stream)
    ind = ind + 1;
end
xlim([0,k_solid * t_steps(end)]);

figure(); hold on; title('Fluid');
xlabel('Z = k_{solid} (t - x / v)'); ylabel('Tg / To');
ind = 1;
for stream = T_fluid_steps
    z = k_solid * (t_steps - time_since_contact(ind));
    plot(z, stream)
    ind = ind + 1;
end
xlim([0,k_solid * t_steps(end)]);
end