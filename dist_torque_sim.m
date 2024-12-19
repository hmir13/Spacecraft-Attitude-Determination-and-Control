close all ; clc ; clear

% Disturbance Torque Simulation
%% Problem Constants
% Earth constants
mu_e = 398600.4418; % km^3/s^2
R_e = 6378; % km
w_e = 7.2921159e-5; % rad/s

% initial state
r_0 = [6930; % km
       0;
       0];
v_0 = [0;    % km/s
       5.3895;
       5.3895];
w_0 = [0;    % rad/s
       0;
       0];

x_0 = [r_0;
       v_0];

% satellite model
I = [50 0   0;  % kg-m^2
     0  100 0;
     0  0   50];
x = 1.7; % m
y = 1;   % m
z = 1.7; % m

r_cp_1 = 0.2; % m
r_cp_2 = 0.1; % m
r_cp_3 = 0; % m

A_1 = y * z; % m^2
A_2 = x * z; % m^2
A_3 = x * y; % m^2

t_0 = 0; % s
dt = 1; % s
t_f = 3.5 * 3600; % s
t_span = t_0:dt:t_f; % s

C_d = 2.2;
P_sol = 4.644e-6; % N/m^2
rho_s = 0.55;
rho_d = 0.2;
B_0 = 3.12e-5; % T
d = [1;     % A-m^2
     1.2;
     0.8];

%% Orbit Propagation
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[~,x] = ode45(@(t,x) keplers_2body(t,x,mu_e), t_span, x_0, opts);

X = x(:,1);
Y = x(:,2);
Z = x(:,3);
X_dot = x(:,4);
Y_dot = x(:,5);
Z_dot = x(:,6);

r_eci = [X, Y, Z];

figure()
sgtitle("Inertial position vector components over time")

subplot(3,1,1)
plot(t_span, X, DisplayName="X(t)")
ylabel("X [km]")

subplot(3,1,2)
plot(t_span, Y, DisplayName="Y(t)")
ylabel("Y [km]")

subplot(3,1,3)
plot(t_span, Z, DisplayName="Z(t)")
xlabel("time [s]")
ylabel("Z [km]")

figure()
plot3(X, Y, Z)
title("3D Orbit of satellite in inertial frame")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
grid

%% Aerodynamic Torque
v_rel = 1000 * [X_dot + w_e * Y Y_dot - w_e * X Z_dot]; % m/s

r_eci_mag = sqrt(X.^2 + Y.^2 + Z.^2);
h = r_eci_mag - R_e;
T = -131.21 + 2.99 * h;
P = 2.488 * ((T + 273.1) / 216.6).^-11.388;
rho = P ./ (0.2869 * (T + 273.1));

F_aero_1 = (1/2) * rho * C_d * (v_rel(1)^2 + v_rel(2)^2 + v_rel(3)^2) * A_1;
F_aero_2 = (1/2) * rho * C_d * (v_rel(1)^2 + v_rel(2)^2 + v_rel(3)^2) * A_2;
F_aero_3 = (1/2) * rho * C_d * (v_rel(1)^2 + v_rel(2)^2 + v_rel(3)^2) * A_3;

T_aero_1 = r_cp_1 * F_aero_1;
T_aero_2 = r_cp_2 * F_aero_2;
T_aero_3 = r_cp_3 * F_aero_3;

figure()
sgtitle("Aerodynamic torque over time")

subplot(3,1,1)
plot(t_span, T_aero_1)
ylabel("\tau_1 [N-m]")

subplot(3,1,2)
plot(t_span, T_aero_2)
ylabel("\tau_2 [N-m]")

subplot(3,1,3)
plot(t_span, T_aero_3)
xlabel("time [s]")
ylabel("\tau_3 [N-m]")

fprintf("Altitude, h(t_0) [km]: %d\n", h(1))
fprintf("Atmospheric density, rho(t_0) [kg/m^3]: %d\n", rho(1))
fprintf("Relative velocity magnitude, v_rel(t_0) [m/s]: %d\n", v_rel(1))
fprintf("\n")

%% Solar Radiation Pressure Torque
F_srp_1 = P_sol * A_1 * (1 + rho_s + (2/3) * rho_d);
F_srp_2 = P_sol * A_2 * (1 + rho_s + (2/3) * rho_d);
F_srp_3 = P_sol * A_3 * (1 + rho_s + (2/3) * rho_d);

T_srp_1 = r_cp_1 * F_srp_1;
T_srp_2 = r_cp_2 * F_srp_2;
T_srp_3 = r_cp_3 * F_srp_3;

% convert single value torque to constant torque array over time
T_srp_1 = zeros(size(t_span)) + T_srp_1;
T_srp_2 = zeros(size(t_span)) + T_srp_2;
T_srp_3 = zeros(size(t_span)) + T_srp_3;

figure()
sgtitle("Solar radiation pressure torque over time")

subplot(3,1,1)
plot(t_span, T_srp_1)
ylabel("\tau_1 [N-m]")

subplot(3,1,2)
plot(t_span, T_srp_2)
ylabel("\tau_2 [N-m]")

subplot(3,1,3)
plot(t_span, T_srp_3)
xlabel("time [s]")
ylabel("\tau_3 [N-m]")

fprintf("T_1^srp(t_0), [N-m]: %d\n", T_srp_1(1))
fprintf("T_2^srp(t_0), [N-m]: %d\n", T_srp_2(1))
fprintf("T_3^srp(t_0), [N-m]: %d\n", T_srp_3(1))
fprintf("\n")

%% Magnetic Field Torque
r_ecef = zeros(3, length(t_span));

for i = 1:length(t_span)
    t = t_span(i);
    dcm_eci_to_ecef = [ cos(w_e*t), sin(w_e*t), 0;
                       -sin(w_e*t), cos(w_e*t), 0;
                        0,          0,          1];
    r_ecef(:, i) = dcm_eci_to_ecef * r_eci(i,:)';
end

r_ecef_mag = sqrt(r_ecef(1).^2 + r_ecef(2).^2 + r_ecef(3).^2);

lambda = atan2(r_ecef(3,:), sqrt(r_ecef(1,:).^2 + r_ecef(2,:).^2));

B = B_0 * (R_e / r_ecef_mag).^3 * sqrt(1 + 3 * cos(lambda).^2);

T_mag_1 = d(1) * B;
T_mag_2 = d(2) * B;
T_mag_3 = d(3) * B;

figure()
sgtitle("Magnetic field torque over time")

subplot(3,1,1)
plot(t_span, T_mag_1)
ylabel("\tau_1 [N-m]")

subplot(3,1,2)
plot(t_span, T_mag_2)
ylabel("\tau_2 [N-m]")

subplot(3,1,3)
plot(t_span, T_mag_3)
xlabel("time [s]")
ylabel("\tau_3 [N-m]")

fprintf("Latitude, lambda(t_0), [deg]: %d\n", rad2deg(lambda(1)));
fprintf("Magnetic field intensity, B(t_0), [T]: %d\n", B(1));
fprintf("\n")
    
%% Gravity-Gradient Torque
R_c = r_eci;
R_c_mag = r_eci_mag;
a = (3 * mu_e ./ R_c_mag.^5);

T_g_1 = a .* (R_c(:,2) .* R_c(:,3) * (I(3,3) - I(2,2)));
T_g_2 = a .* (R_c(:,1) .* R_c(:,3) * (I(1,1) - I(3,3)));
T_g_3 = a .* (R_c(:,1) .* R_c(:,2) * (I(2,2) - I(1,1)));

figure()
sgtitle("Gravity-gradient torque over time")

subplot(3,1,1)
plot(t_span, T_g_1)
ylabel("\tau_1 [N-m]")

subplot(3,1,2)
plot(t_span, T_g_2)
ylabel("\tau_2 [N-m]")

subplot(3,1,3)
plot(t_span, T_g_3)
xlabel("time [s]")
ylabel("\tau_3 [N-m]")

%% Simulation of uncontrolled attitude dynamics
T_aero_1 = T_aero_1';
T_aero_2 = T_aero_2';
T_aero_3 = T_aero_3';

T_g_1 = T_g_1';
T_g_2 = T_g_2';
T_g_3 = T_g_3';

T_1 = T_aero_1 + T_srp_1 + T_mag_1 + T_g_1;
T_2 = T_aero_2 + T_srp_2 + T_mag_2 + T_g_2;
T_3 = T_aero_3 + T_srp_3 + T_mag_3 + T_g_3;

figure()
sgtitle("Total and disturbance torques over time")

subplot(3,1,1)
hold on
plot(t_span, T_1, DisplayName="Total")
plot(t_span, T_aero_1, DisplayName="Aero")
plot(t_span, T_srp_1, DisplayName="SRP")
plot(t_span, T_mag_1, DisplayName="Mag")
plot(t_span, T_g_1, DisplayName="Grav")
ylabel("\tau_1 [N-m]")
legend()
hold off

subplot(3,1,2)
hold on
plot(t_span, T_2, DisplayName="Total")
plot(t_span, T_aero_2, DisplayName="Aero")
plot(t_span, T_srp_2, DisplayName="SRP")
plot(t_span, T_mag_2, DisplayName="Mag")
plot(t_span, T_g_2, DisplayName="Grav")
ylabel("\tau_2 [N-m]")
legend()
hold off

subplot(3,1,3)
hold on
plot(t_span, T_3, DisplayName="Total")
plot(t_span, T_aero_3, DisplayName="Aero")
plot(t_span, T_srp_3, DisplayName="SRP")
plot(t_span, T_mag_3, DisplayName="Mag")
plot(t_span, T_g_3, DisplayName="Grav")
xlabel("time [s]")
ylabel("\tau_3 [N-m]")
legend()
hold off

[~, w] = ode45(@(t,w) euler_rot_gen(w, I(1,1), I(2,2), I(3,3), T_1(find(abs(t_span - t) <= dt / 2, 1)), T_2(find(abs(t_span - t) <= dt / 2, 1)), T_3(find(abs(t_span - t) <= dt / 2, 1))), t_span, w_0);

figure()
sgtitle("Angular velocities over time")

subplot(3,1,1)
plot(t_span, w(:,1))
ylabel("\omega_1 [rad/s]")

subplot(3,1,2)
plot(t_span, w(:,2))
ylabel("\omega_2 [rad/s]")

subplot(3,1,3)
plot(t_span, w(:,3))
xlabel("time [s]")
ylabel("\omega_3 [rad/s]")

%% Functions
function dx = keplers_2body(~,x,mu)
    r = x(1:3);
    v = x(4:6);
 
    r_norm = sqrt(r(1)^2 + r(2)^2 + r(3)^2);

    a = -mu / r_norm^3 * r;

    dx = [v;a];
end

function w_dot = euler_rot_gen(w,I_1,I_2,I_3,T_1,T_2,T_3)
    w_dot = zeros([3,1]);
    w_dot(1) = ((I_2 - I_3) * w(2) * w(3) + T_1) / I_1;
    w_dot(2) = ((I_3 - I_1) * w(3) * w(1) + T_2) / I_2;
    w_dot(3) = ((I_1 - I_2) * w(1) * w(2) + T_3) / I_3;
end
