close all ; clc ; clear

%% Gravity-gradient stability simulation 
a = 7000;   % km
mu = 398600.4418; % km^3/s^2
euler_0 = deg2rad([30, 20, 10]); % rad
euler_dot_0 = [0,0,0]; % rad/s
x_0 = [euler_0 euler_dot_0];
I_big = 400; % kgm^2
I_inter = 300; % kgm^2
I_small = 200; % kgm^2
t_span = 0:0.1:10000; % s
n = sqrt(mu / a^3);

I_1 = [I_big 0 0;
       0 I_inter 0;
       0 0 I_small];
I_2 = [I_inter 0 0;
       0 I_big 0;
       0 0 I_small];
I_3 = [I_small 0 0;
       0 I_inter 0;
       0 0 I_big];

for i = 1:3
    switch i
        case 1
            I = I_1;
        case 2
            I = I_2;
        case 3
            I = I_3;
    end

    [~, x] = ode45(@(t,x) gravity_gradient(x,I,n), t_span, x_0);
    
    x(:, 1:3) = mod(x(:, 1:3) + pi, 2*pi) - pi;

    figure()
    overallTitle = sprintf('Gravity Gradient Satellite Dynamics: MOI Matrix %i', i);
    sgtitle(overallTitle);

    hold on;
    subplot(3,1,1)
    plot(t_span, x(:, 1));
    xlabel('Time (s)')
    ylabel('\phi (rad)')
    title('Roll Angle, \phi')
    
    hold on;
    subplot(3,1,2)
    plot(t_span, x(:, 2));
    xlabel('Time (s)')
    ylabel('\theta (rad)')
    title('Pitch Angle, \theta')
    
    hold on;
    subplot(3,1,3)
    plot(t_span, x(:, 3));
    xlabel('Time (s)')
    ylabel('\psi (rad)')
    title('Yaw Angle, \psi')
    
    % hold on;
    % subplot(2,3,4)
    % plot(t_span, x(:, 4));
    % xlabel('Time (s)')
    % ylabel('\phi Dot (rad/s)')
    % title('Roll Rate, \phi Dot')
    % 
    % hold on;
    % subplot(2,3,5)
    % plot(t_span, x(:, 5));
    % xlabel('Time (s)')
    % ylabel('\theta Dot (rad/s)')
    % title('Pitch Rate, \theta Dot')
    % 
    % hold on;
    % subplot(2,3,6)
    % plot(t_span, x(:, 6));
    % xlabel('Time (s)')
    % ylabel('\psi Dot (rad/s)')
    % title('Yaw Rate, \psi Dot')
end

function x_dot = gravity_gradient(x,I,n)
    ky = (I(2,2)-I(1,1))/I(3,3);
    kr = (I(2,2)-I(3,3))/I(1,1);
    A = [0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         -4*n^2*ky 0 0 0 0 -n*(1-ky);
         0 -3*n^2*((I(1,1)-I(2,2))/I(3,3)) 0 0 0 0;
         0 0 -n^2*kr -n*(kr-1) 0 0];
    x_dot = A * x;
end