% This script solves the dynamics of the two segment leg with a hip and
% knee spring. It follows the lagrange method. It solves the ode using the
% ode45 function for small timesteps so that I can stop the solver based on
% solution info, specifically when the foot lifts off and stance is over.
% So basically, solve for intermediate time intervals and each time the
% solver finishes I check to see if the foot has lifted off yet, and if it
% has I finish solving, and if it hasn't I solve for another interval. 

% The definition of the model is in the XY plane, CCW rotation is positive,
% as shown in the pdf document. I have included damping as a
% nonconservative force in the lagrange method. I also solve for a new K1
% hip stiffness value every time interval that I solve for. I use the
% solveForK1 function.
%% Setup
clear
% Set up the leg model parameters
desiredGRFAngle = 70*pi/180;
mh = 3;
m1 = 0.01;
m2 = 0.01;
L1 = 0.1;
L2 = 0.1;
rHipSphere = 0.1; % m
rLegCylindar = 0.02;
Ih = (2/5)*mh*rHipSphere^2; %0.1;
I1 = (1/4)*m1*rLegCylindar^2 + (1/12)*m1*L1^2; %0.01;
I2 = (1/4)*m1*rLegCylindar^2 + (1/12)*m2*L2^2; %0.01;
theta1home = -45*pi/180;
theta2home = -90*pi/180;
g = 9.81;
K1 = 100;
K2 = 100;
Kg = 1e6;
c = 1000;

syms xh(t) yh(t) theta1(t) theta2(t)
% Set up my leg position equations
x1 = xh + (L1/2)*cos(theta1);
y1 = yh + (L1/2)*sin(theta1);
x2 = xh + L1*cos(theta1) + (L2/2)*cos(theta1+theta2);
y2 = yh + L1*sin(theta1) + (L2/2)*sin(theta1+theta2);
xf = xh + L1*cos(theta1) + L2*cos(theta1+theta2);
yf = yh + L1*sin(theta1) + L2*sin(theta1+theta2);

% Set up my kinetic energy equation
Th = (1/2)*mh*(diff(xh,t)^2+diff(yh,t)^2) + (1/2)*Ih*diff(theta1)^2;
T1 = (1/2)*m1*(diff(x1,t)^2+diff(y1,t)^2) + (1/2)*I1*diff(theta1)^2;
T2 = (1/2)*m2*(diff(x2,t)^2+diff(y2,t)^2) + (1/2)*I2*(diff(theta2,t)+diff(theta1,t))^2;
T = Th+T1+T2;

% Set up my nonconservative forces
rb = [xh;yh]; % vector position of the floating base
r1 = rb+L1*[cos(theta1);sin(theta1)]; % vector position of the end of link 1
r2 = r1+L2*[cos(theta1+theta2);sin(theta1+theta2)]; % vector position of the end of link 2
r2dot = diff(r2,t);
Fdamping = formula(-c*r2dot); % the formula() function is so that it makes a vector of sym variables instead of a sym function with a vector inside it

Qx = Fdamping(1); % Fd dot partial of r2 w.r.t. xh. the partial is just [1;0]
Qy = Fdamping(2); % Fd dot partial of r2 w.r.t. yh. the partial is just [0;1]

dr2dtheta1 = formula(diff(r2,theta1));
Qtheta1 = dr2dtheta1(1)*Fdamping(1)+dr2dtheta1(2)*Fdamping(2); % Fd dot partial of r2 w.r.t.theta1
dr2dtheta2 = formula(diff(r2,theta2));
Qtheta2 = dr2dtheta2(1)*Fdamping(1)+dr2dtheta2(2)*Fdamping(2); % Fd dot partial of r2 w.r.t.theta2


% Set up my initial simulation values 
angleOfTDwrtVertical = 30; % deg
theta1_array = (-45+angleOfTDwrtVertical)*pi/180;
dtheta1_array = 0;
xh_array = -sqrt(L1^2+L2^2)*sind(angleOfTDwrtVertical);
dxh_array = 5;
yh_array = sqrt(L1^2+L2^2)*cosd(angleOfTDwrtVertical);
dyh_array = -3;
theta2_array = -90*pi/180;
dtheta2_array = 0;
K1_array = K1;

timestep = 0.001; % seconds
t_array = 0;

%% Run the numerical simulation
for i = 1:100 
    % Find the current foot position as it presses into the floor contact,
    % modeled as spring dampers in the X and Y directions.
    current_xk = xh_array(end) + L1*cos(theta1_array(end));
    current_yk = yh_array(end) + L1*sin(theta1_array(end));
    current_xf = current_xk + L2*cos(theta1_array(end)+theta2_array(end));
    current_yf = current_yk + L2*sin(theta1_array(end)+theta2_array(end));

    % Solve for which K1 value to use now, given the current state of the
    % system.
    K1 = solveForK1(L1,L2,K2,theta1home,theta2home,desiredGRFAngle,xh_array(end)-current_xf,yh_array(end)-current_yf,theta2_array(end));

    % Error handling
    if isempty(K1)
        K1 = 0;
    end
    if K1<0
        K1=0;
    end

    % Set up my potential energy equation
    V = g*(mh*yh + m1*y1 + m2*y2) + (1/2)*K1*(theta1 - theta1home)^2 + (1/2)*K2*(theta2 - theta2home)^2 ...
        + (1/2)*Kg*(xf^2+yf^2);

    % Set up my Lagrange Equation
    L = T - V;
    
    % Set up my dynamic equations in the x, y, theta1, and theta2
    % directions.
    eq_x = Qx == diff(diff(L,diff(xh,t)),t) - diff(L,xh);
    eq_y = Qy == diff(diff(L,diff(yh,t)),t) - diff(L,yh);
    eq_theta1 = Qtheta1 == diff(diff(L,diff(theta1,t)),t) - diff(L,theta1);
    eq_theta2 = Qtheta2 == diff(diff(L,diff(theta2,t)),t) - diff(L,theta2);
    eq = [eq_x, eq_y, eq_theta1, eq_theta2];

    [F,S] = odeToVectorField(eq);
    M = matlabFunction(F, 'vars', {'t', 'Y'});

    % Set up initial conditions of the next timestep
    init_conds = [theta1_array(end) dtheta1_array(end) xh_array(end) dxh_array(end) theta2_array(end) dtheta2_array(end) yh_array(end) dyh_array(end)];
    t_span = [t_array(end) t_array(end)+timestep]; % Simulate for one timestep

    % Solve the ODE
    [t_sol, Y_sol] = ode45(@(t,Y) M(t,Y), t_span, init_conds);

    % Compile the results
    theta1_array  = [theta1_array; Y_sol(:,1)];
    dtheta1_array = [dtheta1_array; Y_sol(:,2)];
    xh_array      = [xh_array; Y_sol(:,3)];
    dxh_array     = [dxh_array; Y_sol(:,4)];
    theta2_array  = [theta2_array; Y_sol(:,5)];
    dtheta2_array = [dtheta2_array; Y_sol(:,6)];
    yh_array      = [yh_array; Y_sol(:,7)];
    dyh_array     = [dyh_array; Y_sol(:,8)];
    K1_array      = [K1_array; ones(size(t_sol))*K1];

    t_array = [t_array; t_sol];

    % If the current distance from the foot to the hip is larger than the
    % resting length from the foot to the hip then the leg is extended and
    % is lifting off the ground, and the simulation can end.
    if sqrt((xh_array(end)-current_xf)^2+(yh_array(end)-current_yf)^2)>sqrt(L1^2+L2^2)
        break
    end
end
%% GRF Plotting
% Plot the ground reaction force through the stance.
xk = xh_array + L1*cos(theta1_array);
yk = yh_array + L1*sin(theta1_array);
xf = xk + L2*cos(theta1_array+theta2_array);
yf = yk + L2*sin(theta1_array+theta2_array);
GRFx = -Kg*xf;
GRFy = -Kg*yf;
figure()
subplot(2,1,1)
p = plot(t_array/t_array(end)*100,atan2d(GRFy,GRFx),'LineWidth',2);
ylabel('GRF Direction (deg)','Interpreter','latex')
set(gca, 'LineWidth', 1.5);
set(gca, 'FontSize', 14);
subplot(2,1,2)
plot(t_array/t_array(end)*100,K1_array,'LineWidth',2)
ylabel('$K_1$','Interpreter','latex')
xlabel('Portion of Stance (\%)','Interpreter','latex')
set(gca, 'LineWidth', 1.5);
set(gca, 'FontSize', 14);

%% Leg Motion Replay Plotter
% This section plots a replay of the motion of the leg, showing the
% position of the leg at each timestep and a line showing the path the hip
% location took through space as the stance progressed.
% This section is purely for user understandability of what transpired.
for i = 1:length(theta1_array)
    if mod(i,100)==0 || i==length(theta1_array)% this line just speeds up the plotting process if there are a lot of data points in the solution
        % Generate hip, knee, and foot positions in the real/imaginary
        % plane for plotting purposes.
        hip = xh_array(i)+yh_array(i)*1i;

        xk = xh_array(i)+L1*cos(theta1_array(i));
        yk = yh_array(i)+L1*sin(theta1_array(i));
        knee = xk+yk*1i;

        xf = xk+L2*cos(theta1_array(i)+theta2_array(i));
        yf = yk+L2*sin(theta1_array(i)+theta2_array(i));
        foot = xf + yf*1i;

        % These are the x and y components of the trajectory the hip
        % follows through space.
        xpath = xh_array(1:i);
        ypath = yh_array(1:i);

        % Plot the instantaneous configuration of the leg, along with the
        % trajectory the hip has followed.
        figure(2)
        plot([foot knee hip])
        hold on
        plot(hip,'bo')
        plot(0,0,'rx')
        plot(xpath,ypath)
        drawnow;
        hold off
    end
end
figure(2)
title('Leg Motion Replay')