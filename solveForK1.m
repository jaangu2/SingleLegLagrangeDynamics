function Kp = solveForK1(L1, L2, K2, theta10, theta20, desiredGRFAngle, x_in, y_in, kneePosNegConfig)
% This function uses symbolic variables to generate the equations of force
% in the X and Y directions, and then uses that to solve for the necessary
% stiffness value to direct the GRF in the desired direction.

syms x y K1 real

% These equations transform XY coordinates into theta1 theta2 coords.
theta2 = acos((x^2+y^2-L1^2-L2^2)/(2*L1*L2));
theta2 = theta2*sign(kneePosNegConfig);
intermediateValue1 = L1+L2*cos(theta2);
intermediateValue2 = L2*sin(theta2);
theta1 = (atan2(y,x) - atan2(intermediateValue2,intermediateValue1));

% The potential energy of the system is based on the springs and the
% rotation of the leg segements.
E = (1/2)*K1*(theta1-theta10)^2 + (1/2)*K2*(theta2-theta20)^2;

% Force is the derivative of energy w.r.t. position
Fx = diff(E,x);
Fy = diff(E,y);

% I substitute negative values of x and y because the theta1 and theta2
% equations above that are based on x and y are based on x and y being the
% position of the end of link 2, like an end effector for a two link
% manipulator. The x and y that is being sent into this function is
% actually the x and y of the hip, so the "end effector" would actually be
% the negative of the hip x,y values.
finalFx = subs(Fx,[x,y],[-x_in,-y_in]); 
finalFy = subs(Fy,[x,y],[-x_in,-y_in]);

% This is the final equation that sets the forces equal to my desired
% direction, and then solves for the K1 necessary to make that happen.
eqn = tan(desiredGRFAngle-pi) == finalFy/finalFx;
mySolution = solve(eqn,K1); % symvar pulls out all the symvars in the function and thats how the solve function knows to solve for both Ka and Kg
% This is the final result
Kp = double(mySolution);

end