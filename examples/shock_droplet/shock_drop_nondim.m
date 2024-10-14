% Lee    Chase

% Shock drop input file
% Nonlinear system of equations to get relevant pre and post shock
% conditions to complete initialization. Given that post shock density and
% velocity are set to unity, the only changes to the case would be the gas
% Mach number (Ma_g) or the specific heat ratio (gamma). 

% Run this matlab script first with imposed gamma and Ma_g to get the shock
% Mach number, post shock pressure, pre shock density, and pre shock
% pressure. Once the shock Mach number is found and put into the input
% file, all other values will automatically be computed.


% Inputs:

% rho_1 : post shock density
% v_shock: post shock velocity
% Ma_s: gas Mach number

% Outputs:

% x(1): Shock Mach number
% x(2): post shock pressure
% x(3): pre shock density
% x(4); pre shock pressure

rho_1 = 1; v_shock = 1; gamma=1.4; Ma_g = 0.577;

fun = @(x)findval(x,gamma,rho_1,v_shock,Ma_s);
x = fsolve(fun,[1.01 1 1.5 1]);


%% fsolve function

function F = findval(x,gamma,rho,v,Ma_1)
% x(1) = Ma_s
% x(2) = p_1
% x(3) = rho_0
% x(4) = p_0

F(1) = x(4)*(2*gamma*x(1).^2 - (gamma - 1))/(gamma+1) - x(2);
F(2) = x(3)*(x(1).^2*(gamma+1))/((gamma-1)*x(1).^2+2) - rho;
F(3) = rho*v/(gamma*Ma_1.^2) - x(2); 
F(4) = -sqrt(((gamma-1)*x(1).^2+2)/(2*gamma*x(1)^2-(gamma-1)))*sqrt(gamma*x(2)/rho) + x(1)*sqrt(gamma*x(4)/x(3)) - 1;

end