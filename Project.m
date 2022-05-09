%% PROJECT 1 ME 3242
%
%
% Original Script:  Xinyu Zhao, 02/26/2021
% Modification:     Patrick Meagher, 4/18/2021
% Update:           Xinyu Zhao, 03/06/2022
% 
% Student:          Kevin Cayo
% Instructer:       Xinyu Zhao 
% Courese:          ME 3242
% Date:             Wed, April 6, 2022
%%
clc
clear

%% PARAMETERS / KNOWN VARIABLES / EQUATIONS

% Material A : Temparture, Length, Thermal Conductivity, & Heat Source
T_Al = 1700;           % [K], temp on the left boundary
L_A = 0.001;           % [m]
k_A = 0.8;             % [W/mK]
q_A = -100000;         % [W/m]

% Material B : Temparture, Length, Thermal Conductivity, & Heat Source
T_Br = 1500;           % [K], temp on the right boundary
L_B = 0.0008;          % [m]
k_B = 33.0;            % [W/mK]
q_B = 0.0;             % [W/m]

% Variables & Equations
L = L_A + L_B;         % [m], total domain length
N = 18;                 % total amount of subdomains; i-1 grid point
dx = L/N;              % size of the cells
k = zeros(N,1);        % array holds all the k values at each point
qdot = zeros(N,1);     % array holds all the heat source values at each point
i_mid = floor(L_A/dx); % Index of the last point in block A

for i = 1: i_mid
    k(i) = k_A;         % conductivity values from material A
    qdot(i) = q_A;      % heat source terms from material A
end
for i = (i_mid+1): N
    k(i) = k_B;         % conductivity values from material B
    qdot(i) = q_B;      %conductivity values from material B
end

%% SET UP THE LEFT-HAND SIDE of AX=b
A = zeros(N,N)        % initialize a N x N dimension matrix with zero entries

%First set up values for the interior points
for i = 2:N-1
    A(i,i-1) = k(i);             %coefficient in front of point i-1
    A(i,i) = -(k(i)+k(i+1));     %coefficient in front of point i
    A(i,i+1) = k(i+1);           %coefficient in front of point i+1
end

% For the left boundary k values 
A(1,1) = 3*k(1);     % coefficient of point 1
A(1,2) = -1*k(2);    % coefficient of point 2

% For the right boundary k values 
A(N,N) = 3*k(N);         % coefficient of point N
A(N,N-1) = -1*k(N-1);    % coefficient of point N-1

%% SET UP THE RIGHT-HAND SIDE of AX=b
b = zeros(N,1);

for i = 2: N-1
    b(i) = -qdot(i)*(dx*dx); % right-hand side values corresponding to point i
end

b(1) = 2*T_Al*k(1)+qdot(1)*(dx*dx);   % right-hand side values corresponding to point 1
b(N) = 2*T_Br*k(N)+qdot(N)*(dx*dx);   % right-hand side values corresponding to point N

%% FIND THE SOLUTION OF T_I AT EACH GRID POINT

X = A\b;
T = X';

%% Visualize the temperature at each grid point as a function of x coordinates
xgrid = 0.5*dx:dx:L-0.5*dx;
xgrid = [0,xgrid,L];
plot(xgrid,[T_Al,X',T_Br],'x-','Linewidth',1.25)
grid on

xlabel('x (m)')
ylabel('T (K)')
title('Project 1 : N = 18')
xlim([0 L])

fprintf('The temperature @ each grid point\n')
fprintf('  Grid point   Temperature Locataion          Temperature\n')

for i=1: N-1
    fprintf('  %6.0f       ', i)
    fprintf('    T%.0f          ',i)
    fprintf('             %f\n', T(i))
end

hold on


%% Calculate exact solution based on given parameters

% Using the following BCs

% T(0) = Tl
% T(LA+LB) = Tr
% T(LA-) = T(LA+)
% -kA dT/dx|- = -kB dT/dx|+

% We can construct a system of equations in the form of 
% Aexact*C = bexact
% Where C is a vector of the unknown coeffecients from the integration of
% the heat equation, A is a matrix of coeffeceints and b is a vector of
% constants. Using MATLAB we can solve for C with ease.

%Set coeffeceint matrix 
Aexact = [0, 1, 0, 0
          0, 0, L_A+L_B, 1
          L_A, 1, -L_A, -1
          k_A, 0, -k_B, 0]

% Set b vector
bexact = [T_Al, T_Br, q_A*L_A^2/(2*k_A), q_A*L_A]'

% Calculate constants of integration
C = Aexact\bexact

% Temperature distribution in A
T_A_exact = @(x) -q_A*x^2/(2*k_A) + C(1)*x + C(2)

% Temperature distribution in B
T_B_exact = @(x) C(3)*x + C(4)

% Plot temperature in A
fplot(T_A_exact,[0, L_A])

% Plot Temperature in B
fplot(T_B_exact,[L_A, L_A+L_B])

% Add a legend
legend(sprintf('N=%i',N),'Exact, A','Exact, B')
legend boxoff
