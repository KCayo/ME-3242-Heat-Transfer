% PROJECT 2 ME 3242
%
% 
% Student:          Kevin Cayo
% Instructer:       Xinyu Zhao 
% Courese:          ME 3242
% Date:             Wed, May 5, 2022
%

clc
clear

%% PARAMETERS / KNOWN VARIABLES / EQUATIONS

% Variables
    T_inf = 1900;                           % [K]
    T_r = 1500;                             % [K], temp on the right boundary
    h = 250;                                % [W/m^2K]
    N = 18;                                 % total amount of subdomains; i-1 grid point
    E = 0.903085*5.67e-8*(1900^4);          % [W/m^2] radiation Energy

% Material A : Length, Thermal Conductivity, & Heat Source
    L_A = 0.001;                            % [m], length of A
    k_A = 0.8;                              % [W/mK], thermal conductivity for material A
    q_A = -100000 + E/L_A;                  % [W/m], heat source for material A

% Material B : Length, Thermal Conductivity, & Heat Source
    L_B = 0.0008;                           % [m], length of B
    k_B = 33.0;                             % [W/mK], thermal conductivity for material B
    q_B = 0.0;                              % [W/m], heat source for material B

% Equations
    L = L_A + L_B;                          % [m], total domain length
    dx = L/N;                               % size of the cells

% Setup variable arrays
    k = zeros(N,1);                         % array holds all the k values at each point
    qdot = zeros(N,1);                      % array holds all the heat source values at each point
    i_mid = floor(L_A/dx);                  % Index of the last point in block A

    for i = 1: i_mid
        k(i) = k_A;                         % conductivity values from material A
        qdot(i) = q_A;                      % heat source terms from material A
    end
    for i = (i_mid+1): N
        k(i) = k_B;                         % conductivity values from material B
        qdot(i) = q_B;                      % conductivity values from material B
    end

%% SET UP THE LEFT-HAND SIDE of AX=b

    A = zeros(N,N);                         % initialize a N x N dimension matrix with zero entries

%First set up values for the interior points
    for i = 2:N-1
        A(i,i-1) = k(i)+k(i-1);             % coefficient in front of point i-1
        A(i,i) = -(k(i-1)+2*k(i)+k(i+1));   % coefficient in front of point i
        A(i,i+1) = k(i)+k(i+1);             % coefficient in front of point i+1
    end 


% For the left boundary k values
    A(1,1) = -3*k(1) + +2*k(1)^2/(k(1)+h*dx/2);     % coefficient of point 1
    A(1,2) = 1*k(2);                                % coefficient of point 2

% For the right boundary k values
    A(N,N) = -3*k(N-1);                     % coefficient of point N
    A(N,N-1) = 1*k(N);                      % coefficient of point N-1

%% SET UP THE RIGHT-HAND SIDE of AX=b
    b = zeros(N,1);                         % initialize a N x 1 dimension matrix 

    for i = 2: N-1
        b(i) = -2*qdot(i)*(dx*dx);          % right-hand side values corresponding to point i
    end

    b(1) = -qdot(1)*(dx*dx) - 2*k(1)*h*dx/2*T_inf/(k(1) + h*dx/2);      % right-hand side values corresponding to point 1
    b(N) = -2*T_r*k(N) - qdot(N)*(dx*dx);                               % right-hand side values corresponding to point N

%% FIND THE SOLUTION OF T_I AT EACH GRID POINT
    X = A\b;
    T_l = (h*dx/2*T_inf + k(1)*X(1))/(k(1) + h*dx/2);

%% Visualize the temperature at each grid point as a function of x coordinates
    xgrid = 0.5*dx:dx:L-0.5*dx;
    xgrid = [0,xgrid,L];
    plot(xgrid,[T_l,X',T_r],'x-','Linewidth',1.25)
grid on
    
    xlabel('X (Meters)')
    ylabel('T (Kelvin)')
    title('Project 2 : N = 18')
    xlim([0 L])
    
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
              0, 0, L, 1
              L_A, 1, -L_A, -1
              k_A, 0, -k_B, 0];

%Set b vector
    bexact = [T_l, T_r, q_A*L_A^2/(2*k_A), q_A*L_A]';

%Calculate constants of integration
    C = Aexact\bexact;

%Temperature distribution in A
    T_A_exact = @(x) -q_A*x.^2/(2*k_A) + C(1)*x + C(2);

%Temperature distribution in B
    T_B_exact = @(x) C(3)*x + C(4);

%Plot temperature in A
    fplot(T_A_exact,[0, L_A])

%Plot Temperature in B
    fplot(T_B_exact,[L_A, L])

%Add a legend
    legend(sprintf('N=%i',N),'Exact, A','Exact, B')
    legend boxoff