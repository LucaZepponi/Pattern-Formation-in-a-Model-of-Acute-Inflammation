clc
clear
close all

%% Model parameters
Dm = 0.45;   % Macrophage's diffusion coefficient
Dc = 1;      % Chemokines' diffusion coefficient
alpha = 0.5; % Chemokine receptor's saturation parameter
beta = 0.4;  % Anti-inflammatory cytokine's inhibition parameter
m0 = 10;     % Initial macrophage density
tau_vals = [300, 250, 100, 65, 50, 25];  % Time scale of anti-inflammatory cytokines
chi = 3;     % Chemotaxis coefficient

%% Numerical parameters
L = 25;         % Spatial domain length
xspan = [0, L]; % Solution range
k1 = [559.45,664,785.83]; % Values of parameter k1

%% Condizioni iniziali: [C(0), U(0)]
C0 = [0.00145, 0.12, 0.66]; % Initial values of C (scelto in base ai dati del paper)
U0 = 0;    % Boundary condition U(0) = 0

%% Function's color
colors = ['g', 'r', 'k']; % Plot colors

%% Calculation of curve equations
for i = 1:length(k1)
    %% Solving the ODEs system with ode45
    [x, Y] = ode45(@(x, y) system_eqs14(x, y, Dc, Dm, alpha, beta, chi, k1(i)), ...
                    xspan, [C0(i), U0]);

    % Homogeneous Neumann boundary conditions
    Y(1,2) = 0; % C_x(0) = 0
    Y(end,2) = 0; % C_x(L) = 0

    % Plot the results
    hold on;
    plot(x, Y(:, 1), 'LineWidth', 2, 'Color', colors(i));
    xlabel('x');
    ylabel('C(x)');
    title('Solution for C(x) with Neumann conditions');
    grid on;
end

%% Plot legend
legend('k1 = 559.45', 'k1 = 664', 'k1 = 785.83');

%% NON DEVE COMPARIRE NELLA PRESENTAZIONE!!!!!!!
  % Salvataggio della figura in formato PDF
  output_folder = 'C:\Users\UTENTE\OneDrive\Desktop\FigProjMAthBio\Fig5';
  filename = fullfile(output_folder, 'f5_B.png');
  saveas(gcf, filename); % Salva la figura corrente

%% ODEs' system definition
function dydx = system_eqs14(x, y, Dc, Dm, alpha, beta, chi, k1)
    % y(1) = C,
    % y(2) = C_x
    dydx = zeros(2,1);
    
    C = y(1);
    U = y(2);
    
    % Definition of the function m(C)
    mC = k1 * exp(- chi / (Dm * (alpha * (1 + alpha * C)))); % m(C)
    
    % ODEs
    dydx(1) = U;                              % dC/dx = U
    dydx(2) = (C - mC / (1 + beta * C)) / Dc; % dU/dx
    
    % Initial/boundary conditions
    if x == 0
        dydx(2) = 0; % C_x(0) = 0
    elseif x == 25
        dydx(2) = 0; % C_x(L) = 0
    end
end