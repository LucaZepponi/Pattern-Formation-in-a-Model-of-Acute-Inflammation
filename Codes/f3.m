clear;
close all;
clc;

%% Model parameters
Dm = 0.45;   % Macrophage diffusion coefficient
Dc = 1;      % Chemokines' diffusion coefficient
alpha = 0.5; % Chemokine receptor saturation parameter
beta = 0.4;  % Anti-inflammatory cytokine inhibition parameter
rho = 1;     % Exponent for inhibitory effect
m0 = 10;     % Initial macrophage concentration
a0 = (-1+sqrt(1+4*beta*m0))/(2*beta); % Initial anti-inflammatory cytokines concentration
tau_vals = [30, 52, 55, 65, 75, 84, 95, 300];  % Time scale of anti-inflammatory cytokines
chi = 2;     % Chemotaxis coefficient

%% Numerical parameters
L = 200;         % Spatial domain length
dx = 0.5;        % Space step
Nx = L / dx + 1; % Number of spatial points
dt = 0.04;       % Time step
T = 11000;       % Total simulation time
Nt = T / dt + 1; % Number of time steps

%% Spatial and temporal grid
xi = 0:dx:L; % Space nodes
tj = 0:dt:T; % Time nodes

%% Initialisation of variables for each tau
m = zeros(Nx, Nt); % Density of macrophages
c = zeros(Nx, Nt); % Concentration of chemokines
a = zeros(Nx, Nt); % Concentration of anti-inflammatory cytokines
m(:, 1) = m0; % Initial uniform density of macrophages
c(:, 1) = a0 + 0.1 * (2 * rand(Nx, 1) - 1); % Initial concentration of chemokines with perturbation
a(:, 1) = a0; % Initial concentration of anti-inflammatory cytokines

%% Loop to update the time scale
for t = 1:length(tau_vals)
  tau = tau_vals(t); % Current value of tau
    
  % Update initial conditions: Use final fields from previous simulation
  if t ~= 1
    m(:, 1) = m(:, end);
    c(:, 1) = c(:, end);
    a(:, 1) = a(:, end);
  end

  % Time loop
  for j = 1:Nt-1
    % Laplacians
    m_xx = (m(3:end, j) - 2*m(2:end-1, j) + m(1:end-2, j)) / dx^2;
    c_xx = (c(3:end, j) - 2*c(2:end-1, j) + c(1:end-2, j)) / dx^2;
    a_xx = (a(3:end, j) - 2*a(2:end-1, j) + a(1:end-2, j)) / dx^2;

    % Compute chemotaxis term
    %  % Utilisation of vectorised operations to avoid grafted loops
    f = (chi * m(2:end-1, j)) ./ (1 + alpha * c(2:end-1, j)).^2;
    fp = (chi * m(3:end, j)) ./ (1 + alpha * c(3:end, j)).^2;
    fm = (chi * m(1:end-2, j)) ./ (1 + alpha * c(1:end-2, j)).^2;
    chemotaxis_term = ((fp + f) .* (c(3:end, j) - c(2:end-1, j)) - ...
                       (f + fm) .* (c(2:end-1, j) - c(1:end-2, j))) / (2 * dx^2);

    % Update m, c, and a
    m(2:end-1, j+1) = m(2:end-1, j) + dt * (Dm * m_xx - chemotaxis_term);
    c(2:end-1, j+1) = c(2:end-1, j) + dt * (Dc * c_xx - c(2:end-1, j) + ...
                      m(2:end-1, j) ./ (1 + beta * a(2:end-1, j).^rho));
    a(2:end-1, j+1) = a(2:end-1, j) + dt * (Dc * a_xx - a(2:end-1, j) + ...
                      m(2:end-1, j) ./ (1 + beta * a(2:end-1, j).^rho)) / tau;

    % Periodic boundary conditions
    m(1, j+1) = m(end-1, j+1);
    m(end, j+1) = m(2, j+1);
    c(1, j+1) = c(end-1, j+1);
    c(end, j+1) = c(2, j+1);
    a(1, j+1) = a(end-1, j+1);
    a(end, j+1) = a(2, j+1);
  end

  % Plot the results
  figure;
  imagesc(xi, tj, m'); % Transposed by m to have time on the y-axis
  xlabel('Spatial position x');
  ylabel('Time t');
  title(['Time evolution of m(x, t) for tau = ', num2str(tau_vals(t))]);
  colorbar;
  colormap(jet);
  axis tight;
end