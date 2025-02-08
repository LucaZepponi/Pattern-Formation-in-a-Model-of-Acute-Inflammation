clear
close all
clc

%% Model parameters
Dm = 0.45;   % Macrophage's diffusion coefficient
Dc = 1;      % Chemokines' diffusion coefficient
alpha = 0.5; % Chemokine receptor's saturation parameter
beta = 0.6;  % Anti-inflammatory cytokine's inhibition parameter
rho = 1;     % Exponent for inhibitory effect
m0 = 10;     % Initial macrophage density
a0 = (-1+sqrt(1+4*beta*m0))/(2*beta); % Initial nti-inﬂammatory cytokine density
tau = 80;    % Time scale of anti-inflammatory cytokines
chi = 2.8;   % Chemotaxis coefficient

%% Numerical parameters
L = 80;          % Length of the side of the spatial domain square
dx = 0.5;        % Space step 1
dy = 0.5;        % Space step 2
Nx = L / dx + 1; % Number of spatial points 1
Ny = L / dy + 1; % Number of spatial points 2
dt = 0.05;       % Time step
T = 9000;        % Total simulation time
Nt = T / dt + 1; % Number of time steps

%% Spatial and temporal grid
xi = 0:dx:L; % Spatial nodes 1
yk = 0:dy:L; % Spatial nodes 1
tj = 0:dt:T; % Time nodes

%% Variables' initialisation
m = m0 * ones(Nx, Ny);                 % Macrophage concentration
c = a0 + 0.1 * (2 * rand(Nx, Ny) - 1); % Chemokine concentration
a = a0 * ones(Nx, Ny);                 % Concentration of anti-inflammatory cytokines

%% Creating the VideoWriter object
videoFileName = 'evolution_macrophages.mp4'; % Nome del file video
writerObj = VideoWriter(videoFileName, 'MPEG-4');
writerObj.FrameRate = 10;                    % Set the framerate
open(writerObj);                             % Opens the video file for writing

%% Creazione della figura
figure;
h = imagesc(xi, yk, m'); % Initialise display
xlabel('Spatial position x');
ylabel('Spatial position y');
colorbar;
colormap(jet);
axis tight;
titleHandle = title(sprintf('Spatial pattern of m(x,y,t) for t = %.2f', 0));

%% Dynamics
for j = 1:Nt-1 % Time loop
    m_old = m;
    c_old = c;
    a_old = a;
    
    % Space loops to update m, c, and a
    for i = 2:Nx-1 % (boundary excluded)
        for k = 2:Ny-1 % (boundary excluded)
          f_ip = (chi * m_old(i+1, k)) / ((1 + alpha * c_old(i+1, k))^2);
            f_kp = (chi * m_old(i, k+1)) / ((1 + alpha * c_old(i, k+1))^2);
            f = (chi * m_old(i, k)) / ((1 + alpha * c_old(i, k))^2);
            f_im = (chi * m_old(i-1, k)) / ((1 + alpha * c_old(i-1, k))^2);
            f_km = (chi * m_old(i, k-1)) / ((1 + alpha * c_old(i, k-1))^2);
            F = (f_ip + f)/(2*dx^2) * (c_old(i+1,k) - c_old(i,k)) - (f + f_im)/(2*dx^2) * (c_old(i,k) - c_old(i-1,k)) + (f_kp + f)/(2*dx^2) * (c_old(i,k+1) - c_old(i,k)) - (f + f_km)/(2*dx^2) * (c_old(i,k) - c_old(i,k-1));
            g = m_old(i,k)/(1 + beta * (a_old(i,k))^rho );
            
            % Laplacians (m = m_xx + m_yy)
            m_xx_yy = (m_old(i+1, k) - 4*m_old(i, k) + m_old(i-1, k) + m_old(i, k+1) + m_old(i, k-1)) / dx^2; 
            c_xx_yy = (c_old(i+1, k) - 4*c_old(i, k) + c_old(i-1, k) + c_old(i,k+1) + c_old(i,k-1)) / dx^2;
            a_xx_yy = (a_old(i+1, k) - 4*a_old(i, k) + a_old(i-1, k) + a_old(i,k+1) + a_old(i,k-1)) / dx^2;
            
            % Update of m, c and a
            m(i, k) = m_old(i, k) + dt * (Dm * m_xx_yy - F);
            c(i, k) = c_old(i, k) + dt * (Dc * c_xx_yy - c_old(i, k) + g);
            a(i, k) = a_old(i, k) + dt * (Dc * a_xx_yy - a_old(i, k) +g) / tau; % Aggiornamento di a
        end
    end
    
    %% Condizioni al contorno periodiche
    m(1, :)   = m(end-1, :);
    m(end, :) = m(2, :); 
    m(:, 1)   = m(:, end-1); 
    m(:, end) = m(:, 2);

    c(1, :)   = c(end-1, :);
    c(end, :) = c(2, :); 
    c(:, 1)   = c(:, end-1); 
    c(:, end) = c(:, 2);
    
    a(1, :)   = a(end-1, :);
    a(end, :) = a(2, :); 
    a(:, 1)   = a(:, end-1); 
    a(:, end) = a(:, 2);
    
    %% Updating the visualisation
    if mod(j,10) == 0
        set(h, 'CData', m'); 
        set(titleHandle, 'String', sprintf('Spatial pattern of m(x,y,t) for t = %.2f', j * dt));
        
        drawnow;                      % Update figure
        frame = getframe(gcf);        % Capture current frame
        writeVideo(writerObj, frame); % Writes the frame in the video
        
        pause(0.05);                  % Pause for display
    end
end

%% Closing the video file
close(writerObj);
disp('Video successfully saved!')