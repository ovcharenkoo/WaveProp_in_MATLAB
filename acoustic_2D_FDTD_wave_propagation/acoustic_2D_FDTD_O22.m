% Finite-differences in time domain(FDTD) acoustic wave propagation in 2D
% medium surrounded by simple sponge boundaries with exponential decay
% (Cerjan, 1985).
%
% We solve second order wave equation in time domain and displacement 
% formulation getting wavefield in terms of pressure field p.
%
% Acoustic medium is parametrized by acoustic velocity vp only.  We show 
% CFL number and number of points per wavelength prior running loop 
% over time steps.
%
% Conventional FD star-stencils provides accuracy O(2,2)
% Stencils: dx2: [1 -2 1]/dx2
% --------------------------------------------------------------
% The code is intentionally writen in a single file
% to simplify start up.
%
% The program does not save any files, add such option manually if needed.
% Drawing the wavefield is the most computationally demanding. Increase 
% IT_DISPLAY value to reduce output and accelerate computation.
%
% The goal is to provide a simple example of acoustic wave propagation
% in 2D isotropic medium.
%
% --------------------------------------------------------------
% Oleg Ovcharenko and Vladimir Kazei, 2018
%
% oleg.ovcharenko@kaust.edu.sa
% vladimir.kazei@kaust.edu.sa
%
% King Abdullah University of Science and Technology
% Thuwal, Saudi Arabia
% --------------------------------------------------------------

close all;
% Output every ... time steps
IT_DISPLAY = 10;

%% MODEL
% Model dimensions
nx = 401;
nz = 401;
dx = 10;    % [m]
dz = 10;    % [m]

% Elastic parameters
vp = 3300.0 * ones(nz, nx);         % velocity of compressional waves, [m/s]

%% TIME STEPPING
t_total = .55;                      % [sec] recording duration
dt = 0.7 * min(dx,dz)/max(vp(:));   % min grid space / max velocity
nt = round(t_total/dt);             % number of time steps
t = [0:nt]*dt;

CFL = max(vp(:)) * dt / min(dx,dz);
%% SOURCE
f0 = 10.0;                          % dominant frequency of the wavelet
t0 = 1.20 / f0;                     % half Ricker wavelet excitation time
factor = 1e10;                      % amplitude coefficient
angle_force = 90.0;                 % spatial orientation of source
                                    % (not relevant for acoustic case)
jsrc = round(nz/2);                 % source location along OZ
isrc = round(nx/2);                 % source location along OX

a = pi*pi*f0*f0;
dt2 = dt^2;
source_term = factor * exp(-a*(t-t0).^2);                            % Gaussian
% source_term =  -factor*2.0*a*(t-t0)*exp(-a*(t-t0)^2);                % First derivative of a Gaussian:
% source_term = factor * (1.0 - 2.0*a*(t-t0).^2).*exp(-a*(t-t0).^2);        % Ricker source time function (second derivative of a Gaussian):

force_x = sin(angle_force * pi / 180) * source_term * dt2 / (dx * dz);

min_wavelengh = 0.5*min(vp(vp>330))/f0;     % shortest wavelength bounded by velocity in the air

%% ABSORBING BOUNDARY (ABS)
abs_thick = min(floor(0.15*nx), floor(0.15*nz));    % thicknes of the layer
abs_rate = 0.3/abs_thick;                           % decay rate
lmargin = [abs_thick abs_thick];
rmargin = [abs_thick abs_thick];
weights = ones(nz+2,nx+2);
for iz = 1:nz+2
    for ix = 1:nx+2
        i = 0;
        j = 0;
        k = 0;
        if (ix < lmargin(1) + 1)
            i = lmargin(1) + 1 - ix;
        end
        if (iz < lmargin(2) + 1)
            k = lmargin(2) + 1 - iz;
        end
        if (nx - rmargin(1) < ix)
            i = ix - nx + rmargin(1);
        end
        if (nz - rmargin(2) < iz)
            k = iz - nz + rmargin(2);
        end
        if (i == 0 && j == 0 && k == 0)
            continue
        end
        rr = abs_rate * abs_rate * double(i*i + j*j + k*k );
        weights(iz,ix) = exp(-rr);
    end
end

%% SUMMARY
fprintf('#################################################\n');
fprintf('2D acoustic FDTD wave propagation in isotripic \nmedium in displacement formulation with \nCerjan(1985) boundary conditions\n');
fprintf('#################################################\n');
fprintf('Model:\n\t%d x %d\tgrid nz x nx\n\t%.1e x %.1e\t[m] dz x dx\n',nz, nx, dz,dx);
fprintf('\t%.1e x %.1e\t[m] model size\n',nx*dx, nz*dz);
fprintf('\t%.1e...%.1e\t[m/s] vp\n', min(vp(:)), max(vp(:)));
fprintf('Time:\n\t%.1e\t[sec] total\n\t%.1e\tdt\n\t%d\ttime steps\n',t_total,dt,nt);
fprintf('Source:\n\t%.1e\t[Hz] dominant frequency\n\t%.1f\t[sec] index time\n',f0,t0);
fprintf('Other:\n\t%.1f\tCFL number\n', CFL);
fprintf('\t%.2f\t[m] shortest wavelength\n\t%d, %d\t points-per-wavelength OX, OZ\n', min_wavelengh, floor(min_wavelengh/dx), floor(min_wavelengh/dz));
fprintf('#################################################\n');

%% ALLOCATE MEMORY FOR WAVEFIELD
p3 = zeros(nz+2,nx+2);            % Wavefields at t
p2 = zeros(nz+2,nx+2);            % Wavefields at t-1
p1 = zeros(nz+2,nx+2);            % Wavefields at t-2
% Coefficients for derivatives
co_dxx = 1/dx^2;
co_dzz = 1/dz^2;

%% Loop over TIME
tic;
for it = 1:nt
    p3 = zeros(size(p2));
    % Second-order derivatives
    dp_dxx = co_dxx * (p2(2:end-1,1:end-2) - 2*p2(2:end-1,2:end-1) + p2(2:end-1,3:end));
    dp_dzz = co_dzz * (p2(1:end-2,2:end-1) - 2*p2(2:end-1,2:end-1) + p2(3:end,2:end-1));
    % U(t) = 2*U(t-1) - U(t-2) + G dt2/rho;
    p3(2:end-1,2:end-1) = 2.0*p2(2:end-1,2:end-1) - p1(2:end-1,2:end-1) + (vp.^2).*(dp_dxx + dp_dzz).*dt2;
    % Add source term
    p3(jsrc, isrc) = p3(jsrc, isrc) + force_x(it); 
    % Exchange data between t-2 (1), t-1 (2) and t (3) and apply ABS
    p1 = p2 .* weights;
    p2 = p3 .* weights;
    % Output
    if mod(it,IT_DISPLAY) == 0
        fprintf('Time step: %d \t %.4f s\n',it, single(t(it)));
        imagesc(p3); colorbar; 
        axis equal tight; colormap jet;
        title(['Step = ',num2str(it),'/',num2str(nt),', Time: ',sprintf('%.4f',t(it)),' sec']);
        drawnow;
    end
end
toc; disp('End');
