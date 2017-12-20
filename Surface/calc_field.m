function [Bx, By] = calc_field(phi, sigma, sr, d, z)
%calc_field_from_potential Calculate the magnetic field from the given
%electric potential

mu_0 = pi*4e-7; % Permeability of free space

% Calculate electric field
[Ex, Ey] = gradient(phi);
Ex = -1.*Ex;
Ey = -1.*Ey;

% Remove infinities
Ex(isnan(Ex)) = 0.0;
Ey(isnan(Ey)) = 0.0;

% Ohm's law
Jx = sigma.*Ex;
Jy = sigma.*Ey;

% Calculate f(kx, ky, z, d) in Fourier space
dims = size(Jx);
Nx = dims(2);
Ny = dims(1);
kx = (0:1:(Nx-1)).*(sr/Nx);
ky = (0:1:(Ny-1)).*(sr/Ny);

[kx, ky] = meshgrid(kx, ky);
k = sqrt(kx.^2 + ky.^2);
f_hat = (mu_0/2) .* exp(-k.*z); % Calculate f in inverse space

Jx_hat = fft2(Jx);
Jy_hat = fft2(Jy);

Jx_hat = fftshift(Jx_hat);
Jy_hat = fftshift(Jy_hat);

Bx_hat = Jy_hat .* f_hat;
By_hat = -Jx_hat .* f_hat;

Bx_hat(isnan(Bx_hat)) = 0.0; % Remove infinities for inverse transform
By_hat(isnan(By_hat)) = 0.0; % Remove infinities for inverse transform

Bx = abs(ifft2(Bx_hat));
By = abs(ifft2(By_hat));

% Smooth out noise
for i = 1:5
   Bx = smoothdata(Bx, 1);
   Bx = smoothdata(Bx, 2);
   By = smoothdata(By, 1);
   By = smoothdata(By, 2);
end

end

