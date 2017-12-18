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
jx = sigma.*Ex;
jy = sigma.*Ey;

% Calculate f(kx, ky, z, d) in Fourier space
dims = size(jx);
Nx = dims(2);
Ny = dims(1);
kx = (0:1:(Nx-1)).*(sr/Nx);
ky = (0:1:(Ny-1)).*(sr/Ny);

[kx, ky] = meshgrid(kx, ky);
k = sqrt(kx.^2 + ky.^2);
f_hat = (mu_0/2) .* exp(-k.*z);

f_hat(isnan(f_hat)) = 0.0; % Remove infinities for inverse transform
f = real(ifft2(f_hat)); % Calculate f in real space

% Calculate B using the convolution theorem
Bx = conv2(jy, f, 'same');
By = -conv2(jx, f, 'same');

end

