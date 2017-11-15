function [Bx, By] = calc_field(params, phi)
%calc_field_from_potential Calculate the magnetic field from the given
%electric potential

mu_0 = pi*4e-7; % Permeability of free space
sigma = params(1);
d = params(2);
z = params(3);

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
kx = (0:1:length(jx));
ky = (0:1:length(jy));

[kx, ky] = meshgrid(kx, ky);
k = sqrt(kx.^2 + ky.^2);
f_hat = (mu_0/2)*d .* ((1-exp(-d*k))./(k*d)) .* exp(-k.*z);

f_hat(isnan(f_hat)) = 0.0; % Remove infinities for inverse transform
f = real(ifft2(f_hat)); % Calculate f in real space

% Calculate B using the convolution theorem
Bx = conv2(jy, f, 'same');
By = -conv2(jx, f, 'same');

end

