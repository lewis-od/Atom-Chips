function [By, Bz] = eval_B(x,y,z, L, W, J)
%eval_B Calculate the magnetic field of a flat 2D wire

mu_0 = 4e-7 * pi;

By = eval_By(x,y,z, L/2, -W/2);
By = By + eval_By(x,y,z, -L/2, W/2);
By = By - eval_By(x,y,z, L/2, W/2);
By = By - eval_By(x,y,z, -L/2, -W/2);
By = ((mu_0*J)/(4*pi)) .* By;

Bz = eval_Bz(x,y,z, L/2, -W/2);
Bz = Bz + eval_Bz(x,y,z, -L/2, W/2);
Bz = Bz - eval_Bz(x,y,z, L/2, W/2);
Bz = Bz - eval_Bz(x,y,z, -L/2, -W/2);
Bz = ((mu_0*J)/(4*pi)) .* Bz;

end

