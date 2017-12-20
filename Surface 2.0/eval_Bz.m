function [Bz] = eval_Bz(x, y, z, x0, y0)

Bz = log((x + x0) + ((x+x0).^2 + (y+y0).^2 + z.^2).^(1/2));

end


