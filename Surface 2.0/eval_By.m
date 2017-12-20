function [By] = eval_By(x, y, z, x0, y0)

arg = ((x+x0).*(y+y0))./(z.*((x+x0).^2 + (y+y0).^2 + z.^2).^(1/2));
By = atan(arg);

end

