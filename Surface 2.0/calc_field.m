function [Bx, By, Bz] = calc_field(x, y, z, Jx, Jy, dx, dy, xq, yq)
%calc_field Calculate magnetic field using finite element method

% Number of current elements in each dimension
Nx = (max(x(1,:)) - min(x(1,:)))/dx;
Nx = floor(Nx);
Ny = (max(y(:,1)) - min(y(:,1)))/dy;
Ny = floor(Ny);

% Number of data points per current element
dims = size(x);
dxN = floor(dims(2)/Nx);
dyN = floor(dims(1)/Ny);

J_direction = ones(size(x)).*-1;

res = length(xq);
Bx = zeros(res, res);
By = zeros(res, res);
Bz = zeros(res, res);
% Loop over all current elements of conductor
for nx = 1:Nx
    for ny = 1:Ny
        % Centre of current element
        cx = (nx-1)*dx + dx/2;
        cx = cx - 5e-6;;
        cy = (ny-1)*dy + dy/2;
        cy = cy - 5e-6;
    
        % Indices of current density matrix corresponding to this current
        % element
        xN_min = (nx-1)*dxN + 1;
        xN_max = nx*dxN;
        yN_min = (ny-1)*dyN + 1;
        yN_max = ny*dyN;
        
        % Treat the segment as having uniform current density, given by
        % the average of the density across the element
        Jn_x = Jx(yN_min:yN_max, xN_min:xN_max);
        Jn_y = Jy(yN_min:yN_max, xN_min:xN_max);
        
        Jn = sqrt(Jn_x.^2 + Jn_y.^2);
        Jn = mean(mean(Jn));
        
        if Jn == 0
            % No current => no magnetic field
            continue
        end
        
        Jnx_mean = mean(mean(Jn_x));
        Jny_mean = mean(mean(Jn_y));
        
        Jn_dir = abs(Jnx_mean) > abs(Jny_mean);
        J_direction(yN_min:yN_max, xN_min:xN_max) = Jn_dir;
        
        % Calculate the field due to the current segment
        dBx = 0;
        dBy = 0;
        dBz = 0;
        if Jn_dir
            % Current is in x direction - contributes to By
            [dBy, dBz] = eval_B(xq-cx, yq-cy, z, dy, dx, Jn);
        else
            % Current is in y direction - contributes to Bx
            [dBx, dBz] = eval_B(xq-cx, yq-cy, z, dy, dx, Jn);
        end
        
        Bx = Bx + dBx;
        By = By + dBy;
        Bz = Bz + dBz;
    end
end

end

