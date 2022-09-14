function [mat_x,mat_y,mat_z] = compute_l0_linear(param,is_optimization)
% Compute initial spring length in each direction of the space x,y,z. It
% looks if there is a link for each type of connection and computes the
% initial lenght.
%
% Author: David Parent, davidparentalonso@gmail.com
% Last review: 07/02/2021

import casadi.*

i = 0;
nodeInitial = param.nodeInitial;
row = param.row;
col = param.col;

% is_optimization indicates if the output variables will be unknown or not.
% If they are unknouw, we use MX to use
if is_optimization
    mat_x = MX.zeros(16,6);
    mat_y = MX.zeros(16,6);
    mat_z = MX.zeros(16,6);
end

for r=row:-1:1
    nextRow = r + 1;
    prevRow = r - 1;
    for c=1:col
        nextCol = c + 1;
        prevCol = c - 1;
        
        i= i+1;
        
        % Link 1
        if (r < row)
            l0 = nodeInitial{r, c} - nodeInitial{nextRow, c};
            mat_x(i,2) = l0(1);
            mat_y(i,2) = l0(2);
            mat_z(i,2) = l0(3);
        end
        % Link 2
        if (c < col)
            l0 = nodeInitial{r, c} - nodeInitial{r, nextCol};
            mat_x(i,3) = l0(1);
            mat_y(i,3) = l0(2);
            mat_z(i,3) = l0(3);
        end
        % Link 3
        if (r > 1)
            l0 = nodeInitial{r, c} - nodeInitial{prevRow, c};
            mat_x(i,5) = l0(1);
            mat_y(i,5) = l0(2);
            mat_z(i,5) = l0(3);
        end
        % Link 4
        if (c > 1)
            l0 = nodeInitial{r, c} - nodeInitial{r, prevCol};
            mat_x(i,6) = l0(1);
            mat_y(i,6) = l0(2);
            mat_z(i,6) = l0(3);
        end
    end
end
end