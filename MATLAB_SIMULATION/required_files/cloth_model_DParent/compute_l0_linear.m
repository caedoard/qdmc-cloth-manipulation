function [mat_x,mat_y,mat_z] = compute_l0_linear(param,is_optimization)
% Compute initial spring length in each direction of the space x,y,z. It
% looks if there is a link for each type of connection and computes the
% initial length.
%
% - Original Author: David Parent, davidparentalonso@gmail.com
% - Modified by: Adrià Luque, adria.alados@gmail.com
% - Last review: September 2021

import casadi.*

i = 0;
nodeInitial = param.nodeInitial;
row = param.row;
col = param.col;

% is_optimization indicates if the output variables will be unknown
% If they are unknown, we use MX
mat_x = zeros(row*col, 6);
mat_y = zeros(row*col, 6);
if is_optimization
    mat_z = MX.zeros(row*col,6);
else
    mat_z = zeros(row*col, 6);
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
            mat_x(i,2) = nodeInitial{r,c}{1} - nodeInitial{nextRow,c}{1};
            mat_y(i,2) = nodeInitial{r,c}{2} - nodeInitial{nextRow,c}{2};
            mat_z(i,2) = nodeInitial{r,c}{3} - nodeInitial{nextRow,c}{3};
        end
        % Link 2
        if (c < col)
            mat_x(i,3) = nodeInitial{r,c}{1} - nodeInitial{r,nextCol}{1};
            mat_y(i,3) = nodeInitial{r,c}{2} - nodeInitial{r,nextCol}{2};
            mat_z(i,3) = nodeInitial{r,c}{3} - nodeInitial{r,nextCol}{3};
        end
        % Link 3
        if (r > 1)
            mat_x(i,5) = nodeInitial{r,c}{1} - nodeInitial{prevRow,c}{1};
            mat_y(i,5) = nodeInitial{r,c}{2} - nodeInitial{prevRow,c}{2};
            mat_z(i,5) = nodeInitial{r,c}{3} - nodeInitial{prevRow,c}{3};
        end
        % Link 4
        if (c > 1)
            mat_x(i,6) = nodeInitial{r,c}{1} - nodeInitial{r,prevCol}{1};
            mat_y(i,6) = nodeInitial{r,c}{2} - nodeInitial{r,prevCol}{2};
            mat_z(i,6) = nodeInitial{r,c}{3} - nodeInitial{r,prevCol}{3};
        end
    end
end
end