function out_ = lift_z(in_,param)
% It displaces the nodes in the vertical direction (z) a given distance.
%
% Author: David Parent, davidparentalonso@gmail.com
% Last review: 01/02/2021

i=0;
import casadi.*
nx = param.row;
ny = param.col;

for r = nx:-1:1
    for c = 1:1:ny
        i=i+1;
        if r ==1
            puja = 0;
        else
            puja = param.z_sum;
        end
        out_{r,c} = [in_(i,1) in_(i,2) in_(i,3)+puja];
    end
end
end