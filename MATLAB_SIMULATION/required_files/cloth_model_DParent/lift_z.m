function out_ = lift_z(in_,param)
% It displaces the nodes in the vertical direction (z) a given distance.
%
% - Original Author: David Parent, davidparentalonso@gmail.com
% - Revised by: Adrià Luque, adria.alados@gmail.com
% - Last review: September 2021

i=0;
nx = param.row;
ny = param.col;

for r = nx:-1:1
    for c = 1:1:ny
        i=i+1;
        if r == 1
            d_up = 0;
        else
            d_up = param.z_sum;
        end
        out_{r,c} = {in_(i,1) in_(i,2) in_(i,3)+d_up};
    end
end
end