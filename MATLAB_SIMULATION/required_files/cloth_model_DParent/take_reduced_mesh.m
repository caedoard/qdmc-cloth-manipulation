function [phi1,dphi1] = take_reduced_mesh(phi_big,dphi_big)     
% It extracts a 4x4 mesh from a 10x10 mesh.
%
% Author: David Parent, davidparentalonso@gmail.com
% Last review: 01/02/2021

n=4;

phi1 = zeros(n*n*3,1);
j=0;
for coord=[0 1 2]
    for i=[1 4 7 10 31 34 37 40 61 64 67 70 91 94 97 100]
        j=j+1;
        phi1(j,1) = phi_big(i+100*coord,1);
    end  
end

dphi1 = zeros(n*n*3,1);
j=0;
for coord=[0 1 2]
    for i=[1 4 7 10 31 34 37 40 61 64 67 70 91 94 97 100]
        j=j+1;
        dphi1(j,1) = dphi_big(i+100*coord,1);
    end  
end

end