function [phi1,dphi1] = take_reduced_mesh(phi_big,dphi_big, N,n)     
% Extract a n-by-n mesh from a N-by-N mesh.
%
% - Original Author: David Parent, davidparentalonso@gmail.com
% - Modified by: Adrià Luque, adria.alados@gmail.com
% - Last review: September 2021

if nargin < 4
    n = 4;
end

if nargin < 3
    N = 10;
end

i = 0:n^2-1;
redc = mod(i,n)*(N-1)/(n-1) + floor(i/n)*(N-1)/(n-1)*N + 1;

phi1 = zeros(n*n*3,1);
dphi1 = zeros(n*n*3,1);
j=0;
for coord=[0 1 2]
    for i=redc
        j=j+1;
        phi1(j,1) = phi_big(i + N^2*coord,1);
        dphi1(j,1) = dphi_big(i + N^2*coord,1);
    end  
end

end