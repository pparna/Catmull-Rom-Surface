clear;
clc;

numPoints = 10; % control point grid size defined by (numPoints)x(numPoints)

% compute the coordinates of the points
points_x = zeros(numPoints,numPoints);
points_y = zeros(numPoints,numPoints);
points_z = zeros(numPoints,numPoints);
for i = 1 : numPoints
   for j = 1 : numPoints
       points_x(i,j) = i;
       points_y(i,j) = j;
       points_z(i,j) = sin(i) * cos(j); % unsmooth surface
       %if(mod(i,2) == 0) % smooth surface
       %   points_z(i,j) = 1; 
       %end
   end
end
   
% simple function to construct a vector point from the separate arrays
ptIndexed = @(i,j) [points_x(i,j), points_y(i,j), points_z(i,j)];

% construct and render the patches
for j = 1 : (numPoints - 3)
    for i = 1 : (numPoints - 3)
        plotX = [];
        plotY = [];
        plotZ = [];
        plotNx = [];
        plotNy = [];
        plotNz = [];
        
        % load the control points for this patch
        p_ip3j = ptIndexed(i+3,j); p_ip3jp1 = ptIndexed(i+3,j+1); p_ip3jp2 = ptIndexed(i+3,j+2); p_ip3jp3 = ptIndexed(i+3,j+3);
        p_ip2j = ptIndexed(i+2,j); p_ip2jp1 = ptIndexed(i+2,j+1); p_ip2jp2 = ptIndexed(i+2,j+2); p_ip2jp3 = ptIndexed(i+2,j+3);
        p_ip1j = ptIndexed(i+1,j); p_ip1jp1 = ptIndexed(i+1,j+1); p_ip1jp2 = ptIndexed(i+1,j+2); p_ip1jp3 = ptIndexed(i+1,j+3);
        p_ij   = ptIndexed(i,j);   p_ijp1   = ptIndexed(i,j+1);   p_ijp2   = ptIndexed(i,j+2);   p_ijp3   = ptIndexed(i,j+3);
        
        for u = 0 : 0.05 : 1
           for w = 0 : 0.05 : 1
               % parabolic blending matrix
               H = [-0.5, 1.5, -1.5, 0.5;
                    1, -2.5, 2, -0.5;
                    -0.5, 0, 0.5, 0;
                    0, 1, 0, 0];
               
               % construct the separate x,y,z matrices
               P_x = [p_ip3j(1), p_ip3jp1(1), p_ip3jp2(1), p_ip3jp3(1);
                      p_ip2j(1), p_ip2jp1(1), p_ip2jp2(1), p_ip2jp3(1);
                      p_ip1j(1), p_ip1jp1(1), p_ip1jp2(1), p_ip1jp3(1);
                      p_ij(1)  , p_ijp1(1)  , p_ijp2(1)  , p_ijp3(1)];
                   
               P_y = [p_ip3j(2), p_ip3jp1(2), p_ip3jp2(2), p_ip3jp3(2);
                      p_ip2j(2), p_ip2jp1(2), p_ip2jp2(2), p_ip2jp3(2);
                      p_ip1j(2), p_ip1jp1(2), p_ip1jp2(2), p_ip1jp3(2);
                      p_ij(2)  , p_ijp1(2)  , p_ijp2(2)  , p_ijp3(2)];
                
               P_z = [p_ip3j(3), p_ip3jp1(3), p_ip3jp2(3), p_ip3jp3(3);
                      p_ip2j(3), p_ip2jp1(3), p_ip2jp2(3), p_ip2jp3(3);
                      p_ip1j(3), p_ip1jp1(3), p_ip1jp2(3), p_ip1jp3(3);
                      p_ij(3)  , p_ijp1(3)  , p_ijp2(3)  , p_ijp3(3)];
            
               % find the points on the CR surface
               p_x = [u^3, u^2, u, 1] * H * P_x * H' * [w^3; w^2; w; 1];
               p_y = [u^3, u^2, u, 1] * H * P_y * H' * [w^3; w^2; w; 1];
               p_z = [u^3, u^2, u, 1] * H * P_z * H' * [w^3; w^2; w; 1];
               
               % find the tangent in the u direction
               t_x = [3*u^2, 2*u, 1, 0] * H * P_x * H' * [w^3; w^2; w; 1];
               t_y = [3*u^2, 2*u, 1, 0] * H * P_y * H' * [w^3; w^2; w; 1];
               t_z = [3*u^2, 2*u, 1, 0] * H * P_z * H' * [w^3; w^2; w; 1];
               
               % find the tangent in the w direction
               b_x = [u^3, u^2, u, 1] * H * P_x * H' * [3*w^2; 2*w; 1; 0];
               b_y = [u^3, u^2, u, 1] * H * P_y * H' * [3*w^2; 2*w; 1; 0];
               b_z = [u^3, u^2, u, 1] * H * P_z * H' * [3*w^2; 2*w; 1; 0];
               
               % find the surface normal
               normal = cross([t_x, t_y, t_z], [b_x, b_y, b_z]);
               normal = normal / norm(normal);
               
               % store the computed data
               plotX = [plotX, p_x];
               plotY = [plotY, p_y];
               plotZ = [plotZ, p_z];
               
               plotNx = [plotNx, normal(1)];
               plotNy = [plotNy, normal(2)];
               plotNz = [plotNz, normal(3)];
           end
        end
        
        % reformat the data for rendering
        n = sqrt(size(plotX,2));
        plotX = reshape(plotX, n, n);
        plotY = reshape(plotY, n, n);
        plotZ = reshape(plotZ, n, n);
        vertNormals = reshape([plotNx plotNy plotNz], n, n, 3);
        
        % draw the surface patch
        h = surf(plotX, plotY, plotZ, 'VertexNormals', vertNormals, 'EdgeColor', 'none');
        h.FaceLighting = 'gouraud';
        h.AmbientStrength = 0;
        h.DiffuseStrength = 1;
        h.SpecularStrength = 0;
        hold on;
    end
end

% draw the control points
scatterX = reshape(points_x, 1, numPoints * numPoints);
scatterY = reshape(points_y, 1, numPoints * numPoints);
scatterZ = reshape(points_z, 1, numPoints * numPoints);
scatter3(scatterX, scatterY, scatterZ, 3, 'filled', 'MarkerFaceColor', 'red');

% set some general parameters for the plot
colormap white;
axis equal;
view(-26,59);
xlabel('x');
ylabel('y');
zlabel('z');
lightangle(-29,85);
set(gcf, 'color', 'white');
grid off;