
clc
clear
close all 

model = createpde(3)
[vertex, ~, ~] = stlread('Cube.stl');
vertex = vertex.Points';

importGeometry(model,'Cube.stl');


E = 3e6; % Modulus of elasticity in Pascals
nu = 0.3; % Poisson's ratio

m = 10000; % Material density in kg/m^3
c = elasticityC3D(E,nu);
a = 0;

  
specifyCoefficients(model,'m',m,'d',0,'c',c,'a',a,'f',[0;0;0]);
applyBoundaryCondition(model,'mixed','Face',1,'u',2,'EquationIndex',3);

mesh = generateMesh(model, "Hmax", 0.3);

node_surface = findNodes(model.Mesh,"nearest",vertex);

N_surface_node = size(node_surface,2);
N_surface_DOF = 3*N_surface_node;



FEM = assembleFEMatrices(model);
nodes = model.Mesh.Nodes;
N_all = size(nodes,2);
M = FEM.M; % 3 x nodes DOFs
K = FEM.K; % 3 x nodes DOFs


% Rearrange K and M matrices for N_master
idx_surface = [node_surface*3 - 2; node_surface*3 - 1; node_surface*3]; % Indices of N_master DOFs in the 3 x nodes DOFs
idx_surface = idx_surface(:)';

idx_all = 1:size(K,1); % All DOF indices
bool_surface = ismember(idx_all, idx_surface); % Logical indices of specific values in N_all
idx_reorder = [idx_all(bool_surface), idx_all(~bool_surface)];


K_reordered = K(idx_reorder, idx_reorder); % surface DOFs first, the others last.
M_reordered = M(idx_reorder, idx_reorder); 

M_inv = inv(M_reordered);

%% Original method

tic;
% Specify initial conditions
u0 = zeros(size(K_reordered, 1), 1); % Initial displacements (zero)
u0(1:3) = [0.1, 0, 0];
v0 = zeros(size(K_reordered, 1), 1); % Initial velocities (zero)

% Time parameters
tStart = 0; % Start time
tEnd = 0.1; % End time
numSteps = 10000; % Number of time steps
dt = (tEnd - tStart) / numSteps; % Time step size

% Solve the dynamic system using Newmark's method
u = zeros(size(u0, 1), numSteps+1); % Displacement matrix
v = zeros(size(v0, 1), numSteps+1); % Velocity matrix

u(:, 1) = u0; % Set initial displacements
v(:, 1) = v0; % Set initial velocities

M_inv_matmul_K = M_inv*K_reordered;

for i = 1:numSteps
    fprintf("%i\n", i)
    t = tStart + (i - 1) * dt; % Current time
    
    % Calculate acceleration at current time step
    a_current = -M_inv_matmul_K * u(:, i);
    
    % Update displacements and velocities using Newmark's method
    u(:, i+1) = u(:, i) + dt * v(:, i) + ((dt^2) / 2) * (a_current);
    v(:, i+1) = v(:, i) + dt * a_current;
end
u_original = u; % save trajectory

% Plot the displacement of a specific node over time
nodeIndex = 1; % Index of the node for which displacement is plotted
figure(1);
plot(tStart:dt:tEnd, u(nodeIndex, :));
xlabel('Time');
ylabel('Displacement');
title(['Displacement of Node ', num2str(nodeIndex)]);

elapsedTime_original = toc;
fprintf('Elapsed time is %.6f seconds\n', elapsedTime_original);


%%
outputVideo = VideoWriter('Original.avi');
open(outputVideo);

figure(2);

for i = 1:100:numSteps
fprintf("%d\n",i)
[~, inverse_reorder] = sort(idx_reorder);
pdeplot3D(mesh.Nodes + 5*reshape(u(inverse_reorder,i), size(nodes'))', mesh.Elements);
title('3D Mesh Visualization');
view(10,60);
frame = getframe(gcf);
writeVideo(outputVideo, frame);
end

close(outputVideo);

%% Preprocessing for Guyan reduction (without substructuring)

N_total = size(K_reordered,1);

% Matrix partitioned into: surface and internal DOFs.
% FYI, internal DOFs = substructure boundary DOFs + substructure body DOFs
% However, let's do CMS without substructuring as a toy example.
K_ss = K_reordered(1:N_surface_DOF, 1:N_surface_DOF); % Wikipedia notation: master-master 
K_si = K_reordered(1:N_surface_DOF, N_surface_DOF+1:N_total); % master-slave
K_is = K_reordered(N_surface_DOF+1:N_total, 1:N_surface_DOF); % slave-master
K_ii = K_reordered(N_surface_DOF+1:N_total, N_surface_DOF+1:N_total); % slave-slave

M_ss = M_reordered(1:N_surface_DOF, 1:N_surface_DOF);
M_si = M_reordered(1:N_surface_DOF, N_surface_DOF+1:N_total);
M_is = M_reordered(N_surface_DOF+1:N_total, 1:N_surface_DOF);
M_ii = M_reordered(N_surface_DOF+1:N_total, N_surface_DOF+1:N_total);



%% Guyan reduction

% Guyan reduction
% It's bad approximation, but let's disregard you are solving Kx=f problem
% There's no mass term the transformation matrix of Guyan reduction
% but let's just use it as a baseline model.

T_Guyan = [eye(size(K_ss));-K_ii\K_is];
K_Guyan = T_Guyan'*K_reordered*T_Guyan;
M_Guyan = T_Guyan'*M_reordered*T_Guyan;

K_Guyan = full(K_Guyan); % converted to Full, because it's condensed matrix! Using sparse matrix for this to compare with original method is unfair!
M_Guyan = full(M_Guyan);
M_Guyan_inv = full(inv(M_Guyan));



%%
tic;
% Specify initial conditions
u0 = zeros(size(K_Guyan, 1), 1); % Initial displacements (zero)
u0(1:3) = [0.1, 0, 0];
v0 = zeros(size(K_Guyan, 1), 1); % Initial velocities (zero)

% Time parameters
%tStart = 0; % Start time
%tEnd = 0.1; % End time
%numSteps = 1000; % Number of time steps
dt = (tEnd - tStart) / numSteps; % Time step size

% Solve the dynamic system using Newmark's method
u = zeros(size(u0, 1), numSteps+1); % Displacement matrix
v = zeros(size(v0, 1), numSteps+1); % Velocity matrix

u(:, 1) = u0; % Set initial displacements
v(:, 1) = v0; % Set initial velocities

M_inv_matmul_K = M_Guyan_inv*K_Guyan;

for i = 1:numSteps
    fprintf("%i\n", i)
    t = tStart + (i - 1) * dt; % Current time
    
    % Calculate acceleration at current time step
    a_current = -M_inv_matmul_K * u(:, i);
    
    % Update displacements and velocities using Newmark's method
    u(:, i+1) = u(:, i) + dt * v(:, i) + ((dt^2) / 2) * (a_current);
    v(:, i+1) = v(:, i) + dt * a_current;
end

% Plot the displacement of a specific node over time
nodeIndex = 1; % Index of the node for which displacement is plotted
figure(3);
u_Guyan = T_Guyan*u;
plot(tStart:dt:tEnd, u_Guyan(nodeIndex, :));
xlabel('Time');
ylabel('Displacement');
title(['Displacement of Node ', num2str(nodeIndex)]);

elapsedTime_Guyan = toc;
fprintf('Elapsed time is %.6f seconds\n', elapsedTime_Guyan);


%%

outputVideo = VideoWriter('Guyan_Reduced.avi');
open(outputVideo);

figure(4);

for i = 1:100:numSteps
fprintf("%d\n",i)
[~, inverse_reorder] = sort(idx_reorder);
pdeplot3D(mesh.Nodes + 5*reshape(u_Guyan(inverse_reorder,i), size(nodes'))', mesh.Elements);
title('3D Mesh Visualization');
view(10,60);
frame = getframe(gcf);
writeVideo(outputVideo, frame);
end

close(outputVideo);


%% CMS with Sub-Structuring (CMS-SS)
boundaryFacets = model.Mesh.findNodes('region', 'Face', 1);

% Find unique surface nodes
surfaceNodes = unique(boundaryFacets);

% Display the surface nodes
disp('Surface nodes indices:');
disp(idx_surface);

% Extract elements and nodes
elements = mesh.Elements;
nodes = mesh.Nodes;
numElements = size(elements, 2);
numNodes = size(nodes, 2);

% Create edge list for METIS input
edgeList = [];
for i = 1:numElements
    elementNodes = elements(:, i);
    edges = nchoosek(elementNodes, 2);
    edgeList = [edgeList; edges];
end

% Unique edges
uniqueEdges = unique(sort(edgeList, 2), 'rows');
numEdges = size(uniqueEdges, 1);

% Write METIS graph file
fileID = fopen('mesh.graph', 'w');
fprintf(fileID, '%d %d\n', numNodes, numEdges);

% Adjacency list creation
adjacencyList = cell(numNodes, 1);
for i = 1:size(uniqueEdges, 1)
    n1 = uniqueEdges(i, 1);
    n2 = uniqueEdges(i, 2);
    adjacencyList{n1} = [adjacencyList{n1}, n2];
    adjacencyList{n2} = [adjacencyList{n2}, n1];
end

for i = 1:numNodes
    fprintf(fileID, '%d ', adjacencyList{i});
    fprintf(fileID, '\n');
end

fclose(fileID);

% Display message
disp('METIS input file "mesh.graph" created successfully.');


%% Execute metis. (If you don't have matis installed, you can use the partition file I created)
path2metis = "/opt/homebrew/bin/gpmetis";
N_substructure = 4;

partitionFile = 'mesh.graph.part.' + string(N_substructure);

% Partition file이 존재하는지 확인
if ~exist(partitionFile, 'file')
    % Partition file이 없으면 METIS 실행
    system(path2metis + " mesh.graph " + string(N_substructure));
end

% Read partition file
fileID = fopen(partitionFile, 'r');
partitions = fscanf(fileID, '%d');
fclose(fileID);



%%
% new indexing array
Node_membership = [];
Partition_Nodes = [];
for part = 0:max(partitions)
    partitionElements = find(partitions == part);
    Node_membership = [Node_membership, ones(size(partitionElements')) * (part+1)];
    Partition_Nodes = [Partition_Nodes, partitionElements'];
end


idx_membership = [Node_membership'*3-2, Node_membership'*3-1, Node_membership']
idx_membership = idx_membership(:)';

Partition_idx = [Partition_Nodes*3-2; Partition_Nodes*3-1; Partition_Nodes*3];
Partition_idx = Partition_idx(:)';



% Extract boundary nodes
boundaryNodes = [];
for i = 1:numNodes
    neighbors = adjacencyList{i};
    for neighbor = neighbors
        if partitions(i) ~= partitions(neighbor)
            boundaryNodes = [boundaryNodes, i];
            break; % Once you have a neighbor having a different membership, you are a boundary node.
        end
    end
end

% Remove duplication
boundaryNodes = unique(boundaryNodes);
boundary_idx = [boundaryNodes*3-2; boundaryNodes*3-1; boundaryNodes*3];
boundary_idx = boundary_idx(:)';


%%
IsInBoundary = ismember(Partition_idx, boundary_idx);
idx_membership = idx_membership(~IsInBoundary);
Partition_idx = Partition_idx(~IsInBoundary);

IsOnSurface = ismember(Partition_idx, boundary_idx);
idx_membership = idx_membership(~IsOnSurface);
Partition_idx = Partition_idx(~IsOnSurface);


% Matrix indexing array
idx_reorder = [setdiff(idx_surface, Partition_idx), setdiff(boundary_idx, idx_surface), Partition_idx];

% If you want to make [1,2,3] located at the beginning, de-annotate the followings.
%dummy = ismember(idx_reorder, [1,2,3]);
%idx_reorder = [1,2,3, idx_reorder(~dummy)];

%%
% MATLAB's matrix is ordered with all nodes' x-coordinates first, 
% followed by all nodes' y-coordinates, and then all nodes' z-coordinates:
% x1, x2, x3, ..., xn, y1, y2, y3, ..., yn, z1, z2, z3, ..., zn.

% Before substructuring, this might be better if aranged into x1, y1, x1, x2, y2, z2, ...

idx_reorder_xyzbase = [1:3:size(K,1); 2:3:size(K,1); 3:3:size(K,1)]';
idx_reorder_xyzbase = idx_reorder_xyzbase(:)';
[~, idx_reorder_node_base] = sort(idx_reorder_xyzbase);

K_reordered = K(idx_reorder_node_base, idx_reorder_node_base);
K_reordered = K_reordered(idx_reorder,idx_reorder);

M_reordered = M(idx_reorder_node_base, idx_reorder_node_base);
M_reordered = M_reordered(idx_reorder,idx_reorder);


% Visualize
figure(9);
spy(M_reordered);
figure(10);
spy(K_reordered);

%%
N_total = size(K_reordered,2);
N_Surf_N_Boundary = N_total-size(Partition_idx,2);


% Matrix partitioned into: surface and internal DOFs.
% FYI, internal DOFs = substructure boundary DOFs + substructure body DOFs
% However, let's do CMS without substructuring as a toy example.
K_ss = K_reordered(1:N_Surf_N_Boundary, 1:N_Surf_N_Boundary); % Wikipedia notation: master-master 
K_si = K_reordered(1:N_Surf_N_Boundary, N_Surf_N_Boundary+1:N_total); % master-slave
K_is = K_reordered(N_Surf_N_Boundary+1:N_total, 1:N_Surf_N_Boundary); % slave-master
K_ii = K_reordered(N_Surf_N_Boundary+1:N_total, N_Surf_N_Boundary+1:N_total); % slave-slave

M_ss = M_reordered(1:N_Surf_N_Boundary, 1:N_Surf_N_Boundary);
M_si = M_reordered(1:N_Surf_N_Boundary, N_Surf_N_Boundary+1:N_total);
M_is = M_reordered(N_Surf_N_Boundary+1:N_total, 1:N_Surf_N_Boundary);
M_ii = M_reordered(N_Surf_N_Boundary+1:N_total, N_Surf_N_Boundary+1:N_total);


[Phi, Lambda] = eig(full(M_ii), full(K_ii), 'vector');

Nd = 10; % note that this is a hyper parameter

[~, idx] = sort(Lambda);
idx = idx(1:Nd);
Lambda = Lambda(idx);
Phi = Phi(:, idx);

T1 = eye(N_Surf_N_Boundary); % identity matrix for boundary DOFs
T2 = -K_ii \ K_is; % term that links boundary to interior nodes

% Concatenate T1, T2 and Phi to form the full transformation matrix
T_CMS = blkdiag(T1, Phi);
T_CMS(N_Surf_N_Boundary+1:end, 1:N_Surf_N_Boundary) = T2;

K_CMS = T_CMS'*K_reordered*T_CMS;
M_CMS = T_CMS'*M_reordered*T_CMS;
M_CMS_inv = inv(M_CMS);


%%
tic;
% Specify initial conditions
u0 = zeros(size(K_CMS, 1), 1); % Initial displacements (zero)
u0(1:3) = [0.1, 0, 0];
v0 = zeros(size(K_CMS, 1), 1); % Initial velocities (zero)

% Time parameters
%tStart = 0; % Start time
%tEnd = 0.1; % End time
%numSteps = 10000; % Number of time steps
%dt = (tEnd - tStart) / numSteps; % Time step size

% Solve the dynamic system using Newmark's method
u = zeros(size(u0, 1), numSteps+1); % Displacement matrix
v = zeros(size(v0, 1), numSteps+1); % Velocity matrix

u(:, 1) = u0; % Set initial displacements
v(:, 1) = v0; % Set initial velocities

M_inv_matmul_K = M_CMS_inv*K_CMS;

for i = 1:numSteps
    fprintf("%i\n", i)
    t = tStart + (i - 1) * dt; % Current time
    
    % Calculate acceleration at current time step
    a_current = -M_inv_matmul_K * u(:, i);
    
    % Update displacements and velocities using Newmark's method
    u(:, i+1) = u(:, i) + dt * v(:, i) + ((dt^2) / 2) * (a_current);
    v(:, i+1) = v(:, i) + dt * a_current;
end

% Plot the displacement of a specific node over time
nodeIndex = 1; % Index of the node for which displacement is plotted
figure(7);

u_CMS_SS = T_CMS*u;

plot(tStart:dt:tEnd, u_CMS_SS(nodeIndex, :));
xlabel('Time');
ylabel('Displacement');
title(['Displacement of Node ', num2str(nodeIndex)]);

f0_reduced = K_CMS*u0;
elapsedTime_CMS_SS = toc;
fprintf('Elapsed time is %.6f seconds\n', elapsedTime_CMS_SS);

%%
outputVideo = VideoWriter('CMS_reduced.avi');
open(outputVideo);

figure(6);

for i = 1:100:numSteps
fprintf("%d\n",i)

idx = 1:1:size(K,2);
idx = idx(idx_reorder_node_base);
idx = idx(idx_reorder);
[~, inverse_reorder] = sort(idx);

pdeplot3D(mesh.Nodes + 5*reshape(u_CMS_SS(inverse_reorder,i), size(nodes'))', mesh.Elements);
title('3D Mesh Visualization');
view(10,60);
frame = getframe(gcf);
writeVideo(outputVideo, frame);
end

close(outputVideo);


%%
fprintf('\n\nPerformance assessment\n')
fprintf('Elapsed time is %.6f seconds\n', elapsedTime_original);
fprintf('Elapsed time is %.6f seconds\n', elapsedTime_Guyan);
fprintf('Elapsed time is %.6f seconds\n', elapsedTime_CMS_SS);

%%
traj_vector_original = (u_original(1:1,1:end));
traj_vector_original = traj_vector_original(:);
traj_vector_original = traj_vector_original/norm(traj_vector_original);

traj_vector_Guyan = (u_Guyan(1:1,1:end));
traj_vector_Guyan=traj_vector_Guyan(:);
traj_vector_Guyan = traj_vector_Guyan/norm(traj_vector_Guyan);


traj_vector_CMS = (u_CMS_SS(1:1,1:end));
traj_vector_CMS=traj_vector_CMS(:);
traj_vector_CMS = traj_vector_CMS/norm(traj_vector_CMS);

fprintf("Guyan trajectory vector cosine theta with the ground truth : %f \n", dot(traj_vector_original,traj_vector_Guyan)  )
fprintf("CMS trajectory vector cosine thete with the ground truth : %f \n", dot(traj_vector_original,traj_vector_CMS)  )

%%
