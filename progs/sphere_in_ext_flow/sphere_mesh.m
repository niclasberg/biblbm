clc; clear; close all

% Radius and desired edge length
r = 8;
h_desired = 1;

% Compute number of points
area = 4*pi*r^2;
A_tri_desired = sqrt(3)/2 * h_desired^2
N_points = 2*round(area / A_tri_desired / 2)+2

fprintf('#points = %d, distance = %e, area = %e\n', N_points, h_desired, A_tri_desired);

% Generate mesh
[V, Tri] = ParticleSampleSphere('N', N_points, 'Etol', 1e-6);
V = r*V;


tr = TriRep(Tri, V(:, 1), V(:, 2), V(:, 3));
trimesh(tr), axis equal
write_gmsh(sprintf('sphere_r%d.gmsh', r), tr, V)
