function write_gmsh(fname, Tri, X)
    version = 2.1;
    file_type = 0; % Asci
    data_size = 8; % ==sizeof(double)
    element_type = 2; %Triangles
    tags = [0, 1, 0];
    N_points = size(X, 1);
    N_elems = size(Tri, 1);

    handle = fopen(fname, 'w');

    % Write header
    fprintf(handle, '$MeshFormat\n');
    fprintf(handle, '%.1f %d %d\n', version, file_type, data_size);
    fprintf(handle, '$EndMeshFormat\n');

    % Write nodes
    fprintf(handle, '$Nodes\n');
    fprintf(handle, '%d\n', size(X, 1));
    fprintf(handle, '%d %f %f %f\n', [(1:N_points)' X]');
    fprintf(handle, '$EndNodes\n');

    % Write elements
    fprintf(handle, '$Elements\n');
    fprintf(handle, '%d\n', N_elems);

    for i = 1:size(Tri, 1)
        fprintf(handle, '%d %d %d ', i, element_type, length(tags));
        for j = 1:length(tags)
            fprintf(handle, '%d ', tags(j));
        end
        fprintf(handle, '%d %d %d\n', Tri(i, :));
    end
    fprintf(handle, '$EndElements\n');

    % finalize
    fclose(handle);
    
    fprintf('Mesh written to %s\n', fname);
end

