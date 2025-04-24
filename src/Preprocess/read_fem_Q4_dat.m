function fem_data = read_fem_Q4_dat(filename)
    fid = fopen(filename, 'r');

    if fid == -1
        error('Cannot open the file: %s', filename);
    end

    % initialize the FEM data
    fem_data = struct();

    % read the title
    line = fgetl(fid);
    fem_data.title = sscanf(line, '%[^...]');

    fem_data.nbe = 0;
    fem_data.n_bc = zeros(6,1);
    
    % read the model parameters
    while ~feof(fid)  % whether read the end line of the file
        line = fgetl(fid);
        if startsWith(line, "Plane stress")
            fem_data.plane_stress = 1;
            fem_data.plane_strain = 0;
            line = fgetl(fid);
            data = sscanf(line, '%f%f');
            fem_data.num_elements = data(1);
            fem_data.num_nodes = data(2);
            line = fgetl(fid);
            data = sscanf(line, '%f%f');
            fem_data.E = data(1);
            fem_data.nu = data(2);
        elseif startsWith(line, "Plane strain")
            fem_data.plane_stress = 0;
            fem_data.plane_strain = 1;
            line = fgetl(fid);
            data = sscanf(line, '%f%f');
            fem_data.num_elements = data(1);
            fem_data.num_nodes = data(2);
            line = fgetl(fid);
            data = sscanf(line, '%f%f');
            fem_data.E = data(1) / (1.0 - data(2)^2);
            fem_data.nu = data(2) / (1.0 - data(2));
        elseif startsWith(line, "# Nodes:")
            for i = 1:fem_data.num_nodes
                line  = fgetl(fid);
                node_data = sscanf(line, '%f%f%f');
                fem_data.nodes(i, :) = node_data(2:end);
            end
        elseif startsWith(line, "# Element Connectivity:")
            for i = 1:fem_data.num_elements
                line = fgetl(fid);
                ele_data  = sscanf(line, '%f%f%f%f%f');
                fem_data.elements_connectivity(i, :) = ele_data(2:end);
            end
        elseif startsWith(line, "# Boundary conditions:")
            for i = 1:fem_data.num_nodes
                line = fgetl(fid);
                bc_data = sscanf(line, '%f%f%f%f%f');
                fem_data.bc(i, :) = bc_data(2:end);
            end
        elseif startsWith(line, "# Traction:")
            line = fgetl(fid);
            data = sscanf(line, '%f');
            fem_data.nbe = data;
            for i = 1:fem_data.nbe
                line = fgetl(fid);
                nbc_data = sscanf(line, "%f%f%f%f%f%f%f");
                fem_data.n_bc(:, i) = nbc_data(2:end);
            end
        end

        if startsWith(line, "# End of file")
            break;
        end
    end
    fclose(fid);
end