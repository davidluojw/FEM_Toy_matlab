function generateQ8(m, n)
    % generateMeshFile
    % Input parameter:
    %   m - grid number of element along x-coor 
    %   n - grid number of element along y-coor 

    % outputfile name
    filename = sprintf('plate_mesh_%dx%d.dat', m, n);
    fid = fopen(filename, 'w');

    % mesh informations
    fprintf(fid, 'A Plate with Applied Pressure Loads                  %s\n', datestr(now, 'dd-mmm-yy'));
    fprintf(fid, 'Plane stress\n');
    numElements = m * n;
    numNodes1 = 2*m;
    numNodes2 = 2*n;
    numNodes3 = 2*m + 1;
    numNodes4 = m + 1;
    numNodes = (2*m + 1) * (n + 1) + (m + 1) * n;
    fprintf(fid, '%8d%8d                     ! No. of Elements, No. of Nodes\n', numElements, numNodes);
    fprintf(fid, '%10.1e%8.1f               ! Young''s Modulus, Poisson''s Ratio\n', 3e7, 0.3);

    % node information
    fprintf(fid, '# Nodes:\n');
    nodeID = 0;
    for j = 0: numNodes2
        if mod(j, 2) == 0
            for i = 0:numNodes1
                nodeID = nodeID + 1;
                x = 2.0 * i / (2*m);
                y = 1.0 * j / numNodes2 +  (0.5 + 0.5 * j / numNodes2 - 1.0 * j/ numNodes2 ) / numNodes1 * i;  
                fprintf(fid, '%8d%16.8e%16.8e\n', nodeID, x, y);
            end
        else
            for i = 0:m 
                nodeID = nodeID + 1;
                x = 2.0 * i / m;
                y = 1.0 * j / numNodes2 +  (0.5 + 0.5 * j / numNodes2 - 1.0 * j/ numNodes2 ) / m * i;  
                fprintf(fid, '%8d%16.8e%16.8e\n', nodeID, x, y);
            end
        end
    end

    % Element connectivity
    fprintf(fid, '# Element Connectivity:\n');
    for j = 0:n-1
        for i = 0:m-1
            elemID = i + j * m + 1;
            node1 = (numNodes3 + numNodes4) * j + 2 * i + 1;
            node2 = node1 + 2;
            node3 = node2 + (numNodes3 + numNodes4);
            node4 = node1 + (numNodes3 + numNodes4);
            node5 = node1 + 1;
            node6 = (numNodes3 + numNodes4) * j + numNodes3  + i + 2;
            node7 = node4 + 1;
            node8 = (numNodes3 + numNodes4) * j + numNodes3  + i + 1;
            fprintf(fid, '%8d%8d%8d%8d%8d%8d%8d%8d%8d\n', elemID, node1, node2, node3, node4, node5, node6, node7, node8);
        end
    end

    % Boudanry conditions: arbitrarily generating,  we can modify it by hand
    fprintf(fid, '# Boundary conditions:\n');
    nodeID = 0;
    for j = 0: numNodes2
        if mod(j, 2) == 0
            for i = 0:numNodes1
                nodeID = nodeID + 1;
                if i == 0 || i == numNodes1 || j == 0
                    bcx = 1;
                    valx = 0.0;
                    bcy = 1;
                    valy = 0.0;
                else
                    bcx = 2;
                    valx = 0.0;
                    bcy = 2;
                    valy = 0.0;
                end
                
                fprintf(fid, '%8d%8d%12.3e%8d%12.3e\n', nodeID, bcx, valx, bcy, valy);
            end
        else
            for i = 0:m 
                nodeID = nodeID + 1;
                if i == 0 || i == m || j == 0
                    bcx = 1;
                    valx = 0.0;
                    bcy = 1;
                    valy = 0.0;
                else
                    bcx = 2;
                    valx = 0.0;
                    bcy = 2;
                    valy = 0.0;
                end
                
                
                fprintf(fid, '%8d%8d%12.3e%8d%12.3e\n', nodeID, bcx, valx, bcy, valy);
            end
        end
    end

    % Traction
    fprintf(fid, '# Traction:\n');
    % 顶部边界施加载荷
    numEdges = m;
    fprintf(fid, '%8d      ! the number of edge which is applied traction\n', numEdges);
    for i = 0:m-1
        edgeID = i + 1;
        node1 = (numNodes3 + numNodes4) * n + 2 * i + 1;
        node2 = node1 + 1;
        node3 = node2 + 1;
        fprintf(fid, '%8d%8d%8d%8d\n', edgeID, node1, node2, node3);
    end

    % End of file
    fprintf(fid, '# End of file\n');

    fclose(fid);
    fprintf('网格文件 "%s" 已生成，包含 %d 个单元和 %d 个节点。\n', filename, numElements, numNodes);
end