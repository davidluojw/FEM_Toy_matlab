function generateQ4(m, n)
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
    numNodes = (m + 1) * (n + 1);
    fprintf(fid, '%8d%8d                     ! No. of Elements, No. of Nodes\n', numElements, numNodes);
    fprintf(fid, '%10.1e%8.1f               ! Young''s Modulus, Poisson''s Ratio\n', 3e7, 0.3);

    % node information
    fprintf(fid, '# Nodes:\n');
    for j = 0:n
        for i = 0:m
            nodeID = i + j * (m + 1) + 1;
            x = 2.0 * i / m;
            y = 1.0 * j / n +  (0.5 + 0.5 * j / n - 1.0 * j/ n ) / m * i;  
            fprintf(fid, '%8d%16.8e%16.8e\n', nodeID, x, y);
        end
    end

    % Element connectivity
    fprintf(fid, '# Element Connectivity:\n');
    for j = 0:n-1
        for i = 0:m-1
            elemID = i + j * m + 1;
            node1 = i + j * (m + 1) + 1;
            node2 = (i + 1) + j * (m + 1) + 1;
            node3 = (i + 1) + (j + 1) * (m + 1) + 1;
            node4 = i + (j + 1) * (m + 1) + 1;
            fprintf(fid, '%8d%8d%8d%8d%8d\n', elemID, node1, node2, node3, node4);
        end
    end

    % Boudanry conditions: arbitrarily generating,  we can modify it by hand
    fprintf(fid, '# Boundary conditions:\n');
    for j = 0:n
        for i = 0:m
            nodeID = i + j * (m + 1) + 1;
            % 左右边界(x=0和x=2)施加x方向约束
            if i == 0 || i == m
                bcx = 1;
                valx = 0.0;
            else
                bcx = 2;
                valx = 0.0;
            end
            
            % 上下边界(y=0和y=1)施加y方向约束
            if j == 0 || j == n
                bcy = 1;
                valy = 0.0;
            else
                bcy = 2;
                valy = 0.0;
            end
            
            fprintf(fid, '%8d%8d%12.3e%8d%12.3e\n', nodeID, bcx, valx, bcy, valy);
        end
    end

    % Traction
    fprintf(fid, '# Traction:\n');
    % 顶部边界施加载荷
    numEdges = m;
    fprintf(fid, '%8d      ! the number of edge which is applied traction\n', numEdges);
    for i = 0:m-1
        edgeID = i + 1;
        node1 = i + n * (m + 1) + 1;
        node2 = (i + 1) + n * (m + 1) + 1;
        fprintf(fid, '%8d%8d%8d\n', edgeID, node1, node2);
    end

    % End of file
    fprintf(fid, '# End of file\n');

    fclose(fid);
    fprintf('网格文件 "%s" 已生成，包含 %d 个单元和 %d 个节点。\n', filename, numElements, numNodes);
end