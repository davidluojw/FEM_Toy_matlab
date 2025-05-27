function generateFEMData(m, n)
    % 有限元数据文件生成器 - 生成平面应力条件下的四边形网格dat文件
    % 输入参数:
    %   m - x方向单元数量
    %   n - y方向单元数量
    
    % 确保输入为正整数
    if ~(isscalar(m) && isscalar(n) && m > 0 && n > 0 && floor(m) == m && floor(n) == n)
        error('m和n必须为正整数');
    end
    
    % 计算节点数和单元数
    numElements = m * n;              % 单元数量
    numNodesX = m + 1;                % x方向节点数
    numNodesY = n + 1;                % y方向节点数
    numNodes = numNodesX * numNodesY; % 总节点数
    
    % 材料属性
    youngsModulus = 3e7;              % 杨氏模量
    poissonsRatio = 0.3;              % 泊松比
    
    % 板的尺寸 (假设x方向从0到2，y方向从0到1)
    xMin = 0; xMax = 2;
    yMin = 0; yMax = 1;
    
    % 生成节点坐标
    nodes = generateNodeCoordinates(numNodesX, numNodesY, xMin, xMax, yMin, yMax);
    
    % 生成单元连接关系 (Q8单元 - 8节点四边形)
    elements = generateElementConnectivity(numNodesX, numNodesY);
    
    % 生成边界条件
    boundaryConditions = generateBoundaryConditions(numNodes);
    
    % 生成载荷条件 (假设在顶部边缘施加压力)
    traction = generateTractionConditions(numNodesX, numNodesY);
    
    % 文件名
    filename = sprintf('plate_mesh_%dx%d.dat', m, n);
    
    % 写入文件
    writeToFile(filename, numElements, numNodes, youngsModulus, poissonsRatio, ...
                nodes, elements, boundaryConditions, traction);
    
    fprintf('成功生成有限元数据文件: %s\n', filename);
end

function nodes = generateNodeCoordinates(numNodesX, numNodesY, xMin, xMax, yMin, yMax)
    % 生成节点坐标
    % 这里使用线性分布，实际应用中可以根据需要修改为其他分布方式
    x = linspace(xMin, xMax, numNodesX);
    y = linspace(yMin, yMax, numNodesY);
    
    % 预分配节点坐标矩阵
    nodes = zeros(numNodesX * numNodesY, 3);
    
    % 生成节点坐标 (考虑Q8单元的边中点)
    nodeId = 1;
    for j = 1:numNodesY
        for i = 1:numNodesX
            % 角点坐标
            nodes(nodeId, :) = [i-1, j-1, 0];  % 临时坐标，后续调整
            nodeId = nodeId + 1;
        end
    end
    
    % 调整为实际尺寸
    for i = 1:size(nodes, 1)
        % 计算在x和y方向的索引 (从1开始)
        col = mod(i-1, numNodesX) + 1;
        row = ceil(i / numNodesX);
        
        % 线性插值计算实际坐标
        nodes(i, 1) = xMin + (xMax - xMin) * (col - 1) / (numNodesX - 1);
        nodes(i, 2) = yMin + (yMax - yMin) * (row - 1) / (numNodesY - 1);
    end
    
    % 对于Q8单元，还需要生成边中点节点
    % 这里需要重新计算节点坐标，考虑边中点
    totalNodes = numNodesX * numNodesY + (numNodesX-1)*(numNodesY) + (numNodesY-1)*(numNodesX);
    nodes = zeros(totalNodes, 3);
    
    % 首先生成角点节点
    nodeId = 1;
    for j = 1:numNodesY
        for i = 1:numNodesX
            nodes(nodeId, 1) = xMin + (xMax - xMin) * (i-1) / (numNodesX - 1);
            nodes(nodeId, 2) = yMin + (yMax - yMin) * (j-1) / (numNodesY - 1);
            nodeId = nodeId + 1;
        end
    end
    
    % 生成水平边中点节点 (每一行的水平边)
    for j = 1:numNodesY-1
        for i = 1:numNodesX-1
            x1 = nodes((j-1)*numNodesX + i, 1);
            y1 = nodes((j-1)*numNodesX + i, 2);
            x2 = nodes((j-1)*numNodesX + i + 1, 1);
            y2 = nodes((j-1)*numNodesX + i + 1, 2);
            
            nodes(nodeId, 1) = (x1 + x2) / 2;
            nodes(nodeId, 2) = (y1 + y2) / 2;
            nodeId = nodeId + 1;
        end
    end
    
    % 生成垂直边中点节点 (每一列的垂直边)
    for j = 1:numNodesY-1
        for i = 1:numNodesX
            x1 = nodes((j-1)*numNodesX + i, 1);
            y1 = nodes((j-1)*numNodesX + i, 2);
            x2 = nodes(j*numNodesX + i, 1);
            y2 = nodes(j*numNodesX + i, 2);
            
            nodes(nodeId, 1) = (x1 + x2) / 2;
            nodes(nodeId, 2) = (y1 + y2) / 2;
            nodeId = nodeId + 1;
        end
    end
end

function elements = generateElementConnectivity(numNodesX, numNodesY)
    % 生成Q8单元的连接关系
    % Q8单元节点顺序: 4个角点 + 4个边中点 (按顺时针顺序)
    numElementsX = numNodesX - 1;
    numElementsY = numNodesY - 1;
    numElements = numElementsX * numElementsY;
    
    % 预分配单元矩阵 (每个单元8个节点)
    elements = zeros(numElements, 8);
    
    % 角点节点总数
    cornerNodes = numNodesX * numNodesY;
    % 水平边中点节点数
    horizontalMidNodes = (numNodesX - 1) * (numNodesY);
    % 垂直边中点节点数
    verticalMidNodes = (numNodesY - 1) * numNodesX;
    
    elementId = 1;
    for j = 1:numElementsY
        for i = 1:numElementsX
            % 计算当前单元的角点索引
            % 左下角, 右下角, 右上角, 左上角
            node1 = (j-1)*numNodesX + i;
            node2 = (j-1)*numNodesX + i + 1;
            node3 = j*numNodesX + i + 1;
            node4 = j*numNodesX + i;
            
            % 计算水平边中点索引 (左右边)
            node5 = cornerNodes + (j-1)*(numNodesX-1) + i; % 底边中点
            node7 = cornerNodes + (j-1)*(numNodesX-1) + i + (numNodesX-1)*numNodesY; % 顶边中点
            
            % 计算垂直边中点索引 (上下边)
            node6 = cornerNodes + horizontalMidNodes + (j-1)*numNodesX + i; % 右边中点
            node8 = cornerNodes + horizontalMidNodes + (j-1)*numNodesX + i + numNodesX*(numNodesY-1); % 左边中点
            
            % 组装单元节点
            elements(elementId, :) = [node1, node2, node3, node4, node5, node6, node7, node8];
            
            elementId = elementId + 1;
        end
    end
end

function boundaryConditions = generateBoundaryConditions(numNodes)
    % 生成边界条件
    % 这里假设左边界和下边界固定
    boundaryConditions = zeros(numNodes, 4);
    
    % 左边界节点 (x=0) 固定x和y方向位移
    for i = 1:size(boundaryConditions, 1)
        % 假设左边界节点编号为1, numNodesX+1, 2*numNodesX+1, ...
        if mod(i-1, numNodesX) == 0
            boundaryConditions(i, :) = [1, 0, 1, 0];
        % 下边界节点 (除左边界) 只固定y方向位移
        elseif i <= numNodesX
            boundaryConditions(i, :) = [0, 0, 1, 0];
        end
    end
end

function traction = generateTractionConditions(numNodesX, numNodesY)
    % 生成分布载荷条件
    % 这里假设在顶部边缘施加压力
    traction = cell(2, 1);
    traction{1} = 4; % 载荷边数量
    
    % 顶部边缘节点 (y=1)
    topEdgeNodes = zeros(1, numNodesX);
    for i = 1:numNodesX
        topEdgeNodes(i) = (numNodesY-1)*numNodesX + i;
    end
    
    % 分成4段施加压力
    traction{2} = [1, topEdgeNodes(1), topEdgeNodes(2), topEdgeNodes(3)];
    traction{3} = [2, topEdgeNodes(3), topEdgeNodes(4), topEdgeNodes(5)];
    traction{4} = [3, topEdgeNodes(5), topEdgeNodes(6), topEdgeNodes(7)];
    traction{5} = [4, topEdgeNodes(7), topEdgeNodes(8), topEdgeNodes(end)];
end

function writeToFile(filename, numElements, numNodes, youngsModulus, poissonsRatio, ...
                    nodes, elements, boundaryConditions, traction)
    % 将数据写入文件
    fileID = fopen(filename, 'w');
    if fileID == -1
        error('无法创建文件');
    end
    
    % 写入标题行
    dateStr = datestr(now, 'dd-mm-yy');
    fprintf(fileID, 'A Plate with Applied Pressure Loads                  %s\n', dateStr);
    fprintf(fileID, 'Plane stress\n');
    
    % 写入单元数和节点数
    fprintf(fileID, '      %d        %d                     ! No. of Elements, No. of Nodes\n', ...
            numElements, numNodes);
    
    % 写入材料属性
    fprintf(fileID, '      %.1e       %.1f 		       ! Young''s Modulus, Poisson''s Ratio\n', ...
            youngsModulus, poissonsRatio);
    
    % 写入节点部分标题
    fprintf(fileID, '# Nodes:\n');
    
    % 写入节点坐标
    for i = 1:size(nodes, 1)
        fprintf(fileID, '       %d   	  %.1e	  %.1e\n', i, nodes(i, 1), nodes(i, 2));
    end
    
    % 写入单元连接关系标题
    fprintf(fileID, '# Element Connectivity:\n');
    
    % 写入单元连接关系
    for i = 1:size(elements, 1)
        fprintf(fileID, '       %d         %d	  %d	 %d     %d    %d      %d     %d     %d\n', ...
                i, elements(i, 1), elements(i, 2), elements(i, 3), elements(i, 4), ...
                elements(i, 5), elements(i, 6), elements(i, 7), elements(i, 8));
    end
    
    % 写入边界条件标题
    fprintf(fileID, '# Boundary conditions:\n');
    
    % 写入边界条件
    for i = 1:size(boundaryConditions, 1)
        fprintf(fileID, '       %d   	  %d    %.1e	  %d    %.1e\n', ...
                i, boundaryConditions(i, 1), boundaryConditions(i, 2), ...
                boundaryConditions(i, 3), boundaryConditions(i, 4));
    end
    
    % 写入载荷条件标题
    fprintf(fileID, '# Traction:\n');
    
    % 写入载荷条件
    fprintf(fileID, '       %d      ! the number of edge which is applied traction\n', traction{1});
    for i = 2:size(traction, 1)
        fprintf(fileID, '       %d      %d     %d     %d     \n', traction{i}(1), traction{i}(2), ...
                traction{i}(3), traction{i}(4));
    end
    
    % 写入文件结束标记
    fprintf(fileID, '# End of file\n');
    
    % 关闭文件
    fclose(fileID);
end