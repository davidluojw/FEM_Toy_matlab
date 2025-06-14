function model = create_2d_Q4_model(fem_data)

    model = struct();

    model.title = fem_data.title;

    model.nsd = 2;  % space dimension
    model.ndof = 2; % degree of freem of each node
    model.nen = 4;  % number of  nodes in each element
    model.ngp = 2;  % quadrature points for Gauss numberical integration
    model.deg_xi = 1;  % linear polynomial
    model.deg_eta = 1; % linear polynomial
    model.nnp = fem_data.num_nodes; % total number of nodes
    model.nel = fem_data.num_elements; %total numeber of elements
    model.neq = model.ndof * model.nnp;  % total number of equations, (x1,y1,x2,y2,...,xnnp, ynnp)
    model.nee = model.ndof * model.nen;  % number of dofs in each element
    model.bc = fem_data.bc;   % boundary condtion information of the model
    model.nodes = fem_data.nodes;  %the coordinates of nodes
    model.elements_connectivity = fem_data.elements_connectivity;  % element connectivity
    model.E = fem_data.E;   % Young's modulus
    model.nu = fem_data.nu; % poisson's ratio
    model.nint = model.ngp^(model.ndof); % quadrature points in 2D
    model.fact = 1;   % the amplification factor for plot
    model.strain_qdpt = zeros(3, model.nint, model.nel); % each element strain at quadrature point
    model.stress_qdpt = zeros(3, model.nint, model.nel); % each element stress at quadrature point
    model.strain_nodal = zeros(3, model.nnp); % each element strain at nodal point
    model.stress_nodal = zeros(3, model.nnp); % each element stress at nodal point
    model.counter_adjpt = zeros(1, model.nnp); % the number of adjent elements of each node

    model.f = zeros(model.neq, 1); % the global nodal forces
    model.d = zeros(model.neq, 1); % the solution dof
    model.r = zeros(model.neq, 1);    % the reaction nodal forces
    model.K = spalloc(model.neq, model.neq, (model.ndof * model.nen + 1)*model.neq);  % the global stiffness matrix
    model.M = model.K;   % the global mass matrix 
    model.d_dot = zeros(model.neq, 1);  % the first time derivative of solution
    model.d_ddot = zeros(model.neq, 1);  % the second time derivative of solution

    model.k_ele = zeros(model.nee, model.nee);   % dimension of element stiffness is 8 x 8
    model.f_ele = zeros(model.nee, 1);           % element force nodes
    model.m_ele = zeros(model.nee, model.nee);   % element mass matrix 8 x 8


    model.nd = 0;  % total dofs for essential boundary condition
    model.nbe = fem_data.nbe; % total number of edges for natural boundary condtion
    model.e_bc = zeros(model.neq, 1);  % the specified displacement for each dof for essential BC
    model.n_bc = fem_data.n_bc;  % first 2 rows: nodes belonging the edge which is applied natural boundary condition
                                 % next 4 rows: the nodal forces (fx1, fy1, fx2, fy2) applied at each dofs for the two nodes
    model.nbc_nodes = model.n_bc(1:2, :); 
    model.e_bc_dot = zeros(model.neq, 1);  % the specified velocity for each dof for essential BC in elastodynamics
    model.e_bc_ddot = zeros(model.neq, 1); % the specified acceleration for each dof for essential BC in elastodynamics
    model.P = zeros(model.neq, 1); % the external nodal forces applied to each dof
    model.b = zeros(model.nen * model.ndof, model.nel); % the nodal body forces in each element (bx1, by1, bx2, by2, bx3, by3)

    model.flags = zeros(model.neq, 1); % the dofs for essential boundary condition, denoted as 2

    for ii = 1:model.nnp
        if model.bc(ii, 1) == 1  % 1: Dirichlet BC
            index = (ii-1)*model.ndof + 1;
            model.e_bc(index) = model.bc(ii, 2);
            model.nd = model.nd + 1;
            model.flags(index) = 2;
        elseif model.bc(ii, 1) == 2  % 2: Concentrate Force
            index = (ii-1)*model.ndof + 1;
            model.P(index) = model.bc(ii, 2);
        end

        if model.bc(ii, 3) == 1   %  1: Dirichlet BC
            index = (ii-1)*model.ndof + 2;
            model.e_bc(index) = model.bc(ii, 4);
            model.nd = model.nd + 1;
            model.flags(index) = 2;
        elseif model.bc(ii, 3) == 2  % 2: Concentrate Force
            index = (ii-1)*model.ndof + 2;
            model.P(index) = model.bc(ii, 4);
        end
    end

    model.D = [1.0, model.nu, 0.0;
               model.nu, 1.0, 0.0;
               0.0, 0.0, (1.0-model.nu)/2.0] * model.E /(1.0 - model.nu^2);  % elasticity matrix

    model.G = model.E / (2.0 * (1.0 + model.nu));  % shear modulus

    % time solution option
    model.rho = 7850;   % the density of the material
    model.dt = 0.0001;    % time step 
    model.final_t = 1;  % final time 
    model.initial_t = 0;  % initial time
    model.nts = (model.final_t - model.initial_t)/ model.dt + 1;   % number of time steps
    model.time_intervals = 0 : model.dt : model.final_t;  % time sub-intervals
    model.d_ini = zeros(model.neq, 1); % initial condition, displacement
    model.v_ini = zeros(model.neq, 1); % initial condition，velocity
    model.disp = zeros(model.neq, model.nts);  % store all the solutions for all time stepss
    model.gamma = 5/6;  % integration parameter for Newmark method
    model.beta = 4/9;  % integration parameter for Newmark method
    model.rhoinf = 0.5;  % spectral radius
    model.alphaf = model.rhoinf / (model.rhoinf + 1);  % coefficient for Generalized alpha
    model.alpham = (2*model.rhoinf - 1)/ (model.rhoinf + 1); % coefficient for Generalized alpha
    model.beta = (1 - model.alpham + model.alphaf)^2 / 4; 
    model.gamma = 0.5 - model.alpham + model.alphaf; 


    model.IEN = model.elements_connectivity'; % IEN: mapping element node number to global node number
    model.LM = zeros(model.nen * model.ndof, model.nel); %LM: mapping element dof number to global node number
    model.ID = zeros(model.ndof, model.nnp); % ID: mapping the global node number to equation number, -1 for essential bc dofs
    
    


end
