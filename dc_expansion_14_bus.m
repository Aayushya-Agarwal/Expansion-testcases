function dc_expansion_14_bus
%   Detailed explanation goes here

    %mpc = loadcase('case14.m');
    mpc= loadcase('case14.m');
    
    
    %%%get the results from fmincon
    
    mpc.gencost = zeros(length(mpc.gen(:,1)), 7);
% 2 to show the form of gencost
mpc.gencost(:,1) = 2;
% how many coefficients
mpc.gencost(:,4) = 3;
% Make it flat start
mpc.bus(:,8) = 1.0;
mpc.bus(:,9) = 0.0;
mpc.gen(:,2) = 1.0;
%
id = 0;
cost_gen = 0;
gen_buses = [];
for idx = 1:length(mpc.gencost(:,6))
    if (mpc.gen(idx,8) == 1)
        %if mpc.gen(idx,9) == 0. && mpc.gen(idx, 10) == 0.)
            id = id + 1;
            % Quadratic coefficient
            mpc.gencost(idx, 5) = 0;
            % Linear coefficient
            mpc.gencost(idx, 6) = 10-1 - mod(id,10);
            cost_gen = cost_gen + mpc.gen(idx,2)*mpc.gencost(idx,6);
            gen_buses = [gen_buses; mpc.gen(idx, 1)];
            % Constant coefficient
            mpc.gencost(idx, 7) = 0;
        %end
    end
end
    cost_gen

% Set Vhi and Vlo
    mpc.bus(:,12) = 1.2;
    mpc.bus(:,13) = 0.8;
%
% Change if r is below a threshold

% Change if x is below a threshold

%
    pu_current = 60^0.5;
    current_base = 100*1e3./(sqrt(3)*mpc.bus(:,10));
    currentkA = (pu_current/1000).*current_base;
% current at 1pu MVA
    currentMVA = pu_current*100;
% Transformer index
    %xfmr_idx = find(mpc.branch(:,TAP) ~= 0 | mpc.branch(:,SHIFT) ~= 0);
% Branch index
    %branch_idx = find(mpc.branch(:,TAP) == 0 & mpc.branch(:,SHIFT) == 0);
% Ratings for branch
    %mpc.branch(branch_idx,6) = 9999; %currentMVA(1);
    %mpc.branch(branch_idx,7) = 9999; %currentMVA(1);   
    %mpc.branch(branch_idx,8) = 9999; %currentMVA(1);
% Ratings for transformers
    %mpc.branch(xfmr_idx, 6) = 9999;
    %mpc.branch(xfmr_idx, 7) = 9999;
    %mpc.branch(xfmr_idx, 8) = 9999;

%
% Set Pmin to 0
%mpc.gen(:,10) = -9999;
% Result_Store = [];
%load_factor = 2;
% mpc.bus(:,3:4) = mpc.bus(:,3:4)*load_factor;
% mpc.gen(:,2) = mpc.gen(:,2)*load_factor;
% In case you want to change initial conditions
%VM_val = 1.043;
%VA_val = 0;
%mpc.bus(:, VM) = VM_val;
% mpc.bus(2:end, VM) = VM_val;
%mpc.bus(:, VA) = VA_val;
% % mpc.bus(2789:end, VA) = VA_val;
%mpopt = mpoption('pf.alg', 'NR', 'pf.nr.max_it', 300, ...
%    'pf.tol', 1e-10,'pf.enforce_q_lims', 1);
%[RESULTS, SUCCESS, RAW] = MIPSOPF_SOLVER(mpc, MPOPTm);
%mpc2 = loadcase('case_ACTIVSg10k_OPF_created2.m');
%mfile = [name '2' '.m'];
%savecase(mfile, mpc);
    mpopt=mpoption;
%mpopt.opf.ac.solver = "KNITRO";
    mpopt.opf.ac.solver = "fmincon";
%mpopt.opf.current_balance = 1;
%mpopt.opf.v_cartesian = 1;
%mpopt.opf.flow_lim = 'I';
    mpopt.opf.ignore_angle_lim = 1;
    results = opf(mpc, mpopt);
    mpc = results;
    %%%
        
    
    [BBUS, BF, PBUSINJ, PFINJ] = makeBdc(mpc.baseMVA, mpc.bus, mpc.branch);
    B = full(BBUS);
    %%%%%%
    %statevector x = [theta1, ..., theta14, flow1,..., flow40, Pg1, Pg2,..., Pg4, gen1,...,
    num_branches = size(mpc.branch);
    num_bus = size(mpc.bus);
    num_gen = size(mpc.gen);
    num_branches = num_branches(1);
    num_bus = num_bus(1);
    num_gen = num_gen(1);

    %total_cols = num_bus + 2*num_branches + num_gen+num_branches+num_gen;
    %total_rows = num_bus + 2*num_branches + 2*num_branches + 2*num_gen;
    total_cols = num_bus + 2*num_branches +num_gen +num_gen + num_branches;
    %total_cols = num_bus + num_gen;
    total_rows = num_bus + 2*num_branches +  2*num_gen + 2*num_branches;
    %total_rows = num_bus + 2*num_branches ;
    %total_rows = num_bus;
    A = zeros(total_rows, total_cols);
    
    %Susceptance matrix for the first num_bus x num_bus submatrix
    A(1:num_bus, 1:num_bus) = B;
    rhs = zeros(total_rows,1);
    
    %obbjectuve function vector
    obj_c = zeros(total_cols,1);
    
    % define the flows, both from and to flows
    for branch = 1:num_branches
        branch_obj = mpc.branch(branch,:);
        from_bus =branch_obj(1);
        to_bus = branch_obj(2);
        B_pu = -1/branch_obj(4); 
        A(num_bus+branch*2-1, from_bus) = B_pu;
        A(num_bus+branch*2-1, to_bus) = -B_pu;
        A(num_bus+branch*2-1, num_bus+num_gen+branch*2-1) = -1;
        
        A(num_bus+branch*2, from_bus) = -B_pu;
        A(num_bus+branch*2, to_bus) = B_pu;
        A(num_bus+branch*2, num_bus+num_gen+branch*2) = -1;
        
        
    end
    
    %add the generator to the dc flows
    for gen =1:num_gen
        gen_obj = mpc.gen(gen,:);
        gen_bus = gen_obj(1);
        %A(gen_bus, num_bus+2*num_branches+gen) = -1;
        A(gen_bus, num_bus+gen) = -1;
    end
    
    % add loads to the rhs of the flows
    for bus =1:num_bus
        load_P = mpc.bus(bus,3)/100;
        rhs(bus) = rhs(bus) - load_P;
    end
    
    %add min and max generator
    for gen=1:num_gen
       Pmax = mpc.gen(gen, 9)/100;
       Pmin  =0;
       
       A(num_bus+2*num_branches+gen*2-1, num_bus+gen) = 1;
       A(num_bus+2*num_branches+gen*2-1, num_bus+2*num_branches+num_gen+gen) = -Pmax;

       %rhs(num_bus+2*num_branches+gen*2-1) = Pmax;
       
       A(num_bus+2*num_branches+gen*2, num_bus+gen) = -1;
       A(num_bus+2*num_branches+gen*2, num_bus+2*num_branches+num_gen+gen) = -Pmin;

       %rhs(num_bus+2*num_branches+gen*2) = Pmin;
        
    end
%         
    for gen =1:num_gen
       gen_cost = mpc.gencost(gen, 6);
       %obj_c(num_bus + 2*num_branches+gen) = gen_cost;
       obj_c(num_bus +gen) = gen_cost;
       
       obj_c(num_bus+num_gen+2*num_branches+gen) = 10;
        
        
    end
    
    for branch=1:num_branches
        
        obj_c(num_bus+2*num_branches+num_gen+num_gen+branch) =10;
        
    end
    
    
    %%define the upper bounds of each flow
    for branch=1:num_branches

       A(num_bus+num_branches*2+2*num_gen+branch*2-1, num_bus+num_gen+branch*2-1) = 1;
       A(num_bus+num_branches*2+2*num_gen+branch*2-1, num_bus+2*num_branches+num_gen+num_gen+branch) = -10;

       
       
       %rhs(num_bus+num_branches*2+2*num_gen+branch*2-1) = 10;
       
       A(num_bus+num_branches*2+2*num_gen+branch*2, num_bus+num_gen+branch*2) = 1;
       A(num_bus+num_branches*2+2*num_gen+branch*2, num_bus+2*num_branches+num_gen+num_gen+branch) = -10;
       %rhs(num_bus+num_branches*2+2*num_gen+branch*2) = 10;    
        
    end
    
    
    lower_bound = -Inf*ones(total_cols,1);
    lower_bound(1)=0;
    upper_bound = Inf*ones(total_cols,1);
    upper_bound(1)=0;
    lower_bound(15)=0;
    lower_bound(16)=0;
    lower_bound(17)=0;
    lower_bound(18)=0;
    lower_bound(19)=0;
    model.A = sparse(A);
    model.obj = obj_c;
    model.rhs = rhs;
    model.sense='======================================================<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<';
    %model.sense='=';
    model.vtype = 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCBBBBBBBBBBBBBBBBBBBBBBBBB';
    model.modelsense = 'min';
    model.lb= lower_bound;
    model.ub = upper_bound;
    
    gurobi_write(model, 'dcopf.lp');

    params.outputflag = 0;

    result = gurobi(model, params);

    disp(result);
    

end

