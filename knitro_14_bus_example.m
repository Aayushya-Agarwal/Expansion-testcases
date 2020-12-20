function knitro_14_bus_example
%UNTITLED Summary of this function goes here
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
            mpc.gencost(idx, 6) = 1 + mod(id,10);
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
    
    
    
    setGlobalx(mpc);
    
    options = knitro_options('algorithm',3, 'outlev',4, ...
                         'maxit',1000,'xtol',1e-4, ...
                         'feastol', 1e-4, 'opttol', 1e-4,'derivcheck',1);
    
    %%%%STATE VECTOR X=[Vr, Vi, P, Q]
    %Vr=X[bus_number*2-1]
    %Vi=X[bus_numer*2]
    %P=X[14*2+gen_number*2]
    %Q=X[14*2+gen_number*2-1]
    
    %number of constraints = 2*number of buses
    
    %objective = sum of generator cost
    gen_size=size(mpc.gen);
    bus_size=size(mpc.bus);
    sol_size=bus_size(1)*2+gen_size(1)*2;
    %x0 = ones(sol_size,1);
    %x0(2:2:sol_size) = 0;
    lb= -2*ones(sol_size,1);
    ub = 2*ones(sol_size,1);   
    
    for bus =1:bus_size(1)
       
        x0(bus*2-1)=mpc.bus(bus,8)*cos(mpc.bus(bus,9)*pi/180);
        x0(bus*2)=mpc.bus(bus,8)*sin(mpc.bus(bus,9)*pi/180);
        

        
    end
    
    for gen=1:gen_size(1)
       
        x0(bus_size(1)*2 +gen*2-1) = mpc.gen(gen,2)/100;
        x0(bus_size(1)*2 +gen*2) = -mpc.gen(gen,3)/100;
        
        ub(bus_size(1)*2 +gen*2-1) = mpc.gen(gen,9)/100;
        lb(bus_size(1)*2 +gen*2-1) = mpc.gen(gen,10)/100;

        lb(bus_size(1)*2 +gen*2) = -10;%mpc.gen(gen,9)/100;
        ub(bus_size(1)*2 +gen*2) = 10;%mpc.gen(gen,10)/100;        

        
    end
    x_fmincon = x0;
    A=[];
    b=[];

    Aeq=[];
    beq=[];
    
    xType = zeros(sol_size,1);
    x0 = ones(sol_size,1);
    x0(2:2:sol_size) = 0;
    % Call Knitro to solve the optimization model.
    %[sol,fval,exitflag,output,lambda, grad, hess] = ...
     %knitro_nlp(@objfun,x0,A,b,Aeq,beq,lb,ub,@constfun,[],options);
     
    %sol = fmincon(@objfun, x0, A, b, Aeq, beq, lb, ub, @constfun);
    
    [sol,fval,exitflag,output,lambda] = ...
     knitro_minlp(@objfun,x0,xType,A,b,Aeq,beq,lb,ub,@constfun,[],options,'mipoptions.opt');
     
     
    sol
end


function setGlobalx(val)
    global x
    x = val;
end

function r = getGlobalx
    global x
    r = x;
end
function [f] = objfun(X)
    disp '';
    mpc=  getGlobalx;
    f=0;
    bus_size = size(mpc.bus);
    gen_size=size(mpc.gen);
    for gen=1:gen_size(1)
        Pg = X(bus_size(1)*2+gen*2-1);
        f = f + Pg*(1+ mod(gen+1,10));
        
        g(2*bus_size(1)+gen*2-1) =  (1+ mod(gen+1,10));
    end
    
    
end


%function [c,ceq,Gc,Geq]= constfun(V)
function [c,ceq]= constfun(V)

    mpc=  getGlobalx;
    bus_size = size(mpc.bus);
    
    
    Y = zeros(bus_size(1), bus_size(1));

    c=zeros(2*bus_size(1),1);
    %c=[];
    %Gc=[];
    gen_size = size(mpc.gen);
    ceq=zeros(bus_size(1)*2,1);
    Geq =zeros(bus_size(1)*2 + gen_size(1)*2,bus_size(1)*2);

    
    %Bus voltage bounds
    for bus=1:bus_size(1)
       
        Vr_node = bus*2 -1;
        Vi_node = bus*2;
        Vr = V(bus*2-1);
        Vi= V(bus*2);
        c(2*bus-1) = Vr^2 + Vi^2 - mpc.bus(bus,12)^2;
        c(2*bus) = -Vr^2 - Vi^2 + mpc.bus(bus,13)^2;

        if nargout>2
           Gc(Vr_node, 2*bus-1) = 2*Vr;
           Gc(Vr_node, 2*bus) = -2*Vr;
           Gc(Vi_node, 2*bus-1) = 2*Vi;
           Gc(Vi_node, 2*bus) = -2*Vi;
        end
    end
    
    %for each branch
    branch_size = size(mpc.branch);
    for branch =1:branch_size(1)
        from_bus = mpc.branch(branch,1);
        to_bus = mpc.branch(branch,2);
        r=mpc.branch(branch,3);
        x=mpc.branch(branch,4);
        b=mpc.branch(branch,5);
        
        Vr_from_node = from_bus*2-1;
        Vr_to_node = to_bus*2-1;
        Vi_from_node = from_bus*2;
        Vi_to_node = to_bus*2;
        Vr_from = V(Vr_from_node);
        Vr_to =V(Vr_to_node);
        Vi_from = V(Vi_from_node);
        Vi_to = V(Vi_to_node);
        G_pu=r/(r^2 +x^2);
        B_pu = -x/(r^2 +x^2);        
        
        ratio = mpc.branch(branch,9);
        if ratio ==0.0
            %%% Transmission line model
                    


            Ir_from = G_pu*(Vr_from - Vr_to) - B_pu*(Vi_from-Vi_to) -b*Vi_from/2;
            Ir_to = G_pu*(Vr_to - Vr_from) - B_pu*(Vi_to-Vi_from)-b*Vi_to/2;
            Ii_from = G_pu*(Vi_from-Vi_to) + B_pu*(Vr_from-Vr_to) +b*Vr_from/2;
            Ii_to = G_pu*(Vi_to-Vi_from) +B_pu*(Vr_to-Vr_from) + b*Vr_to/2;

            dIr_from_dVr_from = G_pu;
            dIr_from_dVi_from = -B_pu-b/2;
            dIr_from_dVr_to = -G_pu;
            dIr_from_dVi_to = B_pu;
            
            dIr_to_dVr_from = -G_pu;
                dIr_to_dVr_to = G_pu;
                dIr_to_dVi_from = B_pu;
                dIr_to_dVi_to = -B_pu-b/2;
                
                dIi_from_dVr_from = B_pu +b/2;
                dIi_from_dVr_to = -B_pu;
                dIi_from_dVi_from = G_pu;
                dIi_from_dVi_to = -G_pu;
                
                dIi_to_dVr_from = -B_pu;
                dIi_to_dVr_to = B_pu+b/2;
                dIi_to_dVi_from= -G_pu;
                dIi_to_dVi_to = G_pu;
                
                Y(from_bus, from_bus) = Y(from_bus, from_bus) + dIr_from_dVr_from -1i*dIr_from_dVi_from;
                Y(from_bus, to_bus) = Y(from_bus, to_bus) + dIr_from_dVr_to - 1i*dIr_from_dVi_to;
                Y(to_bus, to_bus) = Y(to_bus, to_bus) + dIr_to_dVr_to - 1j*dIr_to_dVi_to; 
                Y(to_bus, from_bus) = Y(to_bus, from_bus) + dIr_to_dVr_from -1j*dIr_to_dVi_from;
                   
            ceq(Vr_from_node) = ceq(Vr_from_node)+Ir_from;
            ceq(Vi_from_node) = ceq(Vi_from_node)+Ii_from;
            ceq(Vr_to_node) = ceq(Vr_to_node)+Ir_to;
            ceq(Vi_to_node) = ceq(Vi_to_node)+Ii_to;
            
            if nargout>2
            
                dIr_from_dVr_from = G_pu;
                dIr_from_dVi_from = -B_pu-b/2;
                dIr_from_dVr_to = -G_pu;
                dIr_from_dVi_to = B_pu;
            
                dIr_to_dVr_from = -G_pu;
                dIr_to_dVr_to = G_pu;
                dIr_to_dVi_from = B_pu;
                dIr_to_dVi_to = -B_pu-b/2;
                
                dIi_from_dVr_from = B_pu +b/2;
                dIi_from_dVr_to = -B_pu;
                dIi_from_dVi_from = G_pu;
                dIi_from_dVi_to = -G_pu;
                
                dIi_to_dVr_from = -B_pu;
                dIi_to_dVr_to = B_pu+b/2;
                dIi_to_dVi_from= -G_pu;
                dIi_to_dVi_to = G_pu;
                


                
                Geq(Vr_from_node, Vr_from_node)= Geq(Vr_from_node, Vr_from_node)+ dIr_from_dVr_from;
                Geq(Vr_from_node, Vr_to_node) = Geq(Vr_from_node, Vr_to_node)+dIr_to_dVr_from;
                Geq(Vr_from_node, Vi_from_node) = Geq(Vr_from_node, Vi_from_node)+dIi_from_dVr_from;
                Geq(Vr_from_node, Vi_to_node) = Geq(Vr_from_node, Vi_to_node)+dIi_to_dVr_from;
                
                Geq(Vr_to_node, Vr_from_node)=Geq(Vr_to_node, Vr_from_node)+ dIr_from_dVr_to;
                Geq(Vr_to_node, Vr_to_node) = Geq(Vr_to_node, Vr_to_node)+dIr_to_dVr_to;
                Geq(Vr_to_node, Vi_from_node) = Geq(Vr_to_node, Vi_from_node)+dIi_from_dVr_to;
                Geq(Vr_to_node, Vi_to_node) =Geq(Vr_to_node, Vi_to_node)+ dIi_to_dVr_to;   
                
                Geq(Vi_from_node, Vr_from_node)= Geq(Vi_from_node, Vr_from_node)+dIr_from_dVi_from;
                Geq(Vi_from_node, Vr_to_node) =Geq(Vi_from_node, Vr_to_node)+ dIr_to_dVi_from;
                Geq(Vi_from_node, Vi_from_node) =Geq(Vi_from_node, Vi_from_node)+ dIi_from_dVi_from;
                Geq(Vi_from_node, Vi_to_node) =Geq(Vi_from_node, Vi_to_node)+ dIi_to_dVi_from;
                
                Geq(Vi_to_node, Vr_from_node)=Geq(Vi_to_node, Vr_from_node)+ dIr_from_dVi_to;
                Geq(Vi_to_node, Vr_to_node) =Geq(Vi_to_node, Vr_to_node)+ dIr_to_dVi_to;
                Geq(Vi_to_node, Vi_from_node) =Geq(Vi_to_node, Vi_from_node)+ dIi_from_dVi_to;
                Geq(Vi_to_node, Vi_to_node) = Geq(Vi_to_node, Vi_to_node)+dIi_to_dVi_to;
            end
            
        else
            %%% Transformer Model
            Gt = G_pu;
            Bt = B_pu+b/2;
            tr = ratio;
            tr2 = tr*tr;
            
            ang = mpc.branch(branch,10);
            
            if ang
               phi_rad =  ang*pi/180;
               cosphi = cos(phi_rad);
               sinphi =sin(phi_rad);
               Gcosphi = G_pu*cosphi;
               Gsinphi = G_pu*sinphi;
               Bcosphi = B_pu*cosphi;
               Bsinphi = B_pu*sinphi;
               
               G_shunt_from = G_pu/tr2;
               B_shunt_from = Bt/tr2;
               zero_factor=1;
               MR_from = zero_factor*(Gcosphi  - Bsinphi)/tr;
               MI_from = zero_factor*(Gsinphi  + Bcosphi)/tr;
               G_to = zero_factor*(Gcosphi + Bsinphi)/tr;
               B_to = zero_factor*(Bcosphi - Gsinphi)/tr;
               MR_to = zero_factor*Gt;
               MI_to = zero_factor*Bt;
        
               dIrfdVrf = G_shunt_from;
               dIrfdVrt = -MR_from;
               dIrfdVif = -B_shunt_from;
               dIrfdVit = MI_from;
        
               dIrtdVrf = -G_to;
               dIrtdVrt = MR_to;
               dIrtdVif = B_to;
               dIrtdVit = -MI_to;
        
               dIifdVrf = B_shunt_from;
               dIifdVrt = -MI_from;
               dIifdVif = G_shunt_from;
               dIifdVit = -MR_from;
        
               dIitdVrf = -B_to;
               dIitdVrt = MI_to;
               dIitdVif = -G_to;
               dIitdVit = MR_to;
   
               Ir_from = Vr_from*dIrfdVrf + Vr_to*dIrfdVrt + Vi_from*dIrfdVif + Vi_to*dIrfdVit;
               Ir_to = Vr_from*dIrtdVrf + Vr_to*dIrtdVrt + Vi_from*dIrtdVif + Vi_to*dIrtdVit;
               Ii_from = Vr_from*dIifdVrf + Vr_to*dIifdVrt + Vi_from*dIifdVif + Vi_to*dIifdVit;
               Ii_to = Vr_from*dIitdVrf + Vr_to*dIitdVrt + Vi_from*dIitdVif + Vi_to*dIitdVit;
               
               
               ceq(Vr_from_node) = ceq(Vr_from_node)+Ir_from;
               ceq(Vi_from_node) = ceq(Vi_from_node)+Ii_from;
               ceq(Vr_to_node) = ceq(Vr_to_node)+Ir_to;
               ceq(Vi_to_node) = ceq(Vi_to_node)+Ii_to;
               
               if nargout>2
                   Geq(Vr_from_node, Vr_from_node) = Geq(Vr_from_node, Vr_from_node) + dIrfdVrf;
                   Geq(Vr_from_node, Vr_to_node) = Geq(Vr_from_node, Vr_to_node) + dIrtdVrf;
                   Geq(Vr_from_node, Vi_from_node) = Geq(Vr_from_node, Vi_from_node) + dIifdVrf;
                   Geq(Vr_from_node, Vi_to_node) = Geq(Vr_from_node, Vi_to_node) + dIitdVrf;
                   
                   Geq(Vr_to_node, Vr_from_node) = Geq(Vr_to_node, Vr_from_node) + dIrfdVrt;
                   Geq(Vr_to_node, Vr_to_node) = Geq(Vr_to_node, Vr_to_node) + dIrtdVrt;
                   Geq(Vr_to_node, Vi_from_node) = Geq(Vr_to_node, Vi_from_node) + dIifdVrt;
                   Geq(Vr_to_node, Vi_to_node) = Geq(Vr_to_node, Vi_to_node) + dIitdVrt;
   
                   
                   Geq(Vi_from_node, Vr_from_node) = Geq(Vi_from_node, Vr_from_node) + dIrfdVif;
                   Geq(Vi_from_node, Vr_to_node) = Geq(Vi_from_node, Vr_to_node) + dIrtdVif;
                   Geq(Vi_from_node, Vi_from_node) = Geq(Vi_from_node, Vi_from_node) + dIifdVif;
                   Geq(Vi_from_node, Vi_to_node) = Geq(Vi_from_node, Vi_to_node) + dIitdVif;
                   
                   
                   Geq(Vi_to_node, Vr_from_node) = Geq(Vi_to_node, Vr_from_node) + dIrfdVit;
                   Geq(Vi_to_node, Vr_to_node) = Geq(Vi_to_node, Vr_to_node) + dIrtdVit;
                   Geq(Vi_to_node, Vi_from_node) = Geq(Vi_to_node, Vi_from_node) + dIifdVit;
                   Geq(Vi_to_node, Vi_to_node) = Geq(Vi_to_node, Vi_to_node) + dIitdVit;
                   
                   
               end
               
            else
                %no angle case
                zero_factor=1;
                G_series = zero_factor*G_pu/tr;
                B_series = zero_factor*B_pu/tr;
                G_shunt_from = zero_factor*(1-tr)/(tr2)*G_pu;
                B_shunt_from = zero_factor*((1-tr)/(tr2)*B_pu + (b/2)/tr2);
                G_shunt_to = zero_factor*G_pu*(1-1/tr);
                B_shunt_to = zero_factor*((1-1/tr)*B_pu + b/2);     
                
                if G_pu~=0
                    dIrfdVrf = G_series +G_shunt_from;
                    dIrfdVrt = -G_series;
                    
                    dIrtdVrf = -G_series;
                    dIrtdVrt = G_series + G_shunt_to;
                    
                    dIifdVif = G_series + G_shunt_from;
                    dIifdVit = -G_series;
                    
                    dIitdVif = -G_series;
                    dIitdVit = G_series + G_shunt_to;
                    
                    
                else
                    dIrfdVrf = 0;
                    dIrfdVrt = 0;
                    
                    dIrtdVrf = 0;
                    dIrtdVrt = 0;
                    
                    dIifdVif = 0;
                    dIifdVit = 0;
                    
                    dIitdVif = 0;
                    dIitdVit = 0;          
                    
                end
                
                dIrfdVif = -B_series -B_shunt_from;
                dIrfdVit = B_series;
                
                dIrtdVif = B_series;
                dIrtdVit = -B_series - B_shunt_to;
                
                dIifdVrf = B_series + B_shunt_from;
                dIifdVrt = -B_series;
                
                dIitdVrt = B_series + B_shunt_to;
                dIitdVrf = -B_series;
                     
               Ir_from = Vr_from*dIrfdVrf + Vr_to*dIrfdVrt + Vi_from*dIrfdVif + Vi_to*dIrfdVit;
               Ir_to = Vr_from*dIrtdVrf + Vr_to*dIrtdVrt + Vi_from*dIrtdVif + Vi_to*dIrtdVit;
               Ii_from = Vr_from*dIifdVrf + Vr_to*dIifdVrt + Vi_from*dIifdVif + Vi_to*dIifdVit;
               Ii_to = Vr_from*dIitdVrf + Vr_to*dIitdVrt + Vi_from*dIitdVif + Vi_to*dIitdVit;
               
              
                Y(from_bus, from_bus) = Y(from_bus, from_bus) + dIrfdVrf -1i*dIrfdVif;
                Y(from_bus, to_bus) = Y(from_bus, to_bus) + dIrfdVrt - 1i*dIrfdVit;
                Y(to_bus, to_bus) = Y(to_bus, to_bus) + dIrtdVrt - 1j*dIrtdVit; 
                Y(to_bus, from_bus) = Y(to_bus, from_bus) + dIrtdVrf -1j*dIrtdVif;
               
               ceq(Vr_from_node) = ceq(Vr_from_node)+Ir_from;
               ceq(Vi_from_node) = ceq(Vi_from_node)+Ii_from;
               ceq(Vr_to_node) = ceq(Vr_to_node)+Ir_to;
               ceq(Vi_to_node) = ceq(Vi_to_node)+Ii_to;
               
               if nargout>2
                   Geq(Vr_from_node, Vr_from_node) = Geq(Vr_from_node, Vr_from_node) + dIrfdVrf;
                   Geq(Vr_from_node, Vr_to_node) = Geq(Vr_from_node, Vr_to_node) + dIrtdVrf;
                   Geq(Vr_from_node, Vi_from_node) = Geq(Vr_from_node, Vi_from_node) + dIifdVrf;
                   Geq(Vr_from_node, Vi_to_node) = Geq(Vr_from_node, Vi_to_node) + dIitdVrf;
                   
                   Geq(Vr_to_node, Vr_from_node) = Geq(Vr_to_node, Vr_from_node) + dIrfdVrt;
                   Geq(Vr_to_node, Vr_to_node) = Geq(Vr_to_node, Vr_to_node) + dIrtdVrt;
                   Geq(Vr_to_node, Vi_from_node) = Geq(Vr_to_node, Vi_from_node) + dIifdVrt;
                   Geq(Vr_to_node, Vi_to_node) = Geq(Vr_to_node, Vi_to_node) + dIitdVrt;
   
                   
                   Geq(Vi_from_node, Vr_from_node) = Geq(Vi_from_node, Vr_from_node) + dIrfdVif;
                   Geq(Vi_from_node, Vr_to_node) = Geq(Vi_from_node, Vr_to_node) + dIrtdVif;
                   Geq(Vi_from_node, Vi_from_node) = Geq(Vi_from_node, Vi_from_node) + dIifdVif;
                   Geq(Vi_from_node, Vi_to_node) = Geq(Vi_from_node, Vi_to_node) + dIitdVif;
                   
                   
                   Geq(Vi_to_node, Vr_from_node) = Geq(Vi_to_node, Vr_from_node) + dIrfdVit;
                   Geq(Vi_to_node, Vr_to_node) = Geq(Vi_to_node, Vr_to_node) + dIrtdVit;
                   Geq(Vi_to_node, Vi_from_node) = Geq(Vi_to_node, Vi_from_node) + dIifdVit;
                   Geq(Vi_to_node, Vi_to_node) = Geq(Vi_to_node, Vi_to_node) + dIitdVit;
           
               end
                
                
            end
            
        end

    
    end

    bus_size = size(mpc.bus);
    %%% GENERATOR MODEL
    gen_size = size(mpc.gen);
    for gen =1:gen_size(1)
        from_bus = mpc.gen(gen,1);
        Vr_node = from_bus*2-1;
        Vi_node = from_bus*2;
        P_node = bus_size(1)*2+gen*2-1;
        Q_node = bus_size(1)*2+gen*2;
        
        Vr = V(Vr_node);
        Vi = V(Vi_node);
        
        Pg = V(P_node);
        Qg = V(Q_node);
    
        Ir = (-Pg*Vr+Qg*Vi)/(Vr^2 + Vi^2);
        Ii = (-Pg*Vi - Qg*Vr)/(Vr^2 +Vi^2);
        
        ceq(Vr_node)=ceq(Vr_node) +Ir;
        ceq(Vi_node)=ceq(Vi_node)+Ii;
        
        if nargout>2
          
            dIr_dVr = (-Pg*(Vi^2 - Vr^2) - 2*Qg*Vr*Vi)/((Vr^2+Vi^2)^2);
            dIr_dVi = (Qg*(Vr^2-Vi^2) + 2*Pg*Vr*Vi)/((Vr^2+Vi^2)^2);
            
            dIr_dP = -Vr/(Vr^2-Vi^2);
            dIr_dQ = Vi/(Vr^2-Vi^2);
            
            dIi_dVr = dIr_dVi;
            dIi_dVi=-dIr_dVr;
            
            dIi_dP = -Vi/(Vr^2 + Vi^2);
            dIi_dQ = -Vr/(Vr^2 + Vi^2);
            
            Geq(Vr_node, Vr_node) =Geq(Vr_node, Vr_node)+ dIr_dVr;
            Geq(Vr_node, Vi_node) =Geq(Vr_node, Vi_node)+ dIi_dVr;
            
            Geq(Vi_node, Vr_node) =Geq(Vi_node, Vr_node)+ dIr_dVi;
            Geq(Vi_node, Vi_node) = Geq(Vi_node, Vi_node) +dIi_dVi;
            
            Geq(P_node, Vr_node) = Geq(P_node, Vr_node) + dIr_dP;
            Geq(P_node, Vi_node) = Geq(P_node, Vi_node) + dIi_dP;
            
            Geq(Q_node, Vr_node) = Geq(Q_node, Vr_node) + dIr_dQ;
            Geq(Q_node, Vi_node) = Geq(Q_node, Vi_node) + dIi_dQ;
            
        end
        
    end    
    
    
    
    %%%LOAD MODEL
    bus_size = size(mpc.bus);
    for bus =  1:bus_size(1)
        from_bus=mpc.bus(bus,1);
        Pd = mpc.bus(bus,3)/100;
        Qd = mpc.bus(bus, 4)/100;
        Vr_node = from_bus*2-1;
        Vi_node = from_bus*2;
        Vr = V(Vr_node);
        Vi = V(Vi_node);     
        Ir = (Pd*Vr+Qd*Vi)/(Vr^2 +Vi^2);
        Ii = (Pd*Vi -Qd*Vr)/(Vr^2 +Vi^2);
        
        ceq(Vr_node)=ceq(Vr_node) +Ir;
        ceq(Vi_node)=ceq(Vi_node)+Ii;
        
        if nargout>2

            Vr2Vi2 = Vr^2 - Vi^2;
            Vmag4 = (Vr^2 + Vi^2)^2;
            dIr_dVr = -(Pd*Vr2Vi2 + 2.0*Qd*Vr*Vi)/Vmag4;
            dIr_dVi = (Qd*Vr2Vi2 - 2.0*Pd*Vr*Vi)/Vmag4;
            
            dIi_dVr = dIr_dVi;
            dIi_dVi = (Pd*Vr2Vi2 + 2.0*Qd*Vr*Vi)/Vmag4;
            
            
            Geq(Vr_node, Vr_node) =Geq(Vr_node, Vr_node)+ dIr_dVr;
            Geq(Vr_node, Vi_node) =Geq(Vr_node, Vi_node)+ dIi_dVr;
            
            Geq(Vi_node, Vr_node) =Geq(Vi_node, Vr_node)+ dIr_dVi;
            Geq(Vi_node, Vi_node) = Geq(Vi_node, Vi_node) +dIi_dVi;
            
        end
    end
    
    
    %%% SHUNT MODEL
    bus_size = size(mpc.bus);
    for bus =  1:bus_size(1)
        from_bus=mpc.bus(bus,1);
        Bsh = mpc.bus(bus, 6)/100;
        Vr_node = from_bus*2-1;
        Vi_node = from_bus*2;
        Vr = V(Vr_node);
        Vi = V(Vi_node); 
        ceq(Vr_node) = ceq(Vr_node) - Vi*Bsh;
        ceq(Vi_node) = ceq(Vi_node) + Vr*Bsh;
        
        Y(bus,bus) = Y(bus,bus) + Bsh*1i;
        if nargout>2
            Geq(Vr_node, Vi_node)= Geq(Vr_node, Vi_node)+Bsh;
            Geq(Vi_node, Vr_node) =Geq(Vi_node, Vr_node) -Bsh;
            
        end
    end
    
 
end



