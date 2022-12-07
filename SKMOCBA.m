%% Copyright 2020 Sebastian Rojas Gonzalez <srojas33@gmail.com><sebastian.rojasgonzalez@ugent.be>
%{
Source code for the paper: 
%--------------------------------------------------------------------------
Gonzalez, S. R., Jalali, H., & Van Nieuwenhuyse, I. (2020). A multiobjective 
stochastic simulation optimization algorithm. European Journal of Operational 
Research, 284(1), 212-226.
%--------------------------------------------------------------------------

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the above copyright notice is
retained.
%}

function SKMOCBA(Global)

    %% Pre-process data
    %Test function
    test_function = func2str(Global.problem);
    
    %LHS size
    lhs_size = 11*Global.D-1;
    
    %%Import LHS and dicretized search space 
    %LHS: Pop_lhs
    %Design space:  seq_set
    if strcmp(test_function,'DTLZ7')
        load('DTLZ7a_1'); %#ok<LOAD>
    else
        load('ZDT1a_1'); %#ok<LOAD>
    end
    %Note: The LHS and design space are scenarios previously generated,
    %in order to control the size of the true Pareto set. 
    
    
    %% Noise parameters 
    %Level: 1 =  low, 2 = medium, 3 = high
    level = 1;
    %Linear case: 1 = best case, 2 = worst case
    caso = 1;
    
    %% Replication budget
    B = 50; 
    b = 25;
    Bmax = 100;
        
    %% Initialize population
    Population = INDIVIDUAL(repmat(Global.upper-Global.lower,[lhs_size,1]).*Pop_lhs+repmat(Global.lower,[lhs_size,1]));
    
    %% Design space      
    design_space = seq_set.decs;
    prt3 = ['Design space size = ', num2str(size(seq_set,2))];
    disp(prt3);
    %pause;
    
    %All candidate points 
    Candidates = [seq_set, Population];
    
    %True PS
    design_nondom = NDSort(Candidates.objs,1);
    design_nondom_size = sum(design_nondom(:) == 1);
    print = ['Size of design space non-dominated set = ',num2str(design_nondom_size)];
    disp(print); 
    %pause;
    
    
    %% Initial phase
    
    %Initialize objective and variance matrices
    Cell_Obj_rep = cell(Global.evaluation+lhs_size); %Matrix all reps
    Mat_Obj = zeros(Global.evaluation+lhs_size,Global.M);
    Mat_Var = zeros(Global.evaluation+lhs_size,Global.M);
    Mat_Obj_scal = zeros(Global.evaluation+lhs_size,1);
    Mat_Var_scal = zeros(Global.evaluation+lhs_size,1);
    
    %Heterogeneous noise
    constants = heter_noise(Global, Candidates, level, caso,B);
    
    %Simulate (expensive) objectives on LHS design points
    for i = 1 : lhs_size
        det_obj = Population(1, i).obj;
        Mat_Obj_rep = zeros(B,Global.M);
        for k = 1 : Global.M
            a = constants(k,1);
            t = constants(k,2);
            for j = 1 : B
                Mat_Obj_rep(j,k) = det_obj(k)+normrnd(0,(a*det_obj(k)+a*t));
            end
            Mat_Obj(i,k) = mean(Mat_Obj_rep(:,k)); %Save mean of the sample
            Mat_Var(i,k) = var(Mat_Obj_rep(:,k))/B; %Save variance on the mean
        end
        Cell_Obj_rep{i} = Mat_Obj_rep;
    end
            
    %% Weight vectors
    n_weights = Global.evaluation;
	[Weights,Global.N] = UniformPoint(n_weights,Global.M);
    W = unique(Weights,'rows');
    
    %% Optimization
    %Record computing time
    tstart = tic;
    PopSize  = lhs_size;
    iter = 0;
    Global.evaluated = 0;
    while Global.NotTermination(Population) 
        eval = ['Number of evaluations = ',num2str(size(Population,2))];
        disp(eval);
        
        %% Scalarization
        
        %Randomly select a weight vector
        lambda  = W(randi(size(W,1)),:); 
        [N,~]  = size(Population.decs);  
        
        %Scalarize objectives per point sampled
        for i = 1:N
            Mat_Obj_rep_i =  Cell_Obj_rep{i};                      
            rep_size = size(Mat_Obj_rep_i,1);
            PCheby_i = max(Mat_Obj_rep_i.*repmat(lambda,[rep_size,1]),[],2)+0.05.*sum(Mat_Obj_rep_i.*repmat(lambda,[rep_size,1]),2);
            Mat_Obj_scal(i) = mean(PCheby_i);
            Mat_Var_scal(i) = var(PCheby_i)/rep_size;
        end
       
       %Normalize objectives
       Mat_Obj_scal = Mat_Obj_scal(1:N,:);
       PCheby = (Mat_Obj_scal-repmat(min(Mat_Obj_scal,[],1),[N,1]))./repmat((max(Mat_Obj_scal,[],1)-min(Mat_Obj_scal,[],1)),[N,1]);
       NCheby = Mat_Var_scal(1:N,:);       
       PDec   = Population.decs;
        
       %% Bayesian optimization using SK metamodels
        dmodel     = SKfit(PDec,PCheby,ones(PopSize,1),NCheby,2);
        
        %Search infill point
        [PopDec] = SearchALG(PCheby,Population.decs,dmodel, design_space);
        Population = [Population,INDIVIDUAL(PopDec)]; %#ok<AGROW>
        PopSize = PopSize + 1;
        
        %Simulate objectives in the new point
        det_obj = Population(1, PopSize).obj;
        Mat_Obj_rep = zeros(b,Global.M);
        for k = 1 : Global.M
            a = constants(k,1);
            t = constants(k,2);
            for j = 1 : b
                Mat_Obj_rep(j,k) = det_obj(k)+normrnd(0,(a*det_obj(k)+a*t));
            end
            Mat_Obj(PopSize,k) = mean(Mat_Obj_rep(:,k)); %Save sample mean
            Mat_Var(PopSize,k) = var(Mat_Obj_rep(:,k))/b; %Save sample variance
        end
        Cell_Obj_rep{PopSize} = Mat_Obj_rep;
        
        iter = iter + 1;
        
        %Remove infill point from design space
        design_space(ismember(design_space,PopDec,'rows'),:)=[];
        
       
       %%  Save data, compute metrics and plot %%
        if iter == Global.evaluation
            
            N = Global.evaluation+lhs_size;
            
          %% Results without MOCBA
            elapsed_time = toc(tstart);
            disp('Results without MOCBA');
            prt3 = ['Design space size = ', num2str(size(seq_set,2))]; disp(prt3);
            print = ['Size of design space Pareto set = ',num2str(design_nondom_size)]; disp(print);
            
            %Performance metrics
            metrics(Global, Population, Mat_Obj, Candidates);
            plots(Global,Population, test_function, Mat_Obj, Candidates);
            pause;

          %% MORS procedure - MOCBA 
            if b < B
                tstart2 = tic;
                disp('Running MOCBA...');
                % Parameters and initialization
                ind = (1:N)';
                Mat_Obj_ind = [Mat_Obj,ind];
                Mat_Var_ind = [Mat_Var,ind];
                total_budget = (B-b)*Global.evaluation;

                %Main loop
                while total_budget > 0
                    points = size(Mat_Obj_ind,1);
                    Mat_Obj_temp = Mat_Obj_ind(1:points,1:Global.M);
                    Mat_Var_temp = Mat_Var_ind(1:points,1:Global.M);                    
                    alphas = mocba(Mat_Obj_temp',Mat_Var_temp');
                    budget = round(alphas*total_budget);
                    if all(budget == 0)
                        [~,max_ind] = max(alphas);
                        budget(max_ind) = 1;
                    end
                    for q = 1:points
                        i = Mat_Obj_ind(q,end);
                        num_reps = size(Cell_Obj_rep{i},1);
                        reps = budget(q);
                        if (num_reps < Bmax) && (reps ~= 0)
                            total_reps = reps + num_reps;
                            if total_reps <= Bmax
                                run_reps = reps;
                            elseif total_reps > Bmax
                                run_reps = Bmax-num_reps;
                            end
                            det_obj = Population(1, i).obj;
                            obj_temp = zeros(run_reps,Global.M);
                            Mat_Obj_i = [Cell_Obj_rep{i};obj_temp];    
                            for k = 1 : Global.M
                                a = constants(k,1);
                                t = constants(k,2); 
                                r1 = num_reps+1;
                                r2 = num_reps+run_reps;
                                for j = r1 : r2
                                    Mat_Obj_i(j,k) = det_obj(k)+normrnd(0,(a*det_obj(k)+a*t));
                                end
                                Mat_Obj(i,k) = mean(Mat_Obj_i(:,k)); %Save new mean response
                            end
                            Cell_Obj_rep{i} = Mat_Obj_i;
                            total_budget = total_budget - run_reps;
                        end
                    end
                    Mat_Obj_ind_temp = Mat_Obj_ind;
                    Mat_Var_ind_temp = Mat_Var_ind;
                    for w = 1:points
                        num_reps = size(Cell_Obj_rep{w},1);
                        if num_reps >= Bmax  
                            Mat_Obj_ind_temp(w,:) = inf;
                            Mat_Var_ind_temp(w,:) = inf;
                        end
                    end
                    Mat_Obj_ind = Mat_Obj_ind_temp(Mat_Obj_ind_temp(:,1)~=inf,:);
                    Mat_Var_ind = Mat_Var_ind_temp(Mat_Var_ind_temp(:,1)~=inf,:);      
                end
                
                disp('Results with MOCBA')
                
                prt3 = ['Design space size = ', num2str(size(seq_set,2))]; disp(prt3);
                print = ['Size of design space Pareto set = ',num2str(design_nondom_size)]; disp(print);
            
                %Performance metrics
                metrics(Global, Population, Mat_Obj, Candidates);
                plots(Global,Population, test_function, Mat_Obj, Candidates);
                %pause;
                
                elapsed_time2 = toc(tstart2);
                elapsed_time = elapsed_time + elapsed_time2;
                print_time = ['Elapsed time = ', num2str(elapsed_time),'s'];
                disp(print_time);
            end
        end
    end
end
