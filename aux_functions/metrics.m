function metrics( Global, Population, Mat_Obj, Candidates)

%% Retrieve data
[~, ~, ~, ~, PF_obs, ~, ~, ~, ~, ~, ~, type1_size, ~, type2_size, ~, NonDominated_obs_size, PS_sampled_size, PS_s, PS_ID, ...
    APS, prec, rec, f1 ] = sets(Global, Population, Mat_Obj, Candidates);

    %% Quality indicators
     %Hypervolume of observed front
     obs_hv = HV(PF_obs,Global.PF);
     %obs_hv2 = stk_dominatedhv(PF_obs,ref); 

     %IGD of observed front
     obs_igd = IGD(PF_obs,Global.PF);

     %Print indicators
     hv_print = ['HV = ', num2str(obs_hv)];
     disp(hv_print);
     %hv_print2 = ['HV2 = ', num2str(obs_hv2)];
     %disp(hv_print2);
     igd_print = ['IGD = ', num2str(obs_igd)];
     disp(igd_print);
     %pause;
    %}
    
    %% Search metrics
    %prt3 = ['Design space size = ', num2str(size_set)];
    %disp(prt3);
    
    %print = ['Size of design space non-dominated set = ',num2str(design_nondom_size)];
    %disp(print);

    test = ['Number of points in observed PF = ',num2str(NonDominated_obs_size)];
     disp(test); 
     
    print1 = ['Number of true Pareto points sampled = ',num2str(PS_sampled_size)];
    disp(print1);
    
    result1 = ['PS_s = ',num2str(PS_s) ,'%'];
    disp(result1);

     prt_t1 = ['Number of MCE error points = ',num2str(type1_size)];
     disp(prt_t1);
     
      prt_t2 = ['Number of MCI error points = ',num2str(type2_size)];
     disp(prt_t2);
     
     psid = ['PS_ID = ',num2str(PS_ID),'%'];
     disp(psid); 
     
    metric1 = ['APS = ', num2str(APS)];
    disp(metric1);
    
    metric2 = ['precision = ', num2str(prec)];
    disp(metric2);
    
    metric3 = ['recall = ', num2str(rec)];
    disp(metric3);
    
    metric4 = ['f1 = ', num2str(f1)];
    disp(metric4);
end

