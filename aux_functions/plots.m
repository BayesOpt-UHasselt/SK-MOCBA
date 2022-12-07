function plots(Global,Population,test_function, Mat_Obj, Candidates)

%Retrieve data
[~,~,~, ~, PF_obs, nonPF_obs, ~, ~, PF_pop, nonPF_pop, type1_obs_obj,type1_size, type2_obs_obj,type2_size,~,~,~,~,~,~,...
    ~,~,~ ] = sets(Global, Population, Mat_Obj, Candidates);

if Global.M == 2
        %Axis
        if strcmp(test_function,'DTLZ7')
            size_axis = [0 0.91 2.3 4.1];
        elseif (strcmp(test_function,'WFG1')) || (strcmp(test_function,'WFG3')) || (strcmp(test_function,'WFG4'))
            size_axis = [0 3 0 5];
        else
            size_axis = [-0.1 1.1 -0.1 1.6];
        end
        
         % Plot observed means 
         %{
         figure;
         Draw(PF_obs);
         axis(size_axis);
         title(sprintf('%s on %s(observed means)',func2str(Global.algorithm),func2str(Global.problem)),'Interpreter','none');
         pause;
         %}

         %PFobs + MCE + MCI
         figure;
         plot(PF_obs(:,1),PF_obs(:,2),'o','MarkerSize',10,'MarkerEdgeColor',[.4 .4 .4],'LineWidth',1,'MarkerFaceColor',[0.7,0.7,0.7]);
         hold on
         plot(nonPF_obs(:,1),nonPF_obs(:,2),'o','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
         plot(type1_obs_obj(:,1),type1_obs_obj(:,2),'^','MarkerSize',10,'MarkerEdgeColor','b','LineWidth',1,'MarkerFaceColor','b');
         plot(type2_obs_obj(:,1),type2_obs_obj(:,2),'s','MarkerSize',10,'MarkerEdgeColor','r','LineWidth',1,'MarkerFaceColor','r');
         titl = ['PF vs. MCE and MCI on ',test_function];
         hTitle = title(titl);
         set(hTitle,'FontSize',20)
         set(gca,'FontSize',15);
         axis(size_axis);
         if (type1_size == 0) && (type2_size ~= 0)
             legend('PF','Dominated points', 'MCI');
         elseif (type1_size ~= 0) && (type2_size == 0)
             legend('PF','Dominated points', 'MCE');
         else
             legend('PF','Dominated points', 'MCE', 'MCI');
         end
         xlabel('f_1'); ylabel('f_2');
         hold off
         %pause;

         %PF sampled vrs PF obs 
         %
         figure;
         plot(PF_pop(:,1),PF_pop(:,2),'o','MarkerSize',10,'MarkerEdgeColor',[0 0 0],'LineWidth',1,'MarkerFaceColor',[0,0,0]);
         hold on
         plot(nonPF_pop(:,1),nonPF_pop(:,2),'o','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
         axis(size_axis);
         titl = ['True performance on ',test_function];
         hTitle = title(titl);
         set(hTitle,'FontSize',20)
         set(gca,'FontSize',15);
         legend('Pareto front','Dominated points','Location','northeast');
         xlabel('f_1'); ylabel('f_2');
         hold off
         %pause;
         %

         %True PF vrs PF sampled vrs PFobs 
         %{
         figure;
         plot(PF_pop(:,1),PF_pop(:,2),'o','MarkerSize',12,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','k');
         hold on
         plot(true_PF(:,1),true_PF(:,2),'o','MarkerSize',12,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
         %plot(PF_obs(:,1),PF_obs(:,2),'o','MarkerSize',12,'MarkerEdgeColor',[.4 .4 .4],'LineWidth',1,'MarkerFaceColor',[0.7,0.7,0.7]);
         title('True PF vrs PF sampled vrs PFobs');
         axis(size_axis);
         legend('Sampled PF','True PF', 'PF_{obs}');
         hold off
         %pause;
         %}

         %PF_obs vrs true obj values 
         %{
         figure;
         plot(PF_obs_true(:,1),PF_obs_true(:,2),'o','MarkerSize',8,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','k');
         hold on
         %plot(true_PF(:,1),true_PF(:,2),'o','MarkerSize',8,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
         plot(PF_obs(:,1),PF_obs(:,2),'o','MarkerSize',8,'MarkerEdgeColor',[.4 .4 .4],'LineWidth',1,'MarkerFaceColor',[0.7,0.7,0.7]);
         title('True values of PF_{obs} vrs. PF_{obs}');
         axis(size_axis);
         legend('True values of PF_{obs}','PF_{obs}');
         xlabel('f^1'); ylabel('f^2');
         hold off
         pause;
         %}
end


end

