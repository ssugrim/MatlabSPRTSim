%find the maximal params index
Alpha_bound = 0.10;
Beta_bound = 0.05;

Alpha_mask = Data_SPRT_error_1 < Alpha_bound;
Beta_mask = Data_SPRT_error_2 < Beta_bound;
Combo_mask = Alpha_mask & Beta_mask;
found_max = max(max(max(max(Data_SPRT_channels_found_H0(Combo_mask)))))
num_found = sum(sum(sum(Data_SPRT_channels_found_H0 == found_max)))
[max_chan,max_lam,max_del,max_split] = ind2sub(size(Data_SPRT_channels_found_H0),find(Data_SPRT_channels_found_H0 == found_max));

for i = 1:num_found
    sprintf('inidicies chan:%d lamb:%d del%d split%d',max_chan(i),max_lam(i),max_del(i),max_split(i))
    sprintf('chan:%2.2f lamb:%2.2f del:%2.2f split:%2.2f',Channels(max_chan(i)),lambda_range(max_lam(i)),delta_range(max_del(i)),split_range(max_split(i)))
end




% plotStyle = {'b+-','ro-','gs-','md-','k^-','c*-'};
% 
% f = figure(1);
% hold;
% clear legendInfo;
% for i = 1:length(lambda_range)
%     plot(delta_range,squeeze(Data_SPRT_channels_found_H0(max_chan,i,:,max_split)),plotStyle{mod(i,length(plotStyle))+1},'LineWidth',2);
%     legendInfo{i} = strcat('\lambda',sprintf('= %2.2f%',lambda_range(i)));
% end
% xlabel('\delta, Width of the indifference region');
% ylabel('Number of channels found');
% legend(legendInfo);
% set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold');
% 
% f = figure(2);
% hold;
% clear legendInfo;
% for i = 1:length(lambda_range)
%     plot(delta_range,squeeze(Data_SPRT_mean_channels_passed_on(max_chan,i,:,max_split)),plotStyle{mod(i,length(plotStyle))+1},'LineWidth',2);
%     legendInfo{i} = strcat('\lambda',sprintf('= %2.2f%',lambda_range(i)));
% end
% xlabel('\delta, Width of the indifference region');
% ylabel('Number of channels Passed on');
% legend(legendInfo);
% set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold');
% 
% f = figure(3);
% hold;
% clear legendInfo;
% for i = 1:length(lambda_range)
%     plot(delta_range,squeeze(Data_SPRT_error_2(max_chan,i,:,max_split)),plotStyle{mod(i,length(plotStyle))+1},'LineWidth',2);
%     legendInfo{i} = strcat('\lambda',sprintf('= %2.2f%',lambda_range(i)));
% end
% xlabel('\delta, Width of the indifference region');
% ylabel('Probabilty of Missclassification');
% legend(legendInfo);
% set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold');
% 
% f = figure(4);
% clear legendInfo;
% plot(delta_range,squeeze(Data_SPRT_mean_samples_make_decision(max_chan,max_lam,:,max_split)));
% xlabel('\delta, Width of the indifference region');
% ylabel('Samples required to categorize one channel');
% set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold');
% 
% f = figure(5);
% clear legendInfo;
% hold;
% plot(delta_range,squeeze(Data_SPRT_error_2(max_chan,max_lam,:,max_split)),plotStyle{mod(1,length(plotStyle))+1},'LineWidth',2);
% legendInfo{i} = sprintf('Type II error');
% plot(delta_range,squeeze(Data_SPRT_error_1(max_chan,max_lam,:,max_split)),plotStyle{mod(2,length(plotStyle))+1},'LineWidth',2);
% legendInfo{i} = sprintf('Type I error');
% xlabel('\delta, Width of the indifference region');
% ylabel('Error Probability');
% set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold');
% 
% f = figure(6);
% clear legendInfo;
% plot(lambda_range,squeeze(Data_SPRT_mean_samples_make_decision(max_chan,:,max_del,max_split)));
% xlabel('\lambda, Sample cost coefficent');
% ylabel('Samples required to categorize one channel');
% set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold');
% 
% f = figure(7);
% clear legendInfo;
% hold;
% plot(lambda_range,squeeze(Data_SPRT_error_2(max_chan,:,max_del,max_split)),plotStyle{mod(1,length(plotStyle))+1},'LineWidth',2);
% legendInfo{i} = sprintf('Type II error');
% plot(lambda_range,squeeze(Data_SPRT_error_1(max_chan,:,max_del,max_split)),plotStyle{mod(2,length(plotStyle))+1},'LineWidth',2);
% legendInfo{i} = sprintf('Type I error');
% xlabel('\delta, Width of the indifference region');
% ylabel('Error Probability');
% set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold');