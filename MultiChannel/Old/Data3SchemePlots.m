Ratios = Channels ./ deadline;
figure(1);
s(1) = subplot(2,1,1);
s(2) = subplot(2,1,2);
% s(3) = subplot(3,1,3);

hold(s(1));

plot(s(1),Ratios,mean(Data_SPRT_error_2'),'--b');
% plot(s(1),Ratios,mean(Data_Simple_error_2'),'-xr');
plot(s(1),Ratios,mean(Data_Tree_error_2'),'-.g');

ylabel(s(1),'Error Rate (Beta)');
xlim(s(1),[min(Ratios) max(Ratios)]);

% hold(s(2));
% plot(s(2),Ratios,mean(Data_SPRT_mean_m_k'),'--b');
% plot(s(2),Ratios,mean(Data_Simple_mean_m_k'),'-xr');
% plot(s(2),Ratios,mean(Data_Tree_mean_m_k'),'-.g');
% 
% ylabel(s(2),'Mean Samples per Channel');
% xlim(s(2),[min(Ratios) max(Ratios)]);

hold(s(2));
plot(s(2),Ratios,mean(Data_SPRT_percent_classified'),'--b');
% plot(s(2),Ratios,mean(Data_Simple_percent_classified'),'-xr');
plot(s(2),Ratios,mean(Data_Tree_percent_classified'),'-.g');

ylabel(s(2),'Fraction of Channels classified.');
xlim(s(2),[min(Ratios) max(Ratios)]);

xlabel(s(2),'Ratio of Channels to Measurement Limit');

legend('SPRT','Tree');
% legend('SPRT','Simple','Tree');

figure(2);
sprt_hist = mean(squeeze(Data_SPRT_histogram(length(Channels),:,:)));
% simple_hist = mean(squeeze(Data_Simple_histogram(length(Channels),:,:)));
tree_hist = mean(squeeze(Data_Tree_histogram(length(Channels),:,:)));
%bar(Histo_bin,[sprt_hist',tree_hist',simple_hist']);
bar(Histo_bin,[sprt_hist',tree_hist']);
%ylim([0,45]);
xlim([1,30]);
% legend('SPRT','Tree','Simple');
legend('SPRT','Tree');
xlabel('Samples Allocated')
ylabel('Number of Channels with that Allocation');