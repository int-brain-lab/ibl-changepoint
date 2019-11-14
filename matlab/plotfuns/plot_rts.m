tt = 1:size(data.tab,1)';

for f = 1:2
    figure(f)
    switch f
        case 1; p_idx = data.p_true == 0.2;
        case 2; p_idx = data.p_true == 0.8;
    end
    
    for i = 1:9
        subplot(3,3,i);
        idx = (data.signed_contrasts_idx == i)' & data.resp_correct & p_idx;
        yy = data.tab(idx,end); 
        histogram(log10(yy(yy>0)),'BinWidth',0.1);
        hold on;
        idx = (data.signed_contrasts_idx == i)' & ~data.resp_correct & p_idx;
        yy = data.tab(idx,end); 
        histogram(log10(yy(yy>0)),'BinWidth',0.1);
        xlim([-log10(100) log10(60)]);
        title(['contrast = ' num2str(signed_contrasts_vec(i))]);
        box off;
    end
end

figure(3)
yy = data.tab(:,end);
yy(yy < 0.01) = 0.01;
t_vec = tt;
plot(t_vec,yy); hold on;
xlim([1 min(1e4,tt(end))]);
ylim([0.01,60]);
idx = find(data.tab(:,1)==1) - 1;
for i = 1:numel(idx)
    plot(idx(i)*[1 1],[0.01 60],'k--','LineWidth',1);
end
% set(gca,'Yscale','log');
    
%     for i = 1:9
%         subplot(3,3,i);
%         idx = (data.signed_contrasts_idx == i)' & data.resp_correct & p_idx;
%         yy = data.tab(idx,end); t_vec = tt(idx);
%         plot(t_vec(yy>0),log10(yy(yy>0)));
%         hold on;
%         idx = (data.signed_contrasts_idx == i)' & ~data.resp_correct & p_idx;
%         yy = data.tab(idx,end); t_vec = tt(idx);
%         plot(t_vec(yy>0),log10(yy(yy>0)));
%         xlim([1 tt(end)]);
%         title(['contrast = ' num2str(signed_contrasts_vec(i))]);
%         box off;
%     end
