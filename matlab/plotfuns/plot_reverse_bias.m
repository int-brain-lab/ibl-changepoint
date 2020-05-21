function [resp_match,tot_resp,model_match] = plot_reverse_bias(data,side,plotflag,p_model)

if nargin < 2; side = []; end
if nargin < 3 || isempty(plotflag); plotflag = true; end
if nargin < 4; p_model = []; end

if isempty(side) && nargout == 0
    
    close all;
    figure(1);
    subplot(1,2,1);
    plot_reverse_bias(data,1);
    subplot(1,2,2);
    plot_reverse_bias(data,-1);    
    plotflag = false;

elseif iscell(data)
    for iData = 1:numel(data)        
        try            
            dataset = data{iData};
            if ischar(dataset); dataset = read_data_from_csv(dataset); end
            
            params = params_new('changepoint_runlength',dataset);
            params.runlength_min = 1;            
%            params.runlength_min = 20;
            [~,output] = changepoint_bayesian_nll(params,dataset,1);
            p_model = output.p_estimate;
            [resp_match(iData,:),tot_resp(iData,:),model_match(iData,:)] = plot_reverse_bias(dataset,side,false,p_model);
        catch
            fprintf('Error with dataset %d, ignoring.\n',iData);
            resp_match(iData,:) = NaN(1,99);
            tot_resp(iData,:) = NaN(1,99);
            model_match(iData,:) = NaN(1,99);
        end
    end
    
    
else
    X = data.tab;
    X(:,4) = X(:,4).*sign(X(:,5));      % Convert contrast to signed contrast

    sessions = unique(X(:,2))';

    resp_match = zeros(1,99);
    tot_resp = zeros(1,99);
    model_match = zeros(1,99);
    
    for iSession = sessions
        
        Xs = X(X(:,2) == iSession,:);   % Individual session data
        Nt = size(Xs,1);                % tot # trials in the session    
        trial_block = zeros(Nt,1);      % trial number inside block

        block_starts = find(diff([0;Xs(:,3);0]))';

        for iBlock = 1:numel(block_starts)-1
            idx = block_starts(iBlock):(block_starts(iBlock+1)-1);
            trial_block(idx) = 1:numel(idx);
        end

        % Find zero-contrast trials
        idx0 = find(Xs(:,4) == 0);

        % Find trials whose previous stimulus is on the opposite side
        % than the current block

        block_side = zeros(Nt,1);
        block_side(Xs(:,3) > 0.5) = -1;
        block_side(Xs(:,3) < 0.5) = 1;

        Xshift = circshift(Xs,1);
        idx_rev = find((sign(Xshift(:,4)) == side*block_side) & (block_side ~= 0));
        % Xshift2 = circshift(Xs,2);
        % idx_rev = find((sign(Xshift2(:,4)) == -block_side) & (sign(Xshift(:,4)) == -block_side) & (block_side ~= 0));

        % Combine the two constraints
        idx = intersect(idx0,idx_rev);

        if ~isempty(p_model)
            p_model_s = p_model(X(:,2) == iSession,:);
        end
        
        for j = 1:numel(idx)
            block_pos = trial_block(idx(j));

            tot_resp(block_pos) = tot_resp(block_pos) + 1;
            if sign(Xs(idx(j),6)) == block_side(idx(j))
                resp_match(block_pos) = resp_match(block_pos) + 1;
            end
            
            if ~isempty(p_model)
                model_match(block_pos) = model_match(block_pos) ...
                    + (block_side(idx(j)) == -1)*p_model_s(idx(j)) ...
                    + (block_side(idx(j)) == 1)*(1-p_model_s(idx(j)));
            end
            
            
        end
        
        
        
    end
end

if plotflag
    fontsize = 18;
    axesfontsize = 14;
    
    
    xx = 1:size(resp_match,2);
    yy = nansum(resp_match,1) ./ nansum(tot_resp,1);    
    idx = isfinite(yy);
    xx = xx(idx);
    yy = yy(idx);
    
    yy_serr = sqrt(yy .* (1 - yy)) ./ sqrt(nansum(tot_resp(:,idx),1));
    
    yyerr_down = yy - 1.96*yy_serr;
    yyerr_up = yy + 1.96*yy_serr;
    xxerr = [xx, fliplr(xx)];
    
    yyerr = [yyerr_down, fliplr(yyerr_up)];
    fill(xxerr, yyerr,'k','FaceAlpha',0.5,'LineStyle','none'); hold on;
    
    h(1) = plot(xx,yy,'k-','LineWidth',2);
    
    xx = 1:size(resp_match,2);
    yy = nansum(model_match,1) ./ nansum(tot_resp,1);
    idx = isfinite(yy);
    xx = xx(idx);
    yy = yy(idx);
    
    yy_serr = sqrt(yy .* (1 - yy)) ./ sqrt(nansum(tot_resp(:,idx),1));

    yyerr_down = yy - 1.96*yy_serr;
    yyerr_up = yy + 1.96*yy_serr;
    xxerr = [xx, fliplr(xx)];
    
    yyerr = [yyerr_down, fliplr(yyerr_up)];
    fill(xxerr, yyerr,'b','FaceAlpha',0.5,'LineStyle','none'); hold on;
    h(2) = plot(xx,yy,'b-','LineWidth',2);
    
    set(gcf,'Color','w');
    set(gca,'TickDir','out','FontSize',axesfontsize);
    box off;
    xlabel('Trials from block start','FontSize',fontsize);
    ylabel('Pr(choice matches block)','FontSize',fontsize);
    
    if side == -1
        title('Zero-contrast trials, after opposite-side trial','FontSize',fontsize);
    else
        title('Zero-contrast trials, after same-side trial','FontSize',fontsize);
    end
    
    plot([0,100],0.8*[1 1],'--k','LineWidth',1);
    
    ylim([0 1]);
    set(gca,'YTick',[0,0.8,1]);
    if params.runlength_min == 1
        hl = legend(h,'Data','Model (no minimum block length)');        
    else
        hl = legend(h,'Data','Model');
    end
    set(hl,'Location','NorthWest','Box','off','FontSize',fontsize);
    
end
















end