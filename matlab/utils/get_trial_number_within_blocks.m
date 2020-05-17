function [trial_block,trials2end_block] = get_trial_number_within_blocks(X)
%GET_TRIAL_NUMBER_WITHIN_BLOCKS Given a session

Nt = size(X,1);                 % tot # trials in the session    
trial_block = zeros(Nt,1);      % trial number inside block
trials2end_block = zeros(Nt,1); % trial to the end of block
block_starts = find(diff([0;X(:,3);0]))';
for iBlock = 1:numel(block_starts)-1
    idx = block_starts(iBlock):(block_starts(iBlock+1)-1);
    trial_block(idx) = 0:numel(idx)-1;
    trials2end_block(idx) = numel(idx):-1:1;
end

end