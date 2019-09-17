function [meanValue] = get_training_mean(data, Tr, Td, Gr, Gd)
%GET_TRAINING_MEAN Summary of this function goes here
% Returns the mean value over all training cells
data = db2pow(data);
filter = ones(2*Tr + 2*Gr + 1, 2*Td + 2*Gd + 1);
filter(Tr+1:(Tr+2*Gr+1), Td+1:(Td+2*Gd+1)) = 0;
meanValue = pow2db(mean(filter.*data,'all'));
end

