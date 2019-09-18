function [meanValue] = get_training_mean(data, Tr, Td, Gr, Gd)
%GET_TRAINING_MEAN Returns the mean value over all training cells
% 
% Arguments:
% data: matrix of size (2*Tr + 2*Gr+1)x(2*Td + 2*Gd+1) from which we want to extract the mean from the training cells.
% Tr: number of training cells (on each side) in the range direction.
% Td: number of training cells (on each side) in the doppler direction 
% Gr: number of guard cells (on each side) in the range direction.
% Gd: number of guard cells (on each side) in the doppler direction 

data = db2pow(data); % the matrix has values in DB, hence we get to convert it
filter = ones(2*Tr + 2*Gr + 1, 2*Td + 2*Gd + 1); % we create a filter with the same size as data
filter(Tr+1:(Tr+2*Gr+1), Td+1:(Td+2*Gd+1)) = 0; % and set to zero the values we don't want to: guard cells and cut
meanValue = pow2db(mean(filter.*data,'all')); % We get the mean and convert to DB the result
end