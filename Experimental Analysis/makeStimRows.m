function S= makeStimRows(Stim, n)
%  S = makeStimRows(Stim, n, flag); not used in paper
%
%  Converts raw stimulus to a matrix where each row is loaded with the full
%  space-time stimulus at a particular moment in time.  The resulting
%  matrix has length equal to the number of time points and width equal to
%  the (number of spatial dimensions) x (kernel size n).
%
%  Inputs: 
%   Stim = stimulus, first dimension is time, other dimensions are spatial
%          dimensions
%   n = size of temporal kernel; number of time samples to include in each
%       row of stimulus matrix.
%
%  Output: S = matrix where each row is the size of the linear kernel
%
%  Last updated:  6/30/2005, J. Pillow
% Made memory efficient

S = hankel(Stim(1:end-n+1), Stim(end-n+1:end));
