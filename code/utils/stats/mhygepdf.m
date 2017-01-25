function x = mhygepdf(m,n)
%MHYGEPDF Multivariate hypergeometric probability density function.
% X = MHYGEPDF(M,N) returns the multivariate hypergeometric probability
% density function at M with integer parameters in N. Note: The density
% function is zero unless the elements in M are integers.
%
% If there are m_i elements of class i in a population and you take n 
% elements at random without replacement, then the number of elements of
% each class in the sample (x_1,x_2,...,x_c) has the multivariate 
% hypergeometric distribution. This has the same relationship to the 
% multinomial distribution that the hypergeometric distribution has to the
% binomial distribution--the multinomial distribution is the 
% 'with-replacement' distribution and the multivariate hypergeometric is
% the 'without-replacement' distribution.
%
% Syntax: function x = mhygepdf(m,n) 
%      
% Inputs:
% m - vector of total (population) elements of class i  
% n - vector of sampled elements of class i  
%
% Output:
% x - multivariate hypergeometric probability for values x_1,x_2,...,x_c
%
% Example. From the example given on the web Wikipedia: The Free 
%  Encyclopedia [http://en.wikipedia.org/wiki/Hypergeometric_distribution].
%  Suppose there are 5 black, 10 white, and 15 red marbles in an urn. You
%  reach in and randomly select six marbles without replacement. What is 
%  the probability that you pick exactly two of each color?
%
%             m = [5,10,15]; n = [2,2,2];
%
% Calling on Matlab the function: 
%             x = mhygepdf(m,n)
%
% Answer is: (in format long)
%
% x =
%
%   0.07957559681698
%
% Created by A. Trujillo-Ortiz, R. Hernandez-Walls, F.A. Trujillo-Perez
%            and N. Castro-Castro
%            Facultad de Ciencias Marinas
%            Universidad Autonoma de Baja California
%            Apdo. Postal 453
%            Ensenada, Baja California
%            Mexico.
%            atrujo@uabc.mx
% Copyright. September 24, 2009.
%
% ---Con cariño para Ney.----
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A., R. Hernandez-Walls, F.A. Trujillo-Perez and 
%   N. Castro-Castro (2009). mhygepdf:Multivariate hypergeometric 
%   probability density function. A MATLAB file. [WWW document].
%   URL http://www.mathworks.com/matlabcentral/fileexchange/25394
%

if nargin < 2, 
    error('mhygepdf:TooFewInputs','Requires two input arguments.'); 
end

[errorcode m n] = distchck(2,m,n);

if ~isvector(m) || ~isvector(n)
    error('mhygepdf:InvalidData',...
        'MHYGEPDF requires vector rather than matrix data.');
elseif numel(m) ~= numel(n)
    error('mhygepdf:InputSizeMismatch',...
'MHYGEPDF requires the data vectors to have the same number of elements.');
end

if errorcode > 0
    error('mhygepdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

m = m(:);
n = n(:);

if any(fix(m) ~= m) || any(m < -1) || ~isa(m,'double') || ~isreal(m)
    error('MATLAB:mhygepdf:NNegativeInt', ...
        'M must be a vector of non-negative and integers.')
end

if any(fix(n) ~= n) || any(n < -1) || ~isa(n,'double') || ~isreal(n)
    error('MATLAB:mhygepdf:NNegativeInt', ...
        'N must be a vector of non-negative and integers.')
end

M = sum(m);
N = sum(n);
o = numel(m);

A = [];
for i=1:o,
    a = exp(gammaln(m(i) + 1) - gammaln(n(i) +1) - gammaln(m(i) - n(i) + 1));
    A = [A a];
end

B = exp(gammaln(M + 1) - gammaln(N +1) - gammaln(M - N + 1));

x = max(cumprod(A))/B;

return,