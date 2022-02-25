%----------------------------------------------------------------
% 
%                       %%%%%%%%%%
%                       % prodBO %
%                       %%%%%%%%%%
% 
% This function calculates the product of the probability to monitor a sequence O with the probability of the backward B distribution.
% 	function prob=prodBO(B,O,nc)
% 
% Entry:	B{Np}(N,M) is the distribution probability symbol matrix  for each parameter.
% B(i,k) is the probability to obtain the kth  symbol when we are at the state i.
% 
% O is a matrix with the sequence. 
% O{Np}{number of vectors, number of the symbols}.
% 
% nc is the sequence from O to estimate.
% 
% Result:	prob is the result of the product. That means prob is the probability for every state to generate the extracted component nc from the sequence O.
%-----------------------------------------------------------------
function prob=prodBO(B,O,nc)

Np=size(B,1);
prob=ones(size(B{1},1),1);
for ip=1:Np
   prob=prob.*(B{ip}*O{ip}(nc,:)');
end
prob=prob+realmin;
return
