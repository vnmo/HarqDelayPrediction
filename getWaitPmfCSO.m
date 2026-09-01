function [waitFunction, pDrop] = getWaitPmfCSO(varargin)
% Input: arguments accepted by getWaitFunctionCSO.
% Outputs: waiting-delay function waitFunction and auxiliary overflow estimate pDrop.
[waitFunction, pDrop] = getWaitFunctionCSO(varargin{:});
end
