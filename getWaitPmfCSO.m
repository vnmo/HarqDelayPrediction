function [waitFunction, pDrop] = getWaitPmfCSO(varargin)
% Input: arguments accepted by getWaitFunctionCSO.
% Outputs: waiting-delay function waitFunction and arrival-discard probability pDrop for q=Qmax.
[waitFunction, pDrop] = getWaitFunctionCSO(varargin{:});
end
