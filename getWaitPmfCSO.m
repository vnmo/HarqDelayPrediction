function [waitFunction, pDrop] = getWaitPmfCSO(varargin)
% Input: arguments accepted by getWaitFunctionCSO.
% Legacy entry point: the first output is the normalized aggregate
% waiting-delay function. pDrop is the q=Qmax arrival-rejection mass.
[waitFunction, pDrop] = getWaitFunctionCSO(varargin{:});
end
