function [betaMinNormBest, lambda,residualNorm,solutionNorm] = sourceLoc(template , data, plotLcurve)
% Source localisation using lcurve regularisation
% INPUT:
% template = templates (electrodes x ROI)
% data = data to fit (electrodes x time or electrodes x sine and cosine amplitudes)
% optional plotLcurve: default 0, 1 for plotting L-curve
% ATTENTION: Make sure that template & data have the same montage and the same reference   
% use createCustomTemplates to fit your EEG montage if needed
% OUTPUT:
% betaMinNormBest = activity in each ROI for the best regularisation value
% other optional output: lambda = regularisation term, 
% residuals & solutions from the minimum-norm estimates 
% USAGE: roiActivity = minNormTemplate_lcurve(template, data)

if (nargin==2), plotLcurve =0; end % do not plot Lcurve

[u,s,v] = csvd(template);
[lambda,residualNorm,solutionNorm,lambdaGridMinNorm] = l_curve_modified(u,s,data,'Tikh',plotLcurve);

if lambda > lambdaGridMinNorm(4) || lambda < lambdaGridMinNorm(end-3) 
    warning(['Risk of overfitting. You might want to test another (higher) regularisation value.' ... 
        'Try plotting L-curve to determine the corner of the curve using minNormFast_lcurve( template , data, 1)'...
        'Then use that corner value in minNormFix'])
end
    
betaMinNormBest = zeros([size(template,2) size(data, 2)]);
for ll = 1:size(data, 2)
    betaMinNormBest(:, ll) = tikhonov(u, s, v, data(:, ll), lambda);
end
