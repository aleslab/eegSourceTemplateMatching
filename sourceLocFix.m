function [betaMinNormBest, rho, eta] = sourceLocFix( template , data, nLambda )
% Source localisation using a pre-set regularisation parameter (nLambda)
% INPUT: 
% template = templates (electrodes x ROI)
% data = data to fit (electrodes x time or electrodes x sine and cosine amplitudes)
% ATTENTION: Make sure that template & data have the same montage and the same reference   
% use createCustomTemplates to fit your EEG montage if needed
% OUTPUT:
% betaMinNormBest = activity in each ROI for the best regularisation term
% other optional output: residuals (rho) & solutions (eta) from the
% USAGE: activityROI = minNormFast_lcurve(avMap, data, 6000)

[u,s,v] = csvd(template);
rho = zeros(size(data,2),1);
eta = zeros(size(data,2),1);
betaMinNormBest = zeros(size(template,2),size(data,2));
for tt = 1:size(data, 2)
    [betaMinNormBest(:, tt),rho(tt),eta(tt)] = tikhonov(u, s, v, data(:, tt), nLambda);
end
