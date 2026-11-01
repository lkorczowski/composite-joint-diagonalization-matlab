function [C,Ct] = convertDat_BAJD(X,Cf)
% [C,Ct] = convertDat_BAJD(X,Cf)
% Build the set of symmetric target matrices for AJD (gp_AJD_GPT) from
% bilinear data X and additional symmetric matrices Cf.
%
% INPUTS
% -----------------
%   X  : observation matrices [N x T x K] (e.g. trials)
%   Cf : symmetric matrices [N x N x F] (e.g. co-spectra)
%
% OUTPUTS
% -----------------
%   C  : [N x N x (F+K)] concatenation of Cf, rescaled so that its total
%        trace equals the total trace of Ct, followed by Ct
%   Ct : [N x N x K] covariance matrix of each X(:,:,k)
%
% see also: gp_AJD_GPT, simBAJD_dat
%
% *** History: 09-Fev-2016
% *** Author: Florent BOUCHARD & Louis KORCZOWSKI, GIPSA-Lab, 2016
% *** Reference: MINING THE BILINEAR STRUCTURE OF DATA WITH APPROXIMATE
%                JOINT DIAGONALIZATION
%                L. Korczowski, F. Bouchard, C. Jutten, M. Congedo
% *** Contact: louis.korczowski@gmail.com
% *** Licence: GNU GPLv3

N = size(X,1);
T = size(X,2);
K = size(X,3);

% get Ct
Ct = zeros(N,N,K);
trat=0;
for k=1:K
    Ct(:,:,k) = cov(X(:,:,k)');
    trat=trat+trace(Ct(:,:,k));
end
traf=0;
for f=1:size(Cf,3)
   traf=traf+trace(Cf(:,:,f));
end
% Ct = mean(Ct,3);

C = cat(3,trat*Cf/traf,Ct);