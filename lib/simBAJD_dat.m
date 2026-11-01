function [X,C,A,E,info] = simBAJD_dat(options)

% [X,C,A,E,info] = simBAJD_dat(options)
% Generate random data for (composite) bilinear AJD: trials following the
% bilinear model and symmetric matrices following the linear model, both
% sharing the same spatial mixing matrix A.
%   X(:,:,k) = A*S_k*E' + noise     (S_k diagonal)
%   C(:,:,f) = A*Cs_f*A' + noise    (Cs_f diagonal, positive)
%
% INPUTS
% -----------------
%   options* is a structure with
%     .N is the number of sources/sensors. default: 8
%     .T is the number of samples. default: 128
%     .K is the number of matrices X. default: 20
%     .F is the number of matrices C. default: 20
%     .condA is 1 for an orthogonal A, or a range [min max] of condition
%     number for a random A. default: 1
%     .condE is 1 for an orthogonal E, or a range [min max] of condition
%     number for a random E (needs mvnrnd, Statistics Toolbox). default: 1
%     .SNR is the signal to noise ratio. default: inf (no noise)
%
% OUTPUTS
% -----------------
%   X : the observation matrices [N x T x K]
%   C : the symmetric matrices [N x N x F]
%   A : the true spatial mixing matrix [N x N]
%   E : the true temporal mixing matrix [T x N]
%   info : structure with the sources S, Cs, the noise No, Noc and options
%
% see also: convertDat_BAJD, gp_CAJD_GPT
%
% *** Author: Florent BOUCHARD, GIPSA-Lab, 2015
% initialize parameters according to options input
%%% default
opt = struct('N',8,'T',128,'K',20,'F',20,'condA',1,'condE',1,'SNR',inf);
%%% get final options
if nargin<1
   options=struct; 
end
options = gp_getOpt(options,opt);
%%% get parameters
N     = options.N;
T     = options.T;
K     = options.K;
F     = options.F;
condA = options.condA;
condE = options.condE;
SNR   = options.SNR;

% mixing matrices
A = aux_spatMix(N,condA);
E = aux_tempMix(N,T,condE);

% generate X
%%% generate random sources
S = zeros(N,N,K);
for k=1:K
    for i=1:N
        while S(i,i,k)<1e-4
            S(i,i,k) = randn^2;
        end
    end
end
%%% generate mixture
X  = zeros(N,T,K);
No = zeros(N,T,K);
for k=1:K
    % mixture without noise
    X(:,:,k)  = A*S(:,:,k)*E';
    % generate noise
    No(:,:,k) = randn(N,T)/sqrt(SNR);
    % get final mixture
    X(:,:,k)  = X(:,:,k) + No(:,:,k);
end

% generate C
%%% generate cospectra sources
Cs = zeros(N,N,F);
for f=1:F
    for i=1:N
        while Cs(i,i,f)<1e-4
            Cs(i,i,f) = randn^2;
        end
    end
end
%%% generate mixture
C   = zeros(N,N,F);
Noc = zeros(N,N,F);
for f=1:F
    C(:,:,f) = A*Cs(:,:,f)*A';
    % generate noise
    tmp = randn(N,N);
    Noc(:,:,f) = .5*(tmp+tmp')/SNR;
    % get final mixture
    C(:,:,f) = C(:,:,f) + Noc(:,:,f);
end

% fill info struct
info = struct('S',S,'No',No,'Cs',Cs,'Noc',Noc,'options',options);

end

function A = aux_spatMix(N,condA)
    % generation of random spatial mixing matrices according to a condition
    % number range
    if condA==1
        A = orth(randn(N));
    else
        A = randn(N);
        while cond(A)<condA(1) || cond(A)>condA(2)
            A = randn(N);
        end
    end
end

function E = aux_tempMix(N,T,condE)
    % generation of random temporal mixing matrices according to a 
    % condition number range
    % warning: if condE not 1, then it should be a large range in order to 
    %          find something good
    if condE ==1
        E = orth(randn(T,N));
    else
        sig = aux_gensymposmat(N);
        while sqrt(cond(sig))<condE(1) || sqrt(cond(sig))>condE(2)
            sig = aux_gensymposmat(N);
        end
        E = mvnrnd(zeros(1,N),sig,T);
    end

end

function sig = aux_gensymposmat(N)
    % generate random SPD matrices, very random condition number, often bad
    % never very good, if you want it to be very good (around 1.5), use
    % D = diag(1+rand(N,1));
    D = diag(randn(N, 1).^2);
    [Q, ~] = qr(randn(N));
    sig = Q*D*Q';
end