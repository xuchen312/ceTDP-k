function varargout = spm_ceTDP_minsize(xSPM,u,alpha,guess)
%
% Find the minimum significant cluster size threshold for a given alpha
%
% USAGE:
%   [k, Pc] = spm_ceTDP_minsize([xSPM, u, alpha, guess])
%
% INPUTS (optional):
%   xSPM  - Structure containing SPM, distribution & filtering details
%   u     - Cluster defining threshold (if < 1, treated as uncorrected P-value)
%   alpha - FWE-corrected level (default: 0.05)
%   guess - Search method or initial guess:
%           NaN (default) - Newton-Raphson search
%           scalar - Report P-value for that cluster size
%           vector - Brute-force search over cluster sizes
%
% OUTPUTS (optional):
%   k     - Minimum significant cluster size (in voxels)
%   Pc    - Corrected P-value for cluster size k
%
% BASED ON:
%   $Id: CorrClusTh.m,v 1.3 2012/10/29 09:32:10 nichols Exp $
%   Thomas Nichols, Marko Wilke (SPM2 CorrClusTh.m,v 1.10 2006/05/08 15:31:02)

% -------------------- settings --------------------
epsP   = 1e-6;
stepDu = 1e-6;
maxi   = 100;

k  = NaN;
Pc = NaN;

% -------------------- input check --------------------
if nargin < 1 || isempty(xSPM)
    [~, xSPM] = spm_getSPM;
end

df   = xSPM.df;
STAT = xSPM.STAT;
n    = 1;
R    = xSPM.R;
S    = xSPM.S;
M    = xSPM.M;

VOX  = sqrt(diag(M(1:3,1:3)' * M(1:3,1:3)))';
FWHM = xSPM.FWHM;
FWHMmm = FWHM .* VOX;

v2r  = 1 / prod(FWHM(~isinf(FWHM)));

sf_ShowVolInfo(R,S,VOX,FWHM,FWHMmm);

% -------------------- threshold u --------------------
if nargin < 2 || isempty(u)
    u = xSPM.u;
end

if u <= 1
    uP = u;
    u  = spm_u(u,df,STAT);
else
    switch STAT
        case 'T'
            uP = 1 - spm_Tcdf(u,df(end));
        case 'F'
            uP = 1 - spm_Fcdf(u,df);
        case 'X'
            uP = 1 - spm_Xcdf(u,df(1));
        otherwise
            error('Unsupported STAT type');
    end
end

% -------------------- alpha --------------------
if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end

% -------------------- guess --------------------
if nargin < 4 || isempty(guess)
    guess = NaN;
end

epsP = alpha * epsP;
Status = 'OK';

% =========================================================
% CASE 1: direct evaluation
% =========================================================
if isscalar(guess) && ~isnan(guess)

    k  = guess;
    Pc = spm_P(1,k*v2r,u,df,STAT,R,n,S);
    Status = 'JustPvalue';

% =========================================================
% CASE 2: trivial 1-voxel significance
% =========================================================
elseif spm_P(1,1*v2r,u,df,STAT,R,n,S) < alpha

    k  = 1;
    Pc = spm_P(1,1*v2r,u,df,STAT,R,n,S);
    Status = 'TooRough';

% =========================================================
% CASE 3: automatic Newton search
% =========================================================
elseif isnan(guess)

    [~,~,~,En] = spm_P(1,0,u,df,STAT,R,n,S);
    kr  = En;
    rad = kr^(1/3);

    % ---- coarse bounds ----
    Pcl = 1;
    Pcu = 0;
    radl = rad;
    radu = rad;

    while Pcl > alpha
        Pcu  = Pcl;
        radl = radu;
        radu = radu * 1.1;
        Pcl  = spm_P(1,radu^3,u,df,STAT,R,n,S);
    end

    % ---- interpolation ----
    if abs(Pcu - Pcl) < 1e-12
        rad = radu;
    else
        rad = radl*(alpha-Pcl)/(Pcu-Pcl) + ...
              radu*(Pcu-alpha)/(Pcu-Pcl);
    end

    % ---- Newton refinement ----
    d = 1;
    iter = 1;

    while abs(d) > epsP

        Pc  = spm_P(1,rad^3,u,df,STAT,R,n,S);
        Pc1 = spm_P(1,(rad+stepDu)^3,u,df,STAT,R,n,S);
        Pc0 = spm_P(1,(rad-stepDu)^3,u,df,STAT,R,n,S);

        dP = (Pc1 - Pc0) / (2*stepDu);

        if abs(dP) < epsP
            break;
        end

        d = (alpha - Pc) / dP;

        d = max(min(d,(radu-radl)/2),-(radu-radl)/2);

        rad = rad + d;
        rad = max(min(rad,radu),radl);

        iter = iter + 1;
        if iter >= maxi
            Status = 'TooManyIter';
            break;
        end
    end

    kr = rad^3;
    k  = ceil(kr / v2r);
    Pc = spm_P(1,k*v2r,u,df,STAT,R,n,S);

% =========================================================
% CASE 4: brute force
% =========================================================
else

    guess = sort(guess(:));

    for i = 1:numel(guess)

        k  = guess(i);
        Pc = spm_P(1,k*v2r,u,df,STAT,R,n,S);

        if Pc <= alpha
            break;
        end
    end

    if Pc > alpha
        Status = 'OutOfRange';
    end
end

% -------------------- output --------------------
print_output(Status,uP,u,alpha,k,Pc,guess);

if nargout==1
    varargout = {k};
elseif nargout>1
    varargout = {k,Pc};
end

end

% =========================================================
function print_output(Status,uP,u,alpha,k,Pc,guess)

switch Status

 case 'JustPvalue'
  fprintf(['\n====================================================\n' ...
           ' Cluster extent evaluation (P-value only)\n' ...
           '====================================================\n\n' ...
           ' Threshold:\n' ...
           '   P = %0.3f  (T = %0.4f)\n\n' ...
           ' Result:\n' ...
           '   k = %d voxels\n' ...
           '   corrected P = %g\n\n'], ...
           uP,u,k,Pc);

 case 'OK'
  fprintf(['\n====================================================\n' ...
           ' Cluster extent threshold (FWE corrected)\n' ...
           '====================================================\n\n' ...
           ' Threshold:\n' ...
           '   P = %0.3f  (T = %0.4f)\n' ...
           '   alpha = %0.3f\n\n' ...
           ' Result:\n' ...
           '   k = %d voxels\n' ...
           '   corrected P = %0.4f\n\n'], ...
           uP,u,alpha,k,Pc);

 case 'TooRough'
  fprintf(['\n====================================================\n' ...
           ' Warning: single voxel already significant\n' ...
           '====================================================\n\n' ...
           ' Threshold:\n' ...
           '   P = %0.3f  (T = %0.4f)\n\n' ...
           ' Result:\n' ...
           '   k = 1\n' ...
           '   corrected P = %g\n\n'], ...
           uP,u,Pc);

 case 'TooManyIter'
  fprintf(['\n====================================================\n' ...
           ' Warning: no convergence\n' ...
           '====================================================\n\n' ...
           ' Newton search failed to converge.\n\n']);

 case 'OutOfRange'
  fprintf(['\n====================================================\n' ...
           ' Warning: range insufficient\n' ...
           '====================================================\n\n' ...
           ' Search range: %g - %g voxels\n' ...
           ' Smallest P found: %g\n\n'], ...
           guess(1),guess(end),Pc);

 otherwise
  error('Unknown status');

end

end

% =========================================================
function sf_ShowVolInfo(R,S,VOX,FWHM,FWHMmm)

fprintf('\n====================================================\n');
fprintf(' Search Volume Summary\n');
fprintf('====================================================\n\n');

fprintf(' Voxels:\n');
fprintf('   N = %d\n', S);
fprintf('   size = [%.2f %.2f %.2f] mm\n\n', VOX);

fprintf(' Physical space:\n');
fprintf('   volume (mm^3) = %.2f\n', S * prod(VOX));
fprintf('   volume (L)    = %.4f\n', S * prod(VOX) / 100^3);
fprintf('   RESELs        = %.4f\n\n', R(end));

fprintf(' Smoothness:\n');
fprintf('   FWHM mm    = [%.2f %.2f %.2f]\n', FWHMmm);
fprintf('   FWHM vox   = [%.2f %.2f %.2f]\n\n', FWHM);

end
