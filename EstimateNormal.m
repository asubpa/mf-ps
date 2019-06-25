function [nmap] = EstimateNormal(L,I,SH_THRESHOLD, mask)
%
%   Normal map estimation using a microfacet-based model
%   For review
%   Input:
%   L: [n_scene x 3] directions of light.
%   I: [n_scene x n_pixel] image intensities.
%   SH_THRESHOLD: [scalar] Cast shadow threshold.
%   mask: [n_pixel] binary mask.
%   Output:
%   nmap : [3 x n_pixel] normal vectors.
%   ======================================

    if nargin < 3
       SH_THRESHOLD = 0;
       mask = true(size(I,2));
    end

    H = [L(:,1) L(:,2) L(:,3)+1];
    H = H./sqrt(sum(H.*H,2));
    
    npix = size(I,2);
    nmap = zeros(3,npix);
    
    for ipix = 1:npix
        if(mask(ipix))
            profile = I(:,ipix);
            iSelected = profile > 0 & profile > median(profile)*SH_THRESHOLD;
            perSampleH = H(iSelected,:);
            perSampleL = L(iSelected,:);
            profile = profile(iSelected);
            [n,~] = Estimate_Normal_CPS_MicroBRDF_General_Trigo(perSampleH,perSampleL,profile);
            nmap(:,ipix) = n;
        else
            nmap(:,ipix) = [0,0,0];
        end
    end

end