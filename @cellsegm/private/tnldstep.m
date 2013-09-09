function [r] = tnldstep(L,a,b,c)

% function for computing div(D*grad(u))
% D = [a b  
%      c d]

Lpc = transim( L, 1, 0, 0 );
        Lpp = transim( L, 1, 1, 0 );
        Lcp = transim( L, 0, 1, 0 );
        Lnp = transim( L, -1, 1, 0 );
        Lnc = transim( L, -1, 0, 0 );
        Lnn = transim( L, -1, -1, 0 );
        Lcn = transim( L, 0, -1, 0 );
        Lpn = transim( L, 1, -1, 0 );
        
        anc = transim( a, -1, 0, 0 );
        apc = transim( a, +1, 0, 0 );
        bnc = transim( b, -1, 0, 0 );
        bcn = transim( b, 0, -1, 0 );
        bpc = transim( b, +1, 0, 0 );
        bcp = transim( b, 0, +1, 0 );
        ccp = transim( c, 0, +1, 0 );
        ccn = transim( c, 0, -1, 0 );
        
        r = -1/4 * (bnc+bcp) .* Lnp + ...
             1/2 * (ccp+c)   .* Lcp + ...
             1/4 * (bpc+bcp) .* Lpp + ...
             1/2 * (anc+a)   .* Lnc - ...
             1/2 * (anc+2*a+apc+ccn+2*c+ccp) .* L + ...
             1/2 * (apc+a)   .* Lpc + ...
             1/4 * (bnc+bcn) .* Lnn + ...
             1/2 * (ccn+c)   .* Lcn - ...
             1/4 * (bpc+bcn) .* Lpn;


