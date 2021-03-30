function x = ChebyshevRoots( n, type, range )
%   USAGE: x = ChebyshevRoots( n [,type [, range]] )
%
% This method returns the values of the roots of the Chebyshev polynomial, 
% either type one or type two of degree n.  In the literature, these are 
% often referred to as T_n(x) for type 1 or U_n(x) for type 2.  If the
% optional parameter type is omitted, it is assumed to be type 1.  This
% function also supports scaling and translating the roots of the
% polynomial to lie within a range specified.  The polynomials T_n(x) and
% U_n(x) have n roots which lie within [1,-1].  It is often useful to scale
% and translate these roots to have the same relative distance but lie over
% a different range so that they may be used as the nodes for interpolation
% of a function which does not lie within this [-1, 1] range.
%
%   PARAMETERS:
%       n - The degree of the Chebyshev polynomial whose roots we wish to
%       find.
%
%       type [optional] - Either 'Tn' for type 1 or 'Un' for type 2
%       depending on whether you wish to generate the roots of the type 1
%       or type 2 polynomial.  'Tn' is the default if this parameter is 
%       omitted.
%
%       range [optional] - The Chebyshev polynomial is defined over [-1, 1]
%       this parameter allows the roots of the polynomial to be translated
%       to be within the range specified.  ie. The relative distance of the
%       Chebyshev nodes to each other will be the same but their values
%       will span the provided range instead of the range [-1, 1].
%
%   RETURNS:
%       A vector with the n roots of the Chebyshev polynomial T_n(x) or
%       U_n(x) of degree n, optionally scaled to lie within the range 
%       specified.
%
%   AUTHOR: 
%       Russell Francis <rf358197@ohio.edu>
%
%   THANKS:
%       John D'Errico - Provided reference to the Abramowitz and Stegun book 
%       which provides a rigorous definition of the two types of Chebyshev
%       polynomials and is suggested reading particularly chapter 22 for
%       interested parties.
%
%   REFERENCES:
%
% [1] Abramowitz, M.; Stegun, I.A. Handbook of Mathematical Functions with 
% Formulas, Graphs, and Mathematical Tables. U.S. Department of Commerce.
% Online version available at:
% http://www.knovel.com/knovel2/Toc.jsp?BookID=528&VerticalID=0
%
% [2] http://en.wikipedia.org/wiki/Chebyshev_polynomials
%
% [3] Burden, Richard L.; Faires, J. Douglas Numerical Analysis 8th ed. 
% pages 503-512.
%

%%
% Verify our parameters.
%%

% The degree must be specified and must be greater than or equal to 1.
if( nargin() < 1 )
    error( 'You must provide the parameter n which is the degree of the polynomial you wish to calculate the roots of.' );
else
    if( n < 1 )
        error( 'The parameter n must be greater than or equal to 1.' );
    end
end

% The type of the Chebyshev polynomial to calculate the roots of, optional
% and defaults to T_n(x)
if( nargin() < 2)
    type = 'Tn';
else
    if( (strcmp( type, 'Un' ) ~= 1 ) && (strcmp( type, 'Tn' ) ~= 1) )
        error( 'The type parameter which was specified is not valid!, Please specify either "Tn" or "Un"' );
    end        
end

% The range which we wish to scale and translate the result to, optional
% and the default is to not scale or translate the result.
if( nargin() < 3 )
    range = [-1 1];
else
    if( length( range ) ~= 2 )
        error( 'The parameter range must contain two values.' );
    end

    if( range(1) == range(2) )
        error( 'The parameter range must contain two distinct values.' );
    end
        
    range = sort(range);    
end

%%
% Begin to compute our Chebyshev node values.
%%
if( n == 1 )
    x = [0];
else
    if( strcmp( type, 'Tn' ) == 1 )
        x = [(pi/(2*n)):(pi/n):pi];
    else
        x = [(pi/(n+1)):(pi/(n+1)):((n*pi)/(n+1))];
    end
    x = sort( cos(x) );
end

%%
% x now contains the roots of the nth degree Chebyshev polynomial,
% we need to scale and translate the result if necessary.
%%
if( (range(1) ~= -1) || (range(2) ~= 1) )
    M = eye(n+1);    
    % Calculate the scaling factor to apply to the nodes.
    sf = (range(2) - range(1)) / 2;
    % Calculate the translation to apply to the nodes.
    tl = range(1) + sf;
    % Generate our transformation matrix.
    M(1:n,1:n) = M(1:n,1:n) * sf;
    M(n+1,1:n) = tl;
    % Apply to our earlier result.
    x = [x 1] * M;
    x = x(1:n);
end
return; % The x values of the Chebyshev nodes.
