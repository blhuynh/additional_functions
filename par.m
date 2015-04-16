function [ req ] = par( r1,r2 )
%PAR Returns equivalent impedance of parallel components.
%   Req = PAR(r1,r2) returns a single equivalent real impedance of two
%   parallel resistive components.

req = r1*r2 / (r1+r2);

end

