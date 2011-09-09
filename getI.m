% 
% IEEE 738-2006 MATLAB implementation
%
% IEEE Standard for Calculating the Current-Temperature of Bare Overhead Conductors
% 
% Copyright (c) 2011 Steven Blair
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%

% returns the current required for a conductor temperature, Tc, and other parameters
function [I] = getI(Tc, R_T_high, R_T_low, T_high, T_low, Ta, rho_f, D, epsilon, alpha, Q_se, theta, area)
    R_Tc = ((R_T_high - R_T_low) / (T_high - T_low)) * (Tc - T_low) + R_T_low;      % replace with calculation
    
    % Conductor heat loss from natural convection, assuming no wind
    q_cn = 0.0205 * rho_f^0.5 * D^0.75 * (Tc - Ta)^1.25;
    
    % Conductor heat loss from radiation
    q_r = 0.0178 * D * epsilon * (((Tc + 273)/100)^4 - ((Ta + 273)/100)^4);
    
    % Solar heat gain in conductor
    q_s = alpha * Q_se * sin(theta) * area;
    
    I=sqrt((q_cn + q_r - q_s) / R_Tc);
end
