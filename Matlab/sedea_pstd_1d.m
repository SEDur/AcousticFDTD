%
% Pseudo-sprectral time domain method object 1d
% 
%
% Made by S Durbridge
%
% Last Edited: 10/01/2017
%
% Next Task: Write method for calculating simulation parameters
%
% Terminology:
% solid: an N dimensional object that does not contain a portion of the
% computed domain within its area i.e. a silouette in the 2d case, whose
% boundaries are used but the space inside is of no consequence to the
% simulation
% 
% Rules:
% Anything tangible within computation domain such as solids and their
% boundaries are referenced in agrogate by address, and by  using 
% the left hand rule, starting from 0
%
% Individual boundary segments (single addresses) are referenced as such
%
% Base of domain is always a rectangle, and the user builds solids onto the
% base to define a custom shape.
%
% The spatial resolution dx directly defines the resolution of solids,
% but an obstacel canot be less than 2*dx
%

classdef sedea_pstd_1d
    properties
        %Simulation Properties
        fs %Acoustic sampling rate
        dx %Spatial sampling rate in x direction
        dy %Spatial sampling rate in y direction
        dt %Time step
        
        %Domain Properties
        Length %Length of domain
        uArraysize %Size of the velocity array
        pArraysize %Size of the pressure array
        boundaryImpedance %Coefficient for the PML at either end of the domain
        c %Speed of sound in uniform medium 
        
        %Source Properties
            %Source Locations
            %Source Signals

        numSources

        
    end
    
    methods
        function obj = sedea_pstd_1d(fs, xLength, yLength)
            obj.fs = fs;
            obj.xLength = xLength;
            obj.yLength = yLength;
        end
        function [coefs] = sedea_rbjM_lpf(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
            a0 = A + alpha;
            b0 = (1 - cos(w0)) / 2;
            b1 = 1 - cos(w0);
            b2 = (1 - cos(w0)) / 2;
            a1 = -2 * cos(w0);
            a2 = 1 - alpha;
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end
        
        function  [coefs] = sedea_rbjM_hpf(obj)
            
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
            a0 = A + alpha;
            b0 = (1 + cos(w0)) / 2;
            b1 = -(1 + cos(w0));
            b2 = (1 + cos(w0)) / 2;
            a1 = -2 * cos(w0);
            a2 = 1 - alpha;
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end
        
        function  [coefs] = sedea_rbjM_bpfcq(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
            a0 = A + alpha;
            b0 = sin(w0)/2;
            b1 = 0.0;
            b2 = -sin(w0)/2;
            a1 = -2 * cos(w0);
            a2 = A - alpha;
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end
        
        function  [coefs] =  sedea_rbjM_bpfcg(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
            a0 = 1.0 + alpha;
            b0 = alpha;
            b1 = 0;
            b2 = -alpha;
            a1 = -2.0*cos(w0);
            a2 =  1.0 - alpha;
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end
        
        function  [coefs] = sedea_rbjM_notch(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
            a0 = A + alpha;
            b0 = 1.0;
            b1 = -2.0*cos(w0);
            b2 = 1 + alpha;
            a1 = -2.0*cos(w0);
            a2 =  1.0 - alpha;
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end
        
        function  [coefs] =  sedea_rbjM_apf(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
            b0 = 1.0 - alpha;
            b1 = -2.0*cos(w0);
            b2 = 1 + alpha;
            a0 = A + alpha;
            a1 = -2.0*cos(w0);
            a2 =  1.0 - alpha;
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end
        
        function  [coefs] = sedea_rbjM_pek(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
            b0 = 1.0 + alpha * A;
            b1 = -2.0*cos(w0);
            b2 = 1 - alpha * A;
            a0 = 1 + alpha * A;
            a1 = -2.0*cos(w0);
            a2 =  1.0 - alpha / A;
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end
        
        function  [coefs] = sedea_rbjM_ls(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
            b0 =      A * ((A+1) - (A-1)*cos(w0) + 2 * sqrt(A) * alpha);
            b1 =  2 * A * ((A-1) - (A+1)*cos(w0));
            b2 =      A * ((A+1) - (A-1)*cos(w0) - 2 * sqrt(A) * alpha);
            a0 =          ((A+1) + (A-1)*cos(w0) + 2 * sqrt(A) * alpha);
            a1 = -2 *     ((A-1) + (A+1)*cos(w0));
            a2 =          ((A+1) + (A-1)*cos(w0) - 2 * sqrt(A) * alpha);
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end
        
        function  [coefs] = sedea_rbjM_hs(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
            b0 =      A *  ((A+1) - (A-1)*cos(w0) + 2 * sqrt(A) * alpha);
            b1 = -2 * A *  ((A-1) - (A+1)*cos(w0));
            b2 =      A *  ((A+1) - (A-1)*cos(w0) - 2 * sqrt(A) * alpha);
            a0 =           ((A+1) + (A-1)*cos(w0) + 2 * sqrt(A) * alpha);
            a1 =  2     *  ((A-1) + (A+1)*cos(w0));
            a2 =           ((A+1) + (A-1)*cos(w0) - 2 * sqrt(A) * alpha);
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end
    end
end