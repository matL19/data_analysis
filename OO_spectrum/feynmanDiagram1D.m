classdef feynmanDiagram1D
    
    properties
        fun; %for function handle
        R; %calculated response
    end
    
    methods
        function obj = feynmanDiagram(fun,~)
            if nargin>0
                obj.fun = fun;
            end
        end
        
        function obj = calcResponseTime(obj,T1,~,~)
            obj.R = obj.fun(T1);
        end
        
        function obj = timeToFreq(obj,n_zp)
            % do the FFT and flips of spectrum and return the real
            % component (used for global fitting)
            
            obj.R = sgrsfft(obj.R(1,:),n_zp);
            obj.R = fftshift(real(obj.R));

        end
        
        function obj = timeToFreqComplex(obj,n_zp)
            % do the FFT and flips of spectrum and return the complex
            % spectrum (not used for global fitting)
            obj.R = sgrsfft(obj.R(1,:),n_zp);
            obj.R = fftshift(obj.R);
        end
        
    end
    
end
