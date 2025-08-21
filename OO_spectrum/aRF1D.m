classdef aRFWAOBnd < aRFBnd
    
    
    properties
        %
        %  Most properties should be defined in the superclass aRF
        %
        
        % enumerating diagrams
        diagram;
        
        dyn = additionalDynamics; %structure/class) for additional dynamics
        dynOther = additionalDynamics; %for dynamics shared paramters
        
        %
        % Additional things that might be fit parameters not in the
        % superclass
        %
        
        % transition dipole moments
        mu01sq;
        
    end
    
    methods
        function obj = aRFWAOBnd(options)
            obj@aRFBnd(options);
        end
        
        function obj = calcSpectrum(obj,p)
            obj = obj.updateFreeFitParams(p);
            obj = obj.setupFreqAxes;
            obj = obj.setupResponseFunctions;

            obj = obj.calcResponseFunctions(obj.paramStruct);                
            obj = obj.calcPhaseShift(obj.paramStruct);
            obj = obj.calcAnhShift(obj.paramStruct);
            obj = obj.calcTDM(obj.paramStruct);
            obj = obj.calcAdditionalDynamics;
            obj=obj.calcDiagramsFreq(obj.n_zp);
            obj = obj.addDiagrams;
            obj=obj.resample(ii);

        end
        
        function obj = makeResponseFunctions(obj,~)
            
        end
        function obj = setupResponseFunctions(obj)
            
            obj.n_diagrams = 1;
            
            %have to initialize the first one without an index (don't know
            %why)
            obj.diagram = feynmanDiagram();
            
        end
        
        function obj = calcResponseFunctions(obj,p)
            
            %update the lineshape function with new parameters
            obj.damping = obj.damping.updateG(p);
            g = obj.damping.g; %shortcut
            
            R =  exp(-g(obj.T1));
            
            obj.diagram(1).R = R;
        end
        
        function obj = calcPhaseShift(obj,~)
            
            % do nothing
            
        end
        
        function obj = calcTDM(obj,p)
            mu_01_2 = p.mu01sq;
            
            obj.diagram(1).R = mu_01_2^2.*obj.diagram(1).R;
        end
        
        function obj = calcAnhShift(obj,~)
            
            % do nothing
        end
        
        function obj = calcAdditionalDynamics(obj)
            p = obj.paramStruct;
            
            for ii = 1:length(obj.dyn)
                for jj = 1:length(obj.dyn(ii).fun_array)
                    f = obj.dyn(ii).fun_array{jj};
                    ind = obj.dyn(ii).ind_array{jj};
                    for kk = 1:length(ind)
                        obj.diagram(ind(kk)).R = f(obj.T1,p).*obj.diagram(ind(kk)).R;
                    end
                end
            end
            
            for ii = 1:length(obj.dynOther)
                for jj = 1:length(obj.dynOther(ii).fun_array)
                    f = obj.dynOther(ii).fun_array{jj};
                    ind = obj.dynOther(ii).ind_array{jj};
                    for kk = 1:length(ind)
                        obj.diagram(ind(kk)).R = f(obj.T1,p).*obj.diagram(ind(kk)).R;
                    end
                end
            end
        end
    end
end

