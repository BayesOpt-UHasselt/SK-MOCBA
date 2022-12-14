function varargout = DTLZ7(Operation,Global,input)
    switch Operation
        case 'init'
            Global.M        = 3;
            Global.D        = Global.M + 19;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            M      = Global.M;
            
            PopObj          = zeros(size(PopDec,1),M);
            g               = 1+9*mean(PopDec(:,M:end),2);
            PopObj(:,1:M-1) = PopDec(:,1:M-1);
            PopObj(:,M)     = (1+g).*(M-sum(PopObj(:,1:M-1)./(1+repmat(g,1,M-1)).*(1+sin(3*pi.*PopObj(:,1:M-1))),2));
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            interval     = [0,0.251412,0.631627,0.859401];
            median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
            X            = ReplicatePoint(input,Global.M-1);
            X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
            X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
            f            = [X,2*(Global.M-sum(X/2.*(1+sin(3*pi.*X)),2))];
            varargout    = {f};
    end
end

function W = ReplicatePoint(SampleNum,M)
    if M > 1
        SampleNum = (ceil(SampleNum^(1/M)))^M;
        %Gap       = 0:1/(SampleNum^(1/M)-1):1;
        eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
        eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
    else
        W = (0:1/(SampleNum-1):1)';
    end
end
