function [n,lambda] = Estimate_Normal_CPS_MicroBRDF_General_Trigo(H,L,I,maxItr,n0) 
    if nargin < 4
        maxItr = 100;
    end
    
    if nargin == 5
        % debug mode
        [x0(2),x0(1),x0(3)] =cart2sph(n0(1),n0(2),n0(3));
        options = optimset('display','off','MaxIter',100); 
        lb = [-2*pi,-2*pi,0.0001];
        ub = [2*pi,2*pi,1]; 
        [xest,res] = lsqnonlin(@(x)objectivefunction(x,H,L,I),x0,lb,ub,options);
        n = [cos(xest(1))*cos(xest(2)); cos(xest(1))*sin(xest(2)); sin(xest(1))];
        if sum(find(L*n>0)) < sum(find(-L*n>0)) 
            n = -n;
        end
        lambda = xest(3);        
        return;
    end  
    
    nLambert = L\I;
    nLambert = nLambert/norm(nLambert); 
    
    %Note the definition of cart2sph, sph2cart in matlab 
    [x0DiffusiveTrigo(2),x0DiffusiveTrigo(1),x0DiffusiveTrigo(3)] = cart2sph(nLambert(1),nLambert(2),nLambert(3));
    
    [nSpecular,~] = Estimate_Normal_CPS_MicroBRDF_Specular(H,L,I); 
    nSpecular = nSpecular/norm(nSpecular);
    [x0SpecularTrigo(2),x0SpecularTrigo(1),x0SpecularTrigo(3)] =cart2sph(nSpecular(1),nSpecular(2),nSpecular(3));
    x0SpecularTrigo(3) = 0; %lambda
    
    if isnan(x0SpecularTrigo) > 0
       xoneshot = x0DiffusiveTrigo; 
    elseif norm(objectivefunction(x0SpecularTrigo,H,L,I)) < norm(objectivefunction(x0DiffusiveTrigo,H,L,I))
       xoneshot = x0SpecularTrigo; 
    else 
       xoneshot = x0DiffusiveTrigo;
    end
        
    options = optimset('display','off','MaxIter',maxItr);    
    lb = [-2*pi,-2*pi,0.0001];
    ub = [2*pi,2*pi,1]; 

    [xdifussive,resnormDiffusive] = lsqnonlin(@(x)objectivefunction(x,H,L,I),x0DiffusiveTrigo,lb,ub,options);
    xest = real(xdifussive); 
       
    if nnz(isnan(x0SpecularTrigo)) == 0
        [xspecular,resnormSpecular] = lsqnonlin(@(x)objectivefunction(x,H,L,I),x0SpecularTrigo,lb,ub,options);
        if resnormSpecular < resnormDiffusive
           xest = real(xspecular);
        end
    end
    
    if norm(objectivefunction(xoneshot,H,L,I)) < norm(objectivefunction(xest,H,L,I))
        xest = xoneshot;
    end
    
    n = [cos(xest(1))*cos(xest(2)); cos(xest(1))*sin(xest(2)); sin(xest(1))];
%     xoneshot = [cos(xoneshot(1))*cos(xoneshot(2)); cos(xoneshot(1))*sin(xoneshot(2));sin(xoneshot(1))];
    
    %double check the sign
    if sum(find(L*n>0)) < sum(find(-L*n>0)) 
        n = -n;
    end

    lambda = xest(3);

%     normal(4:6) = xoneshot(1:3)/norm(xoneshot(1:3));
end

function diff = objectivefunction(x,H,L,I)

    n = [cos(x(1))*cos(x(2));cos(x(1))*sin(x(2));sin(x(1))];
    
    
    lambda = x(3); 
    nl = L*n;
    
%     selected = nl < 0;
    
    nl2 = nl.^2; 
    nh2 = (H*n).^2; 

    vec = nl./(1-(1-lambda)*nh2).^2./sqrt(lambda+(1-lambda)*nl2);
    diff = 100000*(I - vec*(vec.'*I/(vec.'*vec)));
    
%     diff(selected) = 0;
    
end