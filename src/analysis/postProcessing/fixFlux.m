function fluxes = fixFlux(model, solution, substrateRxn)

% Fix model solution by converting very high fluxes
% Substrate flux must be less than 100

    substrFlux = solution.x(strcmp(model.rxns, substrateRxn));
    highThreshold = abs(substrFlux)*10;
    
    tempFlux = solution.x;
    
    for i=1:length(tempFlux)
        if abs(tempFlux(i)) > 1000-highThreshold && abs(tempFlux(i)) ~= 1000
            if tempFlux(i) > 0
                tempFlux(i) = tempFlux(i) - 1000;
            elseif tempFlux(i) < 0
                tempFlux(i) = 1000 + tempFlux(i);
            end
        end
    end

    fluxes = tempFlux;
end