function fluxesMajor = fluxFollow(model, solution, substrateRxn, lowThreshold, highThreshold, excludeRxns)

% Determine significant fluxes and draw map

substrate = findMetsFromRxns(model,substrateRxn);

cont = true;
init = substrate;
mets = substrate;
rxns = substrateRxn;

while cont
    cont = false;
    tempMets = {};
    % Find the reactions from the metabolites
    [list, formulas] = findRxnsFromMets(model, init);
    IDs = findRxnIDs(model,list);
    % Get the reaction fluxes
    found = solution.x(IDs);
    % For each reaction
    for i=1:length(found)
        % If the flux is greater than the threshold and not maximum
        if abs(found(i))>lowThreshold && abs(found(i))<highThreshold
            % Also if it is not already in the rxn list (i.e. already
            % tested) and not among the exlucded rxns
            if all(~ismember(rxns,list(i)))&&all(~ismember(excludeRxns,list(i)))
                % Add the reactions to the reaction list
                rxns(end+1) = list(i);
                % Get the metabolites and create a new list of metabolites
                tempMets = [tempMets;findMetsFromRxns(model,list(i))];
                cont = true;
            end
        end
    end
    init = tempMets;
end


ID2 = findRxnIDs(model,rxns');
fluxes = solution.x(ID2);
formulas = printRxnFormula(model,'rxnAbbrList',rxns');
fluxesMajor = [rxns',formulas, num2cell(fluxes)];

end