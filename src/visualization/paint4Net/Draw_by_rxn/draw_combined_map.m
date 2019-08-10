function [involvedMets,deadEnds,deadRxns]=draw_combined_map(model,rxns,drawMap,direction,initialMet,excludeMets,sols,save,closev)

  if nargin<7 %if the number of arguments < 7, fill flux vector with x
        for v=1:length(model.rxns)
            flux(v)='x';
        end
  end
  
  flux = sols(:,3); % use the glyceraldehyde case for layout
    
  if nargin<8 %if the number of arguments < 8, save=false
        save=false;
    end
    if nargin<9
        closev=false; %if the number of arguments < 9, closev=false
    end

if ~isempty(flux) %if flux vector is not empty
    Start_time2=clock;
    Start_time=clock;
         
    %print start time for task
    Start_time=strcat(num2str(Start_time(1)),'_',num2str(Start_time(2)),'_',num2str(Start_time(3)),'_',num2str(Start_time(4)),'_',num2str(Start_time(5)),'_',num2str(Start_time(6)))
   
    if nargin<3 %if the number of arguments < 3, drawMap=false
        drawMap=false;
    end  

    if nargin<4 %if the number of arguments < 4, direction='struc'
        direction='struc';
    end

    if nargin<6 %if the number of arguments < 6, excludeMets{1}=''
        excludeMets{1}='';
    end

    %declares variables
    not_present_rxns{1}='0';
    rxns_index=1;

    for i=1:length(rxns) %cycle through the list of reactions
        
        exist=false;
        
        for j=1:length(model.rxns) %cycle through the reactions in the model
            
            if strcmp(model.rxns{j},rxns{i}) %if particular reaction exists in the model
               j=length(model.rxns); 
               exist=true;
            end

            if j==length(model.rxns) && exist~=true % if partucular reaction does not exist in the model, add to the list of not present reactions
               not_present_rxns{rxns_index}=rxns{i};
               rxns_index=rxns_index+1;
            end
            
        end
        
    end

    if not_present_rxns{1}=='0'; %if all reactions were present in the model 

        RxIDs=findRxnIDs(model,rxns); %find reaction IDs in the model
        
        %declare variables
        involvedMets='';
        mets_index=1;    

        for y=1:length(model.mets) %cycle through the metabolites in the model
            
            for x=1:length(RxIDs) %cycle through the reaction IDs
                
                switch direction %check the value of the variable direction
                    
                    case 'sub' %in cace of direction = 'sub'
                         
                        if model.S(y,RxIDs(x))<0 && flux(RxIDs(x))>1e-9 || model.S(y,RxIDs(x))>0 && flux(RxIDs(x))<-1e-9 %if the metabolite in the S matrix has negative coefficient, but the flux is positive or the opposite
                            
                            if ~strcmp(excludeMets{1},'') %if the list of excludable metabolites is not empty
                                
                                exclude=false;
                                
                                for v=1:length(excludeMets) %cycle through the metabolites in the list of excludable metabolites
                                    
                                    if strcmp(excludeMets{v},model.mets{y}) %if particular metabolite exists in the list of excludable metabolites
                                        exclude=true;
                                        break;
                                    end
                                    
                                end
                                
                                if exclude==false %if must not be excluded from visualization, add to the list of involved metabolites
                                    involvedMets{mets_index,1}=model.mets{y};
                                    mets_index=mets_index+1;    
                                end
                                
                            else %if the list of excludable metabolites is empty, add to the list of involved metabolites
                                
                                involvedMets{mets_index,1}=model.mets{y};
                                mets_index=mets_index+1;
                                
                            end
                            break;
                            
                        end
                        
                    case 'prod' %in cace of direction = 'prod'
                        
                        if model.S(y,RxIDs(x))>0 && flux(RxIDs(x))>1e-9 || model.S(y,RxIDs(x))<0 && flux(RxIDs(x))<-1e-9 %if the metabolite in the S matrix has positive coefficient and the flux is positive or the opposite
                            
                            if ~strcmp(excludeMets{1},'') %if the list of excludable metabolites is not empty
                                
                                exclude=false;
                                
                                for v=1:length(excludeMets) %cycle through the metabolites in the list of excludable metabolites
                                    
                                    if strcmp(excludeMets{v},model.mets{y}) %if particular metabolite exists in the list of excludable metabolites
                                        exclude=true;
                                        break;
                                    end
                                    
                                end
                                
                                if exclude==false %if must not be excluded from visualization, add to the list of involved metabolites
                                    involvedMets{mets_index,1}=model.mets{y};
                                    mets_index=mets_index+1;    
                                end
                                
                            else %if the list of excludable metabolites is empty, add to the list of involved metabolites
                                
                                involvedMets{mets_index,1}=model.mets{y};
                                mets_index=mets_index+1;
                                
                            end                
                            break;
                            
                        end
                        
                    case 'struc' %in cace of direction = 'struc'
                        
                        if model.S(y,RxIDs(x))<0 || model.S(y,RxIDs(x))>0 %if the metabolite in the S matrix has negative or possitive coefficient  
                            
                            if ~strcmp(excludeMets{1},'') %if the list of excludable metabolites is not empty
                                
                                exclude=false;
                                
                                for v=1:length(excludeMets) %cycle through the metabolites in the list of excludable metabolites
                                    
                                    if strcmp(excludeMets{v},model.mets{y}) %if particular metabolite exists in the list of excludable metabolites
                                        exclude=true;
                                        break;
                                    end
                                    
                                end
                                
                                if exclude==false %if must not be excluded from visualization, add to the list of involved metabolites
                                    involvedMets{mets_index,1}=model.mets{y};
                                    mets_index=mets_index+1;    
                                end
                                
                            else %if the list of excludable metabolites is empty, add to the list of involved metabolites
                                
                                involvedMets{mets_index,1}=model.mets{y};
                                mets_index=mets_index+1;
                                
                            end                
                            break;
                            
                        end
                    case 'both' %in cace of direction = 'both'
                        %if (the metabolite in the S matrix has negative coefficient, but the flux is positive or the opposite) or (the metabolite in the S matrix has positive coefficient and the flux is positive or the opposite)
                        if model.S(y,RxIDs(x))<0 && flux(RxIDs(x))>1e-9 || model.S(y,RxIDs(x))>0 && flux(RxIDs(x))<-1e-9 || model.S(y,RxIDs(x))>0 || model.S(y,RxIDs(x))>0 && flux(RxIDs(x))>1e-9 || model.S(y,RxIDs(x))<0 && flux(RxIDs(x))<-1e-9
                            
                            if ~strcmp(excludeMets{1},'') %if the list of excludable metabolites is not empty
                                
                                exclude=false;
                                
                                for v=1:length(excludeMets) %cycle through the metabolites in the list of excludable metabolites
                                    
                                    if strcmp(excludeMets{v},model.mets{y}) %if particular metabolite exists in the list of excludable metabolites
                                        exclude=true;
                                        break;
                                    end
                                    
                                end
                                
                                if exclude==false %if must not be excluded from visualization, add to the list of involved metabolites
                                    involvedMets{mets_index,1}=model.mets{y};
                                    mets_index=mets_index+1;    
                                end
                                
                            else %if the list of excludable metabolites is empty, add to the list of involved metabolites
                                
                                involvedMets{mets_index,1}=model.mets{y};
                                mets_index=mets_index+1;
                                
                            end                
                            break;
                            
                        end
                        
                end
                
            end
            
        end
       
        if nargin>5 && ~strcmp(initialMet{1},'') %if the number of arguments > 5 and the argument initialMet is not empty
            
            if ~ismember(initialMet{1},involvedMets) %if initialMet is not already in the list of involved metabolites, add initalMet to the list of involved metabolotes
                involvedMets{length(involvedMets)+1,1}=initialMet{1};
            end 
            
        end

        MetIDs=findMetIDs(model,involvedMets); %find metabolite IDs in the model
    
        for x=1:length(RxIDs) %cycle through the reaction IDs
            
            for y=1:length(MetIDs) %cycle through the metabolite IDs
                
                if model.S(MetIDs(y),RxIDs(x))==0 || model.ub(RxIDs(x))==0&&model.lb(RxIDs(x))==0%if the metabolite in the S matrix has coefficient=0 and botj bounds =0, add zeros into conectivity matrix 
                    matrix(y,x+length(MetIDs))=0;
                    matrix(x+length(MetIDs),y)=0;                    
                elseif model.S(MetIDs(y),RxIDs(x))<0 && flux(RxIDs(x))>1e-9 %if the metabolite in the S matrix has negative coefficient, but the flux is positive, add native direction entry into conectivity matrix 
                    matrix(y,x+length(MetIDs))=1;
                    matrix(x+length(MetIDs),y)=0;                    
                elseif model.S(MetIDs(y),RxIDs(x))<0 && flux(RxIDs(x))<-1e-9 %if the metabolite in the S matrix has negative coefficient and the flux is negative, add not-native direction entry into conectivity matrix
                    matrix(y,x+length(MetIDs))=0;
                    matrix(x+length(MetIDs),y)=1;                    
                elseif model.S(MetIDs(y),RxIDs(x))>0 && flux(RxIDs(x))>1e-9 %if the metabolite in the S matrix has positive coefficient and the flux is positive, add native direction entry into conectivity matrix
                    matrix(length(MetIDs)+x,y)=1;
                    matrix(y,length(MetIDs)+x)=0;                    
                elseif model.S(MetIDs(y),RxIDs(x))>0 && flux(RxIDs(x))<-1e-9 %if the metabolite in the S matrix has positive coefficient, but the flux is negative, add not-native direction entry into conectivity matrix
                    matrix(length(MetIDs)+x,y)=0;
                    matrix(y,length(MetIDs)+x)=1;                     
                elseif model.S(MetIDs(y),RxIDs(x))<0 && model.ub(RxIDs(x))>1e-9%if there is no flux and the metabolite in the S matrix has negative coefficient, and u.bound >0, add native direction entry into conectivity matrix
                    matrix(y,x+length(MetIDs))=1;
                    matrix(x+length(MetIDs),y)=0;
                elseif model.S(MetIDs(y),RxIDs(x))<0 && model.ub(RxIDs(x))==0&& model.lb(RxIDs(x))<-1e-9%if there is no flux and the metabolite in the S matrix has negative coefficient, and u.bound=0, but l.bound<0, add non-native direction entry into conectivity matrix
                    matrix(length(MetIDs)+x,y)=1;
                    matrix(y,length(MetIDs)+x)=0;
                elseif model.S(MetIDs(y),RxIDs(x))>0 && model.ub(RxIDs(x))>1e-9%if there is no flux and the metabolite in the S matrix has positive coefficient, and u.bound >0, add native direction entry into conectivity matrix
                    matrix(length(MetIDs)+x,y)=1;
                    matrix(y,length(MetIDs)+x)=0;  
                elseif model.S(MetIDs(y),RxIDs(x))>0 && model.ub(RxIDs(x))==0&& model.lb(RxIDs(x))<-1e-9%if there is no flux and the metabolite in the S matrix has positive coefficient, and u.bound=0, but l.bound<0, add non-native direction entry into conectivity matrix
                    matrix(y,x+length(MetIDs))=1;
                    matrix(x+length(MetIDs),y)=0;
                end
                
            end
            
        end

        %declare variables
        MetsNodesIDs=[1:length(MetIDs)];
        DeadMetsID(1)=-1;
        
        for v=1:length(MetsNodesIDs) %cycle through the metabolite node IDs
            
            if sum(matrix(v,:)) + sum(matrix(:,v))<2 %if the sum in the conectivity matrix < 2 (0 or 1 incoming/outgoing arrow), add to the list of dead end metabolites
                
                if DeadMetsID(1)==-1
                    DeadMetsID(1)=v;
                    deadEnds{1,1}=involvedMets{v,1};
                else
                    DeadMetsID(length(DeadMetsID)+1)=v;
                    deadEnds{length(deadEnds)+1,1}=involvedMets{v,1};
                end
                
            %if the sum in the conectivity matrix > 1 (at least 2 incoming and 0 outgoing arrows or at least 2 outgoing and 0 incoming arrows) it could be dead end metabolite  
            elseif sum(matrix(v,:)) + sum(matrix(:,v))>1 && sum(matrix(v,:))>0 && sum(matrix(:,v))==0 || sum(matrix(v,:)) + sum(matrix(:,v))>1 && sum(matrix(v,:))==0 && sum(matrix(:,v))>0
                
                DeadMet=true;            
                
                for w=length(MetsNodesIDs)+1:length(matrix) %cycle through the conectivity matrix starting from the first reaction index
                    
                    if matrix(v,w)==1 &&model.lb(RxIDs(w-length(MetsNodesIDs)))<-1e-9&&model.ub(RxIDs(w-length(MetsNodesIDs)))>1e-9 %if the coefficient in the conectivity matrix = 1 and this reaction in the model is reversible, then it is not a dead end metabolite                 
                        DeadMet=false;
                        break;
                    end
                    
                end

                if DeadMet==true
                    
                    for w=length(MetsNodesIDs)+1:length(matrix) %cycle through the conectivity matrix starting from the first reaction index
                        
                        if matrix(w,v)==1 && model.lb(RxIDs(w-length(MetsNodesIDs)))<-1e-9&&model.ub(RxIDs(w-length(MetsNodesIDs)))>1e-9  %if the coefficient in the conectivity matrix = 1 and this reaction in the model is reversible, then it is not a dead end metabolite                                
                            DeadMet=false;
                            break;
                        end
                        
                    end

                    if DeadMet==true %if still dead end, add to the list of dead end metabolites
                        
                        if DeadMetsID(1)==-1
                               DeadMetsID(1)=v;
                               deadEnds{1,1}=involvedMets{v,1};
                        else
                               DeadMetsID(length(DeadMetsID)+1)=v;
                               deadEnds{length(deadEnds)+1,1}=involvedMets{v,1};
                        end
                        
                    end
                    
                end  
                
            end
            
        end   

        if DeadMetsID(1)==-1 %if no dead ends
            deadEnds='No dead ends';
        end
        
        DeadRxnsID(1)=-1;%declare variable
        DeadRxnsID2(1)={'-1'};%declare variable
        
        %performe FVA
%         [minFlux,maxFlux] = fluxVariability(model,100,'max',rxns);
%         
%         for v=1:length(minFlux)
%             if minFlux(v)<1e-9&maxFlux(v)==minFlux(v)&minFlux(v)>-1e-9
%                  if strcmp(DeadRxnsID2(1),'-1')
%                      DeadRxnsID2(1)=rxns(v);                                          
%                  else
%                       DeadRxnsID2(length(DeadRxnsID2)+1)=rxns(v);                                              
%                 end
%             end
%         end     
%         
%         deadRxns=DeadRxnsID2';  %return list of dead reactions 

        if strcmp(drawMap,'true') %if drawMap=true
            
% CREATE NODES STARTING WITH METABOLITES            

            ids=involvedMets;

            if nargin>5 && ~strcmp(initialMet{1},'') %if the number of arguments > 5 and the argument initialMet is not empty
                
                for v=1:length(ids) %cycle through the node ids
                    
                    if strcmp(initialMet,ids{v})
                        initialMetID=v; %get the initial met ID
                    end
                    
                end
                
            end

            
% ADD REACTION NODES            
            for v=1:length(RxIDs) %cycle through the reaction ids
                %ids{length(ids)+1}=num2str(flux(RxIDs(v))); %node id=reaction abbr+flux rate           
                ids{length(ids)+1}=strcat(model.rxns{RxIDs(v)},' (',num2str(flux(RxIDs(v))),')'); %node id=reaction abbr+flux rate           
                RxnsNodesIDs(v)=length(ids); %get all reaction nodes IDs
            end

            ids=ids';      

            map = biograph(matrix,ids,'ShowTextInNodes','Label','LayoutType', 'hierarchical',...
                'LayoutScale',0.1,'Scale',1,'EdgeFontSize',12,'NodeAutoSize','off') %create a biograph object EdgeType/ 'LayoutType', 'evuilibrium'

% SET GENERAL METABOLITE REACTION NODE PROPERTIES            
            
            for m=1:length(involvedMets) %cycle through the involved mets
                
                for n=1:length(model.mets) %cycle through the mets in the model
                    
                    if strcmp(involvedMets{m},model.mets{n}) %if particular metabolite has been found in the model
                        
                        try %try to set label=metabolite name
                            set(map.nodes(m), 'Label', model.metNames{n},'Shape','box','Size',[150,15],'Color',[1 1 0.9],'LineColor',[0.1 0.1 0.1])                           
                        catch info
                        end
                        
                        try %try to set description=metabolite charged formula                                                   
                            set(map.nodes(m), 'Description', strcat('Charged formula: ',model.metFormulas{n}))
                        catch info
                        end
                        
                    end
                    
                end
                
            end

            
% SET GENERAL REACTION NODE PROPERTIES

            for m=1:length(RxIDs) %cycle through the reaction IDs
                
                for n=1:length(model.rxns) %cycle through the reactions in the model
                    
                    if strcmp(model.rxns{RxIDs(m)},model.rxns{n}) %if particular reaction has been found in the model
                        
                        try %try to set label=reaction name
                            if abs(sols(RxIDs(m),1)) < 0.001
                                flux1 = '0';
                            else
                                flux1 = num2str(abs(sols(RxIDs(m),1)),3);
                            end
                            if abs(sols(RxIDs(m),2)) < 0.001
                                flux2 = '0';
                            else
                                flux2 = num2str(abs(sols(RxIDs(m),2)),3);
                            end
                            if abs(sols(RxIDs(m),3)) < 0.001
                                flux3 = '0';
                            else
                                flux3 = num2str(abs(sols(RxIDs(m),3)),3);
                            end
                            reactionlabel = strcat({'(1) '}, flux1,{', (2) '}, flux2,{', (3) '},flux3);
                            set(map.nodes(length(involvedMets)+m), 'Label', string(reactionlabel),'Color',[1,1,1],'LineColor',[1,1,1],'Size',[150,10])
                        catch info
                        end
                        
                        try %try to set description=reaction evuation 
                            strtmp=printRxnFormula(model,model.rxns{RxIDs(m)},false);                    
                            set(map.nodes(length(involvedMets)+m), 'Description', strcat('Reaction: ',strtmp{1},' lb=',num2str(model.lb(RxIDs(m))),' ub=',num2str(model.ub(RxIDs(m)))))
                        catch info
                        end
                        
                    end
                    
                end
                
            end

% SET EDGE PROPERTIES
            
            max_flux=max(abs(flux(RxIDs))); %get max absolute flux value from the input vector

            for p=length(MetsNodesIDs)+1:length(ids) %cycle through the nodes IDs starting from the first reaction node
                    edges1 = getedgesbynodeid(map,ids(p),[]); %get all source edges connected to particular reaction node 
                    edges2 = getedgesbynodeid(map,[],ids(p)); %get all sink edges connected to particular reaction node  
                    tot_edges=[edges1;edges2]; %concatenate edges                    
                    set(tot_edges,'LineColor',[0 0 0]); %set green color to found edges
%                    set(tot_edges,'LineWidth',abs(6*((flux(RxIDs(p-length(MetsNodesIDs))))/max_flux)+1)); %set thichness of edge proportional to max flux 
                    set(tot_edges,'LineWidth',2);          
%                 if flux(RxIDs(p-length(MetsNodesIDs)))>1e-9 && flux(RxIDs(p-length(MetsNodesIDs)))~='x' %if flux is positive and not evual to x                    
%                     edges1 = getedgesbynodeid(map,ids(p),[]); %get all source edges connected to particular reaction node 
%                     edges2 = getedgesbynodeid(map,[],ids(p)); %get all sink edges connected to particular reaction node  
%                     tot_edges=[edges1;edges2]; %concatenate edges                    
%                     set(tot_edges,'LineColor',[0 0 0]); %set green color to found edges
%                     set(tot_edges,'LineWidth',abs(6*((flux(RxIDs(p-length(MetsNodesIDs))))/max_flux)+1)); %set thichness of edge proportional to max flux 
%                 elseif flux(RxIDs(p-length(MetsNodesIDs)))<-1e-9 && flux(RxIDs(p-length(MetsNodesIDs)))~='x' %if flux is negative and not evual to x 
%                     edges1 = getedgesbynodeid(map,ids(p),[]); %get all source edges connected to particular reaction node
%                     edges2 = getedgesbynodeid(map,[],ids(p)); %get all sink edges connected to particular reaction node                 
%                     tot_edges=[edges1;edges2]; %concatenate edges 
%                     set(tot_edges,'LineColor',[0 0 0]); %set green color to found edges
%                     set(tot_edges,'LineWidth',abs(6*((flux(RxIDs(p-length(MetsNodesIDs))))/max_flux)-1)); %set thichness of edge proportional to max flux 
%                 end
            end  
            
               h=view(map); %generate layout            
%               set(h.nodes(RxnsNodesIDs),'Color',[1,1,1],'LineColor',[1,1,1],'Size',[0.1,0.1]); %set color to reaction nodes
%               set(h.nodes(MetsNodesIDs),'Shape','ellipse','Size',[1000,10]); %set shape of ellipse to metabolite nodes
               dolayout(map)
                           
                  
             for v=1:length(DeadRxnsID2)
                 for z=length(involvedMets)+1:length(ids)
                    tmpR=[];
                    for q=1:length(ids{z})
                        if ~strcmp(ids{z}(q),' ')
                            tmpR=[tmpR ids{z}(q)];
                        else
                            break
                        end
                    end
                      if strcmp(DeadRxnsID2{v},tmpR)
                             if DeadRxnsID(1)==-1
                                 DeadRxnsID(1)=z;     
                                 break
                             else
                                  DeadRxnsID(length(DeadRxnsID)+1)=z;     
                                  break
                             end
                      end
                end
             end   

            if DeadMetsID(1)~=-1 %if dead end mets exist
                set(h.nodes(DeadMetsID),'Color',[1 1 1]); %set red color to dead end metabolite nodes
            end           
            
            ExRxnsID=-1;%declare variable
            
            for q=length(MetsNodesIDs)+1:size(matrix,2) %cycle through the nodes IDs starting from the first reaction node
                
                if sum(matrix(q,:)) + sum(matrix(:,q))<2 %if sum of edges to one reactio node is less then 2, then it could be exchange reaction
                    
                    if ExRxnsID(1)==-1
                        ExRxnsID(1)=q;
                    else
                        ExRxnsID(length(ExRxnsID)+1)=q;
                    end
                    
                end
                
            end
            
            if ExRxnsID(1)~=-1 %if exchange reactions exist
                set(h.nodes(ExRxnsID),'Color',[173/255, 1, 47/255]); %set red color to exchange reaction nodes                
            end
   
            if nargin>5 && ~strcmp(initialMet{1},'') %if the number of arguments > 5 and the argument initialMet is not empty
                set(h.nodes(initialMetID),'Color',[0 1 0]); %set green color to initial met node          
            end            
        
%             if DeadRxnsID(1)~=-1 %if dead reactions exist
%                 set(h.nodes(DeadRxnsID),'Color',[1 0 0]); %set red color to dead reaction nodes
%             end
            
            if ExRxnsID(1)~=-1 %if exchange reactions exist                
                set(h.nodes(ExRxnsID),'Shape','circle'); %set shape of ellipse to metabolite nodes
            end
            
            if strcmp(save,'true') %if save=true
                f = get(h.hgAxes, 'Parent');
                clock_tmp=clock;        
                print(f, '-djpeg', strcat(num2str(clock_tmp(1)),'_',num2str(clock_tmp(2)),'_',num2str(clock_tmp(3)),'_',num2str(clock_tmp(4)),'_',num2str(clock_tmp(5)),'_',num2str(clock_tmp(6)),'.jpg'));
        
            end
            
            if strcmp(closev,'true') %if closev=true
                child_handles = allchild(0);
                names = get(child_handles,'Name');
                k = find(strncmp('Biograph Viewer', names, 15));
                close(child_handles(k))
            end
            
        end   

    else %if in the input list was at least one unexisting reaction in the model
        
        disp('The following reactions are not present in the model:')
        not_present=not_present_rxns'
        involvedMets='No involved mets';
        deadEnds='No dead ends';
        
    end
    
    %print end and total time for task
    End_time=clock;
         Total_time=End_time-Start_time2;
         End_time=strcat(num2str(End_time(1)),'_',num2str(End_time(2)),'_',num2str(End_time(3)),'_',num2str(End_time(4)),'_',num2str(End_time(5)),'_',num2str(End_time(6)))
        Total_time=strcat(num2str(Total_time(1)),'_',num2str(Total_time(2)),'_',num2str(Total_time(3)),'_',num2str(Total_time(4)),'_',num2str(Total_time(5)),'_',num2str(Total_time(6)))
else %if flux vector is empty
    
        disp('The flux vector is empty')
        involvedMets='No involved mets';
        deadEnds='No dead ends';
        End_time=clock;
         Total_time=End_time-Start_time2;
         End_time=strcat(num2str(End_time(1)),'_',num2str(End_time(2)),'_',num2str(End_time(3)),'_',num2str(End_time(4)),'_',num2str(End_time(5)),'_',num2str(End_time(6)))
        Total_time=strcat(num2str(Total_time(1)),'_',num2str(Total_time(2)),'_',num2str(Total_time(3)),'_',num2str(Total_time(4)),'_',num2str(Total_time(5)),'_',num2str(Total_time(6)))
        
end
end