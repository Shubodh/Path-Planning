function fn_DstarLite( start,goal,grid,H,NC,BIG )
%Main D star Lite function
% BIG=100; %% an impportant constant, working at the place of infinity or highest number possible
nn=length(grid); %% nn = number of nodes
% the cost grid act as an imaginary grid that can be doctored according to use
c=grid;
pt = 0; % coordinates of mouse click
SIZEob=1;obs=zeros(SIZEob,nn); obhndl=zeros(1,nn);
slow=0.1;% should be less than 0.01
%% 
last=start;
[k,rhs,g,U]=fn_Initialize(nn,goal,BIG,H(start,goal));
% U(:,s)=CalKey(g(goal),rhs(s),k,heuristic(start,s));
ax=axes('XLim',[min(NC(1,1)-1,NC(nn,1)-1) max(NC(1,1)+1,NC(nn,1)+1)],'YLim',[min(NC(1,2)-1,NC(nn,2)-1) max(NC(1,2)+1,NC(nn,2)+1)]);
hold on;
axis square;
[U,g,rhs] = fn_CompShortestPath( U,g,rhs,k );
fn_SimulateGraph;
prev=start;
while start~=goal
    tic;
    if g(start)>=BIG        %% then there is no known path thus exit the code
        disp('path not found');
        break               
    end
%     prev=start;
    %%% start = arg min (s' = Succ(start)(c(start,s')+g(s')) )
    ss=find(c(start,:)<BIG); % make local list of all neighbours of current start node
    ss=setdiff(ss,prev); % remove the predecessor
    [~,i]=min(g(ss)+c(start,ss)); %% find the neighbour with min c+g in local list 'ss'
    start=ss(i);            %% convert from local list to global list
    path=start;
%     path=FindPath(path);  % find the path from current start to goal
%     pth=plot(NC(path,1),NC(path,2),'g'); % plot the path in green
    plot(NC([prev start],1),NC([prev start],2),'r'); % plot the path traveresed in red
%faltu%  start = find(Grid(start,:)==min(Grid(start,Grid(start,:)>0)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     x=NC([start goal],1); y=NC([start goal],2);
%     x=NC(start,1); y=NC(start,2);
    x=linspace(NC(prev,1),NC(start,1),slow*100);y=linspace(NC(prev,2),NC(start,2),slow*100);
    prev=start;
    %     UpdateGraph;
    for i=1:length(x)
    X=bsxfun(@plus,x(i),rc); Y=bsxfun(@plus,y(i),rs);
    set(h(1),'XData',X(1,:),'YData',Y(1,:));
    pause(0.0005);
%     pause(slow);
    end
    %%% Scan graph for changed edge costs
%     pobs=obs;
% keyboard;
   [ change,obs ] = fn_obstacle( NC,obs);
   change=unique(change);
    %%%%%%%%%%%%%%%%
    %%% If any edge cost changed
    if length(change)>1 % i.e. node shud exist, shudnt be at start or goal
        k=k+H(last,start);
        last=start;
        c=grid;
        ob=find(obs);
        c(ob,:)=BIG; c(:,ob)=BIG;        %% Update the edge costs c(u,v)
        for i=2:length(change)
        ss=find(grid(change(i),:)<BIG);    %% find all the neighbours of changed node
%             uu=change(i);
            uu=[change(i),ss,change(i)]; % ,find(grid(change(i),:)<BIG)
            for j1=1:length(uu)
                [rhs,U] = fn_UpdateVertex(U,g,rhs,uu(j1));
            end
        end
        [U,g,rhs] = fn_CompShortestPath( U,g,rhs,k );
    end
%     delete(pth);
    toc;
    
    scatter(NC(:,1),NC(:,2),10,'b','fill'); % display the Network
    incosis=find(g~=rhs);
    scatter(NC(incosis,1),NC(incosis,2),30,'r','fill'); 
    consis=setdiff(find(g==rhs),find(rhs>=100));
    scatter(NC(consis,1),NC(consis,2),30,'g','fill');
    
%     keyboard;
end
%%%%%%%%
function [U,g,rhs] = fn_CompShortestPath( U,g,rhs,k )
while (fn_LexComp(fn_BestKey(U),fn_CalKey(g(start),rhs(start),k,0)) || rhs(start)~=g(start))
    u=fn_Pop(U); % find the node with best key
    kold=U(:,u); % find the key of that node
    U(:,u)=BIG;   % i.e. delete the key of node u from the U
    key=fn_CalKey(g(u),rhs(u),H(u,start),k);
    if fn_LexComp(kold,key)
        U(:,u)=key; % i.e. insert the new key in 'u' node
    elseif g(u)>rhs(u)
           g(u)=rhs(u);        
           nhbr=find(grid(u,:)<BIG);  % ss vector will include all neighbours nodes of vertex u
%           ss=find(c(u,:)<BIG);
           for ii=1:length(nhbr)     % update all neighbpurs of node 'u'
              [rhs,U]= fn_UpdateVertex(U,g,rhs,nhbr(ii));
           end
     else
           g(u)=BIG;
           nhbr=[find(grid(u,:)<BIG),u];
           for ii=1:length(nhbr)
              [rhs,U]= fn_UpdateVertex(U,g,rhs,nhbr(ii));
           end
     end
end
end

function [rhs,U]= fn_UpdateVertex(U,g,rhs,uu)
if uu~=goal
    sss=find(grid(uu,:)<BIG);               % ss = all neighbours of uu
    rhs(uu)=min(g(sss)+c(sss,uu)');     % here g(sss) is row vector and gd(ss,uu) is column vector thus transpose of one of them is taken for their addition
    if rhs(uu)>BIG
        rhs(uu)=BIG; end
end
%%%%%%%%%%%%%%%%%%
%%% if u is from U den remove it
if any(U(:,uu))<BIG
    U(:,uu)=BIG;
end
%%%%%%%%%%%%%%%%%%
%%% if g(u)~=rhs(u) den update key (u) in U
key=fn_CalKey(g(uu),rhs(uu),k,H(start,uu));
if g(uu)~=rhs(uu) && fn_LexComp(key,U(:,uu))
    U(:,uu)=key;
end
end

    function fn_SimulateGraph()
        ang=0:0.01:2*pi;
        r=[0.5;0.25];    % radius of source and goal
        x=NC([start goal],1); y=NC([start goal],2);
        rc=r*cos(ang); rs=r*sin(ang);
        X=bsxfun(@plus,x,rc); Y=bsxfun(@plus,y,rs);
%         f=figure;
        
        scatter(NC(:,1),NC(:,2),10,'b','fill'); % display the Network
        h(1)=fill(X(1,:),Y(1,:),'r');   % red for robot
        fill(X(2,:),Y(2,:),'k');   % black for goal
    end
function [ change,obs ] = fn_obstacle( NC,obs )
%UNTITLED Summary of this function goes here
% Detailed explanation goes here
    change=0;
    pt=0;
    while 1
%         set(ax, 'ButtonDownFcn',@FindClick );
        pt=ginput(1);
        pt=round(pt);
%         display(pt);
        if pt(1)<NC(1,1) || pt(1)>NC(end,1) || pt(2)>NC(1,2) || pt(1)<NC(end,2)
           return 
        end
        node=0;
        if pt~=0
            [~,node]=ismember(pt,NC,'rows');   
            pt=0;
            if ~any([0 start goal]==node) % i.e. node shud exist, shudnt be at start or goal
                change=[change node];
                if ismember(node,obs) % if obstacle already exist den remove it
                    [~,j]=find(obs==node,1);% finds only one obs (1st one)
                    delete(obhndl(j));
                    obs(:,j)=0;    % remove from obs list and add to cost grid
                else    % else add the obs
                    obs(node)=node;
                    obhndl(node)=rectangle('Position',[NC(node,1)-SIZEob+0.5,NC(node,2)-SIZEob+0.5,SIZEob,SIZEob],'faceColor','k');
                end
            end
        end
    end
end
hold off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [k,rhs,g,U] = fn_Initialize( nn,goal,BIG,H )
rhs(1:nn)=BIG;
U=ones(2,nn)*BIG;
g=rhs;
k=0;
rhs(goal)=0;
s=goal;
U(:,s)=fn_CalKey(g(s),rhs(s),k,H);  % insert the goal key in U
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [key]=fn_CalKey(g,rhs,k,h)
key=[min(g,rhs)+h+k; min(g,rhs)]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [key]=fn_BestKey(U)
%% find the smallest key in 'U', if it is empty then return [infi;infi]
u=find(U(1,:)==min(U(1,:))); %% list of indices of smallest elements in 1st row
[~,i]=min(U(2,u)); %% find the index of 1st sallest element in local iist u
key=U(:,u(i));      %% here u(i) is the node with smallest key

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [c]=fn_LexComp(k1,k2)
% check if k1 is lexicographically smaller (i.e. better) than k2 or not
% if yes then it returns 1 else 0;- (1 means 'Yes' and 0 means 'No')
if k1(1)<k2(1) || (k1(1)==k2(1) && k1(2)<k2(2))
    c=1;        % yes k1<k2 (lexicographically)
    return
end
c=0;            % k1>k2 thus false
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [i]=fn_Pop(U)
%find the node with smallest key in U
u = find( U(1,:)==min(U(1,:)) ); %% list of indices of smallest elements in 1st row
[~,i]=min(U(2,u)); %% find the index of 1st smallest element in local iist u
i=u(i);
end