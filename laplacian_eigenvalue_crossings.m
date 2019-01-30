function [E, LpLn, crossings, F] = laplacian_eigenvalue_crossings(A, figs, epsilon, threshold)
    % for a given signed graph with adjacency matrix A, this finds all
    % possible inertias with maximum nullity two.
    % inputs:
    %   A: the signed adjacency matrix
    %   figs: boolean, whether to output figure handles, default: 1
    %   epsilon: the step size for t, defauly: 0.01
    %   threshold: to distinguish small numbers from zero, default: 10^(-10)
    
    % outputs:
    %   E: list of eigenvalues, each column is the n eigenvalues for one t
    %   LpLn: the positive and the negative matrix used. first n columns
    %       are positive matrix and last n columns are LN
    %       for example, if you want to print all spectra with double zeros:
    %
    %       LP = LpLn(:,1:size(A,1));
    %       LN = LpLn(:,size(A,1)+1:end);
    %       for t = crossings
    %           L = LP - t * LN;
    %           eig(L)
    %       end
    %
    %   crossings: t's for which the crossings happen 
    %   F: figure handles
    
    if nargin < 2
        figs = 1; % whether or not output figures
        if nargin < 3
            epsilon = .01; % I should be able to figure this out by using the min gap result
            if nargin < 4
                threshold = 10^(-10); % to avoid numerical errors around 0
            end
        end
    end
    
    n = size(A,1);
    
    AP = A > 0; % fix positive edges with weight 1
    AP = AP - diag(diag(AP)); % get rid of the possible diagonal entries
    
    R = rand(size(A)); % generate a random matrix for weights of negative edges
    R = 1 + (R + R') / 2; % symmetrize it and shift to stay away from zero
    AN = R .* (A < 0); % put the random weights on the negative edges
    AN = AN - diag(diag(AN)); % get rid of any possible diagonal entries
    
    GP = graph(AP); % construct the positive graph
    GN = graph(AN); % construct the negative graph
    
    CP = length(unique(conncomp(GP))); % number of connected components of the positive graph
    CN = length(unique(conncomp(GN))); % number of connected components of the negative graph
    tau = n + 1 - CP - CN; % the "felxibility" of graph, this is the number of corssings
    
    LP = laplacian(AP); % the Laplacian of the positive graph
    LN = laplacian(AN); % the Laplacian of the negative graph
    
    f1 = figure(); % draw graph
        G = graph(AP - AN);
        edge_colors = G.Edges.Weight;
        edge_colors(edge_colors < 0) = 2;
        cmp = distinguishable_colors(2);
        plot(G,'NodeColor','k','EdgeColor',cmp(edge_colors,:),'linewidth',2, 'layout','circle')
        title({['n = ' num2str(n) ', c_+ = ' num2str(CP) ', c_- = ' num2str(CN) ', \tau = ' num2str(tau)]})

    t = 0; % start with no negative edges
    ncrossed = 0; % number of crossed eigenvalues
    crossings = []; % t's for which crossings happen
    E = eig(LP); % list of eigenvalues, each t gets a column of E
    while ncrossed < tau
        t = t + epsilon;
        L = LP - t * LN;
        E(:,end+1) = eig(L);
        if sum(E(:,end-1) > threshold) > sum(E(:,end) > threshold) % if something crossed
            ncrossed = ncrossed + 1;
            crossings(end+1) = t - epsilon/2;
        end
    end
    buffer = 0;
    while buffer < 100 % keep going for another 100 steps
        t = t + epsilon;
        L = LP - t * LN;
        E(:,end+1) = eig(L);
        buffer = buffer + 1;
    end

    f2 = figure(); % draw eigenvalues
        plot(0:epsilon:t+epsilon,-E','linewidth',1.5)
        xlim([0,t])
        ylim([-max(E(:)),max(E(:))])
        set(gca,'XTick', crossings);
        set(gca,'YTick', [0]);
        xtickangle(90);
        title({['\tau = ' num2str(tau)],[num2str(CN-1) ' \leq n_- \leq ' num2str(n-CP)], [num2str(CP-1) ' \leq n_+ \leq ' num2str(n-CN)]})
        hold on
        vline(crossings, 'k')
        
    if figs
        F = [f1, f2];
    else
        F = [];
    end
    LpLn = [LP, LN];
    E = -E; % this is to match their convention
    
end