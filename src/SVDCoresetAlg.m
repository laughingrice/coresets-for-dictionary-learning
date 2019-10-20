classdef SVDCoresetAlg < handle
    % SVDCORESET Construct the coreset for computing the SVD
    %
    % Steps:
    % 1. Randomly select as set B of beta points to use as the vectors
    %    spanning an approximate to the best low dimensional subspace S
    %
    % 2. Extract the set Q of n/2 closest points to S
    %    The subspace S approximates well it's closest n/2 points Q
    %    - compute the reduced QR decomposition on B to get an
    %      orthonormal basis to this subspace S
    %      [QS, ~] = qr(B, 0);
    %    - sqdists = sum((P - P*QS*QS').^2, 2)
    %
    % 3. Select from Q - Q*QS*QS' (the part of Q that is orthogonal to S)
    %    a set of k / epsilon^2 elements with probability
    %    proportional to squared distance from the subspace S
    %    - pr \propto sum((Q - Q*QS*QS').^2, 2)
    %
    % 4. Project Q onto the subspace S (Q*QS) and select another
    %    k / epsilon^2 points from Q with a probability proportional to either:
    %    a) The squared norm of the rows of U from the SVD decomposition
    %       [U, ~, ~] = svd(Q*QS, 0);
    %       Note that this can be computed faster, albeit with a somewhat lower
    %       stability (but good enough for our purposes) based on the PCA
    %       [V, D] = eig((Q*QS)'*(Q*QS));
    %       U = Q*QS*V*diag(1./sqrt(diag(D)))
    %    b) The squared norm of the rows of Q from the QR decomposition
    %       [Q, ~] = qr(Q*QS, 0);
    %
    %    - pr \propto sum(U.^2, 2) or pr \propto sum(Q.^2, 2)
    %
    %    - Note: This select based on normalized distance from center
    %      (circularize the distribution ellipse)
    %
    % 5. Throw out the Q points from P and continue with the rest of the points
    %    until the stopping criterion is reached

    properties
        beta = 5; % Number of samples for the bi-criteria stage. Should be O(k), where k is the dimension of the requested subspace
        minFunctionSize = 10; % Stopping term for the bi-criteria iteration
        bicriteriaRatio = 0.5; % How many elements to remove at each iteration of the bicriteria
        sampleSize = 50; % The number of elements to take into the coreset at each iteration. Should be equal to k / epsilon^2
    end
    
    methods
        % Constructor
        function obj = SVDCoresetAlg()
        end

        % Compute coreset for elements in matrix P
        % Input is a Matlab matrix with elements in rows
        % P is n x d, i.e n points in R^d
        function coreset = computeCoreset(obj, P)
            n = size(P, 1);
            coreset = [];

            while n >= obj.minFunctionSize % Loop terminating condition
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 1. Randomly sample beta elements
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                ind = randsample(n, obj.beta, true);
                ind = sort(ind);

                % accumulate duplicates
                % just removing duplicates is actually enough at this stage, or
                % we can leave this to the QR stage later on
                h = histc(ind, ind);
                ind = ind(h > 0);

                B = P(ind, :);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 2. Select the set Q of the obj.sampleSize closest points
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Project the set P onto the subspace spanned by B
                
                % a) Compute an orthogonal basis for the subspace
                [QS, ~] = qr(B', 0);
                % b) Project the set of all points onto the subspace
                PQS = P*QS;
                % c) Take the part orthogonal to the the subspace, this is
                %    needed both for distance calculation and for the points to select from
                OPQS = P - PQS*QS';
                % d) Compute the squared distances from the subspace 
                sqdists = sum(OPQS .* OPQS, 2);
                % e) Use obj.sampleSize closest and then remove them from the
                %    set for the next iteration
                [sqdists, mIndexes] = sort(sqdists);
                m = ceil(n * obj.bicriteriaRatio);
                % f) store only the data for the closest points
                %    note that sqdists are already sorted so no need to use mIndexes
                PQS = PQS(mIndexes(1:m), :);
                sqdists = sqdists(1:m);

                %%%%%%%%%%%%%%%%%%%%
                % 5. Remove Q from P
                %%%%%%%%%%%%%%%%%%%%
                
                P(mIndexes(1:m), :) = [];
                n = size(P, 1);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 3. Select k / epsilon^2 points with probability
                %    proportional to the squared distance from S
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % This selects from the sorted list
                ind = randsample(m, obj.sampleSize, true, sqdists);
                
                % Count duplicates
                ind = sort(ind);
                h = histc(ind, ind);
                ind = ind(h > 0);
                h(h <= 0) = []; % how many times each index appeared
                mIndexes = mIndexes(ind);

                % Square root of the weights, for L2 norm we can apply wieghts by
                % multiplying the selected vectors by this
                % weight = 1/(c*pr)
                % pr = sqdists / sum(sqdists)
                % c = obj.sampleSize
                % remember to multiply by multiplicity (h)

                sqrtw = sqrt(sum(sqdists) / obj.sampleSize ./ sqdists(ind) .* h);

                % We need to take the part of the vectors perpendicular to S and not
                % the vectors from Q itself
                C1 = bsxfun(@times, sqrtw, OPQS(mIndexes, :));

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 4. Select another k / epsilon^2 points from S with probability
                %    proportional to the rows of U from the SVD decomposition
                %    or to the rows of Q from the QR decomposition
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
                % Sample based on SVD
%                 [U, ~, ~] = svd(PQS, 0);
                % Sample based on QR
                [U, ~] = qr(PQS, 0);
                % Sample based on SVD by using PCA - less stable but should be
                % good enough
%                 [V, D] = eig(PQS'*PQS);
%                 U = PQS*bsxfun(@times, V, 1./sqrt(diag(D)));

                sqdists = sum(U .* U, 2);

                ind = randsample(m, obj.sampleSize, true, sqdists);
                
                % Count duplicates
                ind = sort(ind);
                h = histc(ind, ind);
                ind = ind(h > 0);
                h(h <= 0) = []; % how many times each index appeared

                sqrtw = sqrt(sum(sqdists) / obj.sampleSize ./ sqdists(ind) .* h);

                % We take the part of the vectors projected onto the subspace
                C2 = bsxfun(@times, sqrtw, PQS(ind, :) * QS');

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Add selected elements into coreset
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                coreset = [coreset ; C1 ; C2];
            end
        end
        
        % Compute the optimal cost of the set P
        function optCost = computeOptCost(obj, ~)
            error('not implemented');
        end
        
        % copmute the cost of the center Q to the set P
        function cost = computeCost(obj, ~, ~)
            error('not implemented');
        end
        
        % Compute the optimal solution of P
        function opt = computeOpt(obj, ~)
            error('not implemented');
        end
    end
    
end

