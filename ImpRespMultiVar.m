% System Identificataion Algorithm Based on the Multivariable Impulse
% Response
% Author: Luis Fernando Segalla
% Based on the lecture from professor Mateus Giesbrecht

% Input: y = Impulse response (can be single ou multi variable)
% Output A,B,C,D = Estimated State Space matrixes 
% 
function [A, B, C, D] = ImpRespMultiVar(y,gamma)
    %Start by checking if the input is valid
    if(size(y) == 0)
        error('Invalid data. Input must be not empty.');
    end
    if(gamma == 0)
        error('Invalid data. Gamma must be positive greater than zero');
    end
    
    % The input y is a stack of LxM matrixes
    % containing gamma elements. These elements will be used to create the
    % block Hankel matrix o syze gamma*L x gamma*M. 
    % As the elements are already stacked in y the size L already is
    % gamma*L (or gamma times the lines of a single impulse response
    % matrix)
    [l,m] = size(y);
    if(m>l)
        disp('Transposed the Input Vector')
        [l,m] = size(y');
    end
    
    % block Hankel matrix of the impulse response 
    %H = blkhank(y,l,m*gamma);
    H = blkhank(y,10,10);
    disp(size(H))
    % Singular Value Decomposition of the H matrix
    [U,S,V] = svd(H);
    
    % Using only the non-null singular values of H
    S_sing = diag(S);
    % Number of non-null singular values of H
    n = size(S_sing);
    %disp(S_sing)
    
    S_singMatrix = zeros(n(1),n(1));
    U1 = zeros(l,n(1));
    %size(U1)
    %size(U)
    V1 = zeros(n(1),m*gamma);
    %size(V1)
    for i=1:n(1)
        S_singMatrix(i,i) = S_sing(i);
        U1(:,i)           = U(:,i); 
        V1(i,:)           = V(i,:);
    end
    
    % Defining the observability matrix
    O = U1*sqrt(S_singMatrix);
    % Defining the controllability matrix
    C = sqrt(S_singMatrix)*V1;
    
    % Defining the A matrix
    A = pinv(O)*H*pinv(C);
    
    % Defining the B matrix
    H_col = H(:,1:m);
    B = pinv(O)*H_col;
    
    % Defining the C matrix
    H_lin = H(1:(l/gamma),:);
    C = H_lin*pinv(C);
    
    % Defining the D matrix
    D = H(1:(l/gamma),1:m);
    
end