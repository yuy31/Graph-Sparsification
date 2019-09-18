% read edge list from file and construct a graph
input = fopen('facebook_combined.txt');
C = textscan(input, '%s%s');
fclose(input);

G = graph(C{1}, C{2});
n = numnodes(G);
m = numedges(G);

% Begin the computation of Effective Resistances
B_trans = incidence(G);
L = laplacian(G);

B_trans = full(B_trans);
L = full(L);

% B is the incidence matrix of the graph
B = B_trans';
L_inv = pinv(L);

% compute the k*n transformation matrix Q where k = 24log/epsilon^2 and
% compute the approximations of effective resistances
%{
e = 0.5;
k = 24 * log(m) / (e^2);
k = floor(k);
Q = [];
M = [];

for i = 1:k
    for j = 1:m;
        M(i, j) = randn;
    end
end

Q = 1/(sqrt(k)) * M;
Z = Q * B * L_inv;
r = [];

for i = 1:m
    v = B(i, :);
    v = v';
    v = Z * v;
    r(i) = norm(v);
end
r = r';
%}

% R(e, e) now is the effective resistance of each edge
R = zeros(1,m);
for i = 1:m
    R(i) = B(i,:) * L_inv * B(i,:)';
end
%R = B * L_inv * B_trans;

% begin sampling according to effective resistances
P = R';

%clear B;
clear B_trans;
clear L_inv;
clear R;

% number of samples
% QUESTION: how many samples to take?
k = 50000;

% 'test' contains all the selected edges
test = randsample(1:m,k,true,P);

% get the corresponding vertices according to the selected edges and the
% incidence matrix, convert into a new adjacency matrix and include the
% updated weights
new_A = [];
for i = 1:n
    for j = 1:n
        new_A{i, j} = 0;
    end
end

% for each selected edge, find the vertices it is incident on and calculate
% its new weight for the new adjacency matrix
for i = 1:length(test)
    row = test(i);
    v = B(row,:);
    for j = 1:length(v)
        if v(j) == 1
            head = j;
        end
        if v(j) == -1
            tail = j;
        end
    end
    % calculate the new weight of the edge
    new_A{head, tail} = new_A{head, tail} + 1/(k*P(row));
    new_A{tail, head} = new_A{tail, head} + 1/(k*P(row));
end

X = cell2mat(new_A);
tf = issymmetric(X);

% H is the approximation graph
H = graph(X);
%LWidths = 5 * H.Edges.Weight/max(H.Edges.Weight);
%t = plot(H, 'NodeLabel', nLabels, 'LineWidth', LWidths);

% calculate the quadratic forms and epsilon
L_new = laplacian(H);
L_new = full(L_new);
% generate a random vector for computing quadratic forms
x = rand(n, 1);
Q_G = x' * L * x;
Q_H = x' * L_new * x;
% calculate epsilon to see how good the approximation is
epsilon = (Q_G - Q_H) / Q_G;