% Optimal SOR results
dim = [5 6 11 21 31];
dim2 = dim.*dim;
w = [1.32 1.45 1.65 1.806 1.866];

% dim = [6 11 21 31];
% dim2 = dim.*dim;
% w = [1.45 1.65 1.806 1.866];

x = dim2;                                                    % Independent Variable
y = w;                                                    % Dependent Variable
fcn1 = @(b,x) b(1).*exp(b(2).*x.^b(3)) + b(4);                       % Objective Function #1
fcn2 = @(b,x) b(1).*exp(b(2).*x.^b(3));                              % Objective Function #2
SSECF = @(b) sum((y - fcn1(b,x)).^2);                       % Sum-Squared-Error Cost Function (Use ‘fcn2’ Here)
B0 = [1; 1];                                                % Initial Parameter Estimates
options = optimset('MaxFunEvals',100000);
[B,SSE] = fminsearch(SSECF, [1; -0.2; -1; 1],options);                        % Estimate Parameters
xv = linspace(min(x), max(x));
figure(1)
plot(x, y, 'bp')
hold on
plot(xv, fcn1(B,xv), '-r')
%hold off
grid

title('Optimal w for SOR')
xlabel('Dimension of Laplacian Matrix')
ylabel('w')


xa = 5:1:40000;
ya = 13.4523570058092 * exp(-0.206450260650164 * xa.^-0.434163866503769) - 11.4497834449085;
plot(xa,ya)

%
% Final result
%
%  w = 13.4523570058092 * exp(-0.206450260650164 * dim^-0.434163866503769) - 11.4497834449085
%

