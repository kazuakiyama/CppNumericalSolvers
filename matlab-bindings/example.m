% we minimize the rosenbrock function
%clc
checkrosenbrock = @(s) disp(['|err|=' num2str(norm(s-[1;1]))]);

fprintf('\n------------------------------------------------------------\n\n');
fprintf(' conj. grad-descent starting in [-1.2; 1] \n');
tic
solution = cppsolver([-1.2; 1],@rosenbrock,'gradient',@rosenbrock_grad,'solver','cg');
checkrosenbrock(solution);
toc

fprintf('\n------------------------------------------------------------\n\n');
fprintf(' BFGS starting in [-1.2; 1] \n');
tic
solution = cppsolver([-1.2; 1],@rosenbrock,'gradient',@rosenbrock_grad,'solver','bfgs');
checkrosenbrock(solution);
toc

fprintf('\n------------------------------------------------------------\n\n');
fprintf(' L-BFGS starting in [-1.2; 1]\n');
tic
solution = cppsolver([-1.2; 1],@rosenbrock,'gradient',@rosenbrock_grad);
checkrosenbrock(solution);
toc

fprintf('\n------------------------------------------------------------\n\n');
fprintf(' L-BFGS starting in [-1.2; 1] (without gradient, this is a WARNING!)\n');
tic
solution = cppsolver([-1.2; 1],@rosenbrock);
checkrosenbrock(solution);
toc

fprintf('\n------------------------------------------------------------\n\n');
fprintf(' NEWTON with given hessian and gradient starting in [-1.2; 1]\n');
tic
solution = cppsolver([-1.2; 1],@rosenbrock,'gradient',@rosenbrock_grad,'hessian',@rosenbrock_hessian,'solver','newton');
checkrosenbrock(solution);
toc

fprintf('\n------------------------------------------------------------\n\n');
fprintf(' fminsearch starting in [-1.2; 1]\n');
tic
solution = fminsearch(@rosenbrock,[-1.2; 1]);
checkrosenbrock(solution);
toc

fprintf('\n------------------------------------------------------------\n\n');
fprintf(' fminunc starting in [-1.2; 1] \n');
tic
solution = fminunc(@rosenbrock,[-1.2; 1]);
checkrosenbrock(solution);
toc

