

a = rand(5000,5000,'gpuArray');
b = rand(5000,5000,'gpuArray');

% a = rand(5000,5000);
% b = rand(4000,4000);


parfor i = 1:10
    c = a/b;
    disp(i)
end