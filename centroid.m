function [comX,comY] = centroid(A)
x = 1:size(A,2); 
y = 1:size(A,1); 
[X, Y] = meshgrid(x,y);
meanA = mean(A(:));
comX = mean(A(:).*X(:))/meanA;
comY = mean(A(:).*Y(:))/meanA;
end

