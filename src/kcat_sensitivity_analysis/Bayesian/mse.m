function err = mse(y,x)
% mse  Mean squared error between two arrays.
%
% Computes the mean squared error between y and x, defined as the squared
% Euclidean norm of their difference divided by the number of elements.
%
% Parameters
% ----------
% y : double
%     first array of values.
% x : double
%     second array of values, the same size as y.
%
% Returns
% -------
% err : double
%     the mean squared error between y and x.
err = (norm(x(:)-y(:),2).^2)/numel(x);
end