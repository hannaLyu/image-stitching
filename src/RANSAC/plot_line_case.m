
function plot_line_case(X,result)
    figure;
    hold on
    % ind = results.CS;
    plot(X(1, result.inliers), X(2, result.inliers), '.g')
    plot(X(1, ~result.inliers), X(2, ~result.inliers), '.r')

    xmin = min(X(1,:));
    xmax = max(X(1,:));
    xx = linspace(xmin,xmax,100);
    yy = -(result.params(1).*xx+result.params(3))./(result.params(2)+1e-6);
    plot(xx, yy, 'm--','LineWidth',2);
    xlabel('x')
    ylabel('y')
    title('RANSAC results for 2D line estimation')
    legend('Inliers', 'Outliers')
    axis equal tight
end