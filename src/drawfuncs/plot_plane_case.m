
function plot_plane_case(X,result,Xi,Xo)
    figure;
    hold on
    plot3(Xi(1, :), Xi(2, :), Xi(3, :), '+g','MarkerSize', 15)
    plot3(Xo(1, :), Xo(2, :), Xo(3, :), '+r','MarkerSize', 15)
    plot3(X(1, result.inliers), X(2, result.inliers), X(3, result.inliers), 'sg','MarkerSize', 15)
    plot3(X(1, ~result.inliers), X(2, ~result.inliers), X(3, ~result.inliers), 'sr','MarkerSize', 15)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('RANSAC results for 3D plane estimation')
    legend('Inliers', 'Outliers', 'Estimated Iniliers', 'Estimated Outliers')
    axis equal tight
    view(3)
end

