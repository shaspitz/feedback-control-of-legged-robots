% Animate the Three-Link Walker
function animateThreeLink(tData, qData)
    figure(1000)

    for i =1:length(tData) 
        clf ;
        drawThreeLink(qData(i, :)');
        line([-1, 20],[0;0],'Color', 'k', 'LineWidth', 2)
        axis([-1 15 -1 5]) ; 
        grid on ;
        drawnow ;
        pause(0.001) ;
    end
end


% Draw one frame of the Three-Link Walker
function drawThreeLink(q)
    x = q(1);
    y = q(2);
    q1 = q(3);
    q2 = q(4);
    q3 = q(5);

    pSw = pSw_gen([q;zeros(5,1)]);
    pSt = pSt_gen([q;zeros(5,1)]);
    pT = pComTorso_gen([q;zeros(5,1)]);
    l1 = line([x;pT(1)], [y;pT(2)], 'Color', 'k', 'LineWidth', 2);
    hold on
    l2 = line([x;pSw(1)], [y;pSw(2)], 'Color', 'b', 'LineWidth', 2);
    l3 = line([x;pSt(1)], [y;pSt(2)], 'Color', 'r', 'LineWidth', 2);
    plot(pSw(1), pSw(2), 'bo', 'MarkerSize',7,'MarkerEdgeColor','b','MarkerFaceColor','g')
    plot(pSt(1), pSt(2), 'ro', 'MarkerSize',7,'MarkerEdgeColor','r','MarkerFaceColor','g')
    plot(pT(1), pT(2), 'ko', 'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','g')
    plot(x, y, 'ko', 'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','g')
end