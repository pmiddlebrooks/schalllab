function betAngleArray = get_bet_angle_array(betAngle, betDistribution)


betAngleArray = nan(1, 2);

switch betDistribution
    case 'mirror'
        betAngleArray(1, 1) = betAngle;
        betAngleArray(1, 2) = 180 - betAngle;
    case 'maximize'
        betAngleArray(1, 1) = betAngle;
        betAngleArray(1, 2) = betAngle + 180;
end
