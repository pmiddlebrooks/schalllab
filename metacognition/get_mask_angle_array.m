function targetAngleArray = get_mask_angle_array(angles, nTarget, maskDistribution)


targetAngleArray = nan(1, nTarget);

switch maskDistribution
    case 'mirror'
        targetAngleArray(1, 1) = angles(1);
        targetAngleArray(1, 2) = angles(2);
        targetAngleArray(1, 3) = 180 - angles(1);
        targetAngleArray(1, 4) = 180 - angles(2);
    case 'maximize'
        targetAngleArray(1, 1) = angles(1);
        degreesPerMask = 360 / nTarget;
        for iTarget = 2 : nTarget
            targetAngleArray(1, iTarget) = targetAngleArray(1, iTarget-1) + degreesPerMask;
        end
end
