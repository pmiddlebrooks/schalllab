function inWindow = inAcceptWindow(eyeX, eyeY, acceptWindow)
%
% determine if gx and gy are within a defined window

inWindow = eyeX > acceptWindow(1) &&  eyeX <  acceptWindow(3) && ...
    eyeY > acceptWindow(2) && eyeY < acceptWindow(4) ;
end


