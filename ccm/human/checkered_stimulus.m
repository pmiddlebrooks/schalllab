function checkered_stimulus(stimulusEccentricity, stimulusAngle)

% stimulusEccentricity: in some units (e.g. pixels), the distance from the
%       center of the screen to the center of the stimulus
% stimulusAngle: in degrees, the angle from the center of the screen to the
%       center of the stimulus

% Format (numbering) of the squares in the stimulus (a 3 X 3 example):
%             0   1   2
%             3   4   5
%             6   7   8

% Define constants
iterations = 40;
stimulusColumns = 10;
stimulusRows = 10;
nSquares = stimulusColumns * stimulusRows;
squareWidth = 10;
% proportionGreens = [0 1];
proportionGreens = [.38 .42 .46 .54 .58 .62];

% Initialize variables
testCollection = zeros(nSquares, 2);
testCollectionUL = zeros(nSquares, 2);
testCollectionLR = zeros(nSquares, 2);
% testCollectionXs = zeros(stimulusColumns);
% testCollectionYs = zeros(stimulusRows);

% define the units you're using:
scrsz = get(0, 'ScreenSize');
screenWidth = scrsz(3);
screenHeight = scrsz(4);
% First determine the center of the 10X10 stimulus in x,y coordinates w.r.t the center of the
% screen
stimulusOffsetX = stimulusEccentricity * cosd(stimulusAngle);
stimulusOffsetY = stimulusEccentricity * sind(stimulusAngle);

figure('position', [0 0 screenHeight screenHeight]);
set(gcf, 'color', [1 1 1])
xlim([-screenHeight/2 screenHeight/2]);
ylim([-screenHeight/2 screenHeight/2]);
hold all;
for trial = 1 : iterations
    randomizedGreens = int16(randperm(length(proportionGreens)));
    randomGreenProportion = randomizedGreens(1);
    nGreen = int16(proportionGreens(randomGreenProportion) * nSquares);
    randomGreenIndices = randperm(nSquares);
    randomGreen = randomGreenIndices(1 : nGreen);
    colorArray = zeros(nSquares, 1);
    colorArray(randomGreen) = 1;
    for y = 1 : stimulusRows
        for x = 1 : stimulusColumns
            stimulusIndex = x + (stimulusColumns * (y - 1));
            % Get the centers of each square
            iSquareCenterX = stimulusOffsetX + squareWidth * ((x-1) - (stimulusColumns/2 - .5));
            iSquareCenterY = stimulusOffsetY + squareWidth * ((stimulusRows/2 - .5) - (y-1));
            iSquareEccentricity = sqrt(iSquareCenterX^2 + iSquareCenterY^2);
            iSquareAngle = atand(iSquareCenterY / iSquareCenterX);
            
            % Get the upper left and lower right corners of each square
            iSquareULX = iSquareCenterX - squareWidth/2;
            iSquareULY = iSquareCenterY + squareWidth/2;
            iSquareLRX = iSquareCenterX + squareWidth/2;
            iSquareLRY = iSquareCenterY - squareWidth/2;
            
            if colorArray(stimulusIndex) == 0
                iSquareColor = 'r';
            else
                iSquareColor = 'b';
            end
            rectangle('position', [iSquareCenterX - squareWidth/2, iSquareCenterY - squareWidth/2, squareWidth, squareWidth], ...
                'facecolor', iSquareColor, 'edgecolor', iSquareColor)
            % code from tempo:
            % 	ulx       = round((stim_ecc_x - half_size)*conversion_X);
            % 	uly       = round((stim_ecc_y + half_size)*conversion_Y);
            % 	lrx       = round((stim_ecc_x + half_size)*conversion_X);
            % 	lry       = round((stim_ecc_y - half_size)*conversion_Y);
            testCollection(stimulusIndex, :) = [iSquareCenterX iSquareCenterY];
            testCollectionUL(stimulusIndex, :) = [iSquareULX iSquareULY];
            testCollectionLR(stimulusIndex, :) = [iSquareLRX iSquareLRY];
            % testCollectionXs(x, y) = iSquareCenterX;
            % testCollectionYs(x, y) = iSquareCenterY;
            
            
        end
    end
    pause
end
scatter(stimulusOffsetX, stimulusOffsetY, '*', 'k')
hold on;
% scatter(testCollection(:, 1), testCollection(:, 2), 'o', 'k')
% scatter(testCollectionUL(:, 1), testCollectionUL(:, 2), '.', 'r')
% scatter(testCollectionLR(:, 1), testCollectionLR(:, 2), '.', 'b')
% testCollectionXs
% testCollectionYs
% scatter(testCollectionXs(:), testCollectionYs(:))
