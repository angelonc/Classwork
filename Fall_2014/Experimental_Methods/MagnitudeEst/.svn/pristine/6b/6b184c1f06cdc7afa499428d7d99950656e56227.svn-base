function win = OpenEXPWindowME(bgRGB)
% win = OpenEXPWindowME(bgRGB)
%
% Open the experimental window for the magnitude estimation
% experiment, and set the background to the specified RGB
% values.

% Choose the last attached screen as our target screen, and figure out its
% screen dimensions in pixels.
d = mglDescribeDisplays;
screenDims = d(end).screenSizePixel;

% Open the window.
win = GLWindow('SceneDimensions', screenDims, 'BackgroundColor', bgRGB);
win.open;
win.draw;
