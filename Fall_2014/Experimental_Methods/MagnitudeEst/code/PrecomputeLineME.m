function lineInfo = PrecomputeLineME(stimLength, stimAngle, arrowLength, arrowAngle, horizontalOffset)
% lineInfo = PrecomputeLineME(stimLength, stimAngle, arrowLength, arrowAngle, [horizontalOffset])
%
% Description:
% Precompute all of the icky little line coordinates and put them into a
% structure for us to use when we actually draw the lines. 
%
% Output:
% lineInfo (struct) - Struct containing the (x,y) pairs for the line and arrows.

% Optional horizontal offset, could be used to jitter line positions
if nargin < 5 || isempty(horizontalOffset)
	horizontalOffset = 0;
end

% Calculate line from and to coordinates depending length and angle.
hLength = (stimLength/2)*cos(2*pi*stimAngle/360);
vLength = (stimLength/2)*sin(2*pi*stimAngle/360);
fromH = -hLength;
fromV = -vLength;
toH = hLength;
toV = vLength;

% Now compute arrows.  Angle of zero is special case
% of no arrow, independent of length.
if arrowAngle == 0
	aFromUpH = fromH;
	aFromUpV = fromV;
	aFromDwnH = fromH;
	aFromDwnV = fromV;
	aToUpH = toH;
	aToUpV = toV;
	aToDwnH = toH;
	aToDwnV = toV;
else
	fromUpAngle = stimAngle + arrowAngle;
	aFromUpH = fromH + arrowLength*cos(2*pi*fromUpAngle/360);
	aFromUpV = fromV + arrowLength*sin(2*pi*fromUpAngle/360);
	
	fromDwnAngle = stimAngle - arrowAngle;
	aFromDwnH = fromH + arrowLength*cos(2*pi*fromDwnAngle/360);
	aFromDwnV = fromV + arrowLength*sin(2*pi*fromDwnAngle/360);
	
	toUpAngle = stimAngle + 180 - arrowAngle;
	aToUpH = toH + arrowLength*cos(2*pi*toUpAngle/360);
	aToUpV = toV + arrowLength*sin(2*pi*toUpAngle/360);
	
	toDwnAngle = stimAngle  + 180 + arrowAngle;
	aToDwnH = toH + arrowLength*cos(2*pi*toDwnAngle/360);
	aToDwnV = toV + arrowLength*sin(2*pi*toDwnAngle/360);
end

% Stuff point info for return
lineInfo.from = [fromH+horizontalOffset, fromV];
lineInfo.to = [toH+horizontalOffset, toV];
lineInfo.upArrowFrom = [aFromUpH+horizontalOffset, aFromUpV];
lineInfo.upArrowTo = [aToUpH+horizontalOffset, aToUpV];
lineInfo.downArrowFrom = [aFromDwnH+horizontalOffset, aFromDwnV];
lineInfo.downArrowTo = [aToDwnH+horizontalOffset, aToDwnV];
