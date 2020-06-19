function eventVideoBehav = eventVideo(All, event, neuron, bef, aft,varargin)

% Returns struct with analysed video data per indicated event type
% (correct, incorrect, etc.)

% Input:
% - eventChannelNumber = 6, TTL channel number of event (correct = 6,
% Incorrect = 8, Omission = 10, Premature = 12, Perserverative = 13)
% - bef = 6, time in seconds before event
% - aft = 0, time in seconds after event
% - firstCuelightFrameNr = 12080, Indicate visually inspected first frame
% number where a cue light went on. Needed for video/ephys offset correction.
%
% ouput:
% - eventVideoBehav, Struct with x-y coordinates, time, speed, time between event and treshold crossing, treshold value in meter from border above cue light and
% firstCuelightFrameNr used for offset correction.
%
%
% V 1.0
% 31-12-2017
% h.terra@vu.nl
%

% Set fixed variables: Framerate and Cuelight treshold
frameRate = 25;
treshold = All{12, 1}(neuron).length*3/8;

% Get event times
eventsAll = All{22,1}{neuron,1};
%CH = sprintf('CH%d', eventChannelNumber);
if strcmp(class(event),'char')
    events = eventsAll.(event);
    display('video for PSTH')
else
    events = event;
    display('video for speed')
end

if nargin
for iEvent=1:size(varargin,2)
    switch varargin{iEvent}
        case 'PSTH_EVENT'
            if strcmp(event, 'corITI125Cue') || strcmp(event, 'incITI125Cue') || strcmp(event, 'omITI125Cue')
                events = eventsAll.(event)-12.5;
            elseif strcmp(event, 'corITI75Cue') || strcmp(event, 'incITI75Cue') || strcmp(event, 'omITI75Cue')
                events = eventsAll.(event)-7.5;
            elseif strcmp(event, 'corITI50Cue') || strcmp(event, 'incITI50Cue') || strcmp(event, 'omITI50Cue')
                events = eventsAll.(event)-5;
            elseif strcmp(event, 'premResp')
                events = eventsAll.(event)-eventsAll.preLatency;
            elseif strcmp(event, 'premITI50Resp')
                events = eventsAll.(event)-eventsAll.preLatencyITI5;
            elseif strcmp(event, 'premITI75Resp')
                events = eventsAll.(event)-eventsAll.preLatencyITI75;
            elseif strcmp(event, 'premITI125Resp')
                events = eventsAll.(event)-eventsAll.preLatencyITI125;
            elseif strcmp(event, 'premSD2Resp')
                events = eventsAll.(event)-eventsAll.preLatencySD2;
            elseif strcmp(event, 'premSD5Resp')
                events = eventsAll.(event)-eventsAll.preLatencySD5;
            elseif strcmp(event, 'premSD10Resp')
                events = eventsAll.(event)-eventsAll.preLatencySD10;
            elseif strcmp(event, 'omCue') || strcmp(event, 'incCue') || strcmp(event, 'corCue') || strcmp(event, 'corSD10Cue') || strcmp(event, 'incSD10Cue') || strcmp(event, 'omSD10Cue') || strcmp(event, 'corSD5Cue') || strcmp(event, 'incSD5Cue') || strcmp(event, 'omSD5Cue') || strcmp(event, 'corSD2Cue') || strcmp(event, 'incSD2Cue') || strcmp(event, 'omSD2Cue')
                events = eventsAll.(event)-5;
            end
    end
end
end
                

% Convert seconds to frame numbers
bef = bef*frameRate;
aft = (aft*frameRate)-1;

% Get x, y, speed etc per event and store in eventVideoBehav struct.
for i = 1:size(events,2)
    % Check which index TTL time and video time match the most (they have a
    % small mismatch).
    [c, index] = min(abs(events(i)-All{12, 1}(neuron).time));

    try
        % Get x, y time and speed per event and store in struct
        eventVideoBehav.x(i,:) = All{12, 1}(neuron).x(index-bef:index+aft);

        % Check to replace any 0 entry errors by NaN before clalculating speed and
        % treshod
        zeroIx = eventVideoBehav.x(i,:) == 0;
        [r c] = find(zeroIx);
        eventVideoBehav.x(i,c) = nan;

        eventVideoBehav.y(i,:) = All{12, 1}(neuron).y(index-bef:index+aft);

        % Check to replace any 0 entry errors by NaN before clalculating speed and
        % treshod
        zeroIx = eventVideoBehav.y(i,:) == 0;
        [r c] = find(zeroIx);
        eventVideoBehav.y(i,c) = nan;

        eventVideoBehav.time(i,:) = All{12, 1}(neuron).time(index-bef:index+aft);
        eventVideoBehav.speed(i,:) = All{12, 1}(neuron).speed(index-bef:index+aft);
    catch % If initial indexing did not work replace missing values by NaN (usually the case with first event).
        indexL = bef+aft+1;
        eventVideoBehav.x(i,1:indexL) = nan;
        eventVideoBehav.y(i,1:indexL) = nan;
        eventVideoBehav.time(i,1:indexL) = nan;
        eventVideoBehav.speed(i,1:indexL) = nan;
%         indexR = bef-index+1;
%         eventVideoBehav.x(i,indexR+1:end) = All{12, 1}(neuron).x(index-bef+indexR:index+aft);
%         eventVideoBehav.y(i,indexR+1:end) = All{12, 1}(neuron).y(index-bef+indexR:index+aft);
%         eventVideoBehav.time(i,indexR+1:end) = All{12, 1}(neuron).time(index-bef+indexR:index+aft);
%         eventVideoBehav.speed(i,indexR+1:end) = All{12, 1}(neuron).speed(index-bef+indexR:index+aft);
    end

    % Calculate time (seconds) in takes from treshold crossing to event.

    if numel(find(eventVideoBehav.y(i,:)<treshold))>=1
        tresholdIx = find(eventVideoBehav.y(i,:)<treshold);
        eventVideoBehav.treshold(i,:) = tresholdIx(1)/25; %Time in seconds between passing treshold and behav event
        if  numel(find(eventVideoBehav.y(i,tresholdIx(1):end)>treshold))>=1
            eventVideoBehav.tresholdSecondPass(i,:) = 1; %Indicate if animal passed the treshold again after first passing
        else
            eventVideoBehav.tresholdSecondPass(i,:) = 0;
        end
        eventVideoBehav.tresholdLateralDist(i,:) = max(eventVideoBehav.x(i,tresholdIx(1):end)-eventVideoBehav.x(i,tresholdIx(1)))+...
        abs(min(eventVideoBehav.x(i,tresholdIx(1):end)-eventVideoBehav.x(i,tresholdIx(1))));
    elseif numel(find(eventVideoBehav.y(i,:)<treshold))<1;
        eventVideoBehav.treshold(i,:) = nan;
        eventVideoBehav.tresholdSecondPass(i,:) = nan;
        eventVideoBehav.tresholdLateralDist(i,:) = nan;
    end


end

% Set all trials with treshold value of 6 to NaN. These trials seem to be
% reversed. I don't know why. Maybe Bonsai bug
% six = eventVideoBehav.treshold == 6;
% eventVideoBehav.x(six,:) = nan;
% eventVideoBehav.y(six,:) = nan;
% eventVideoBehav.treshold(six,:) = nan;
% eventVideoBehav.speed(six,:) = nan;
% eventVideoBehav.time(six,:) = nan;
eventVideoBehav.tresholdMeter = treshold;
eventVideoBehav.firstCuelightFrameNr = All{12, 1}(neuron).firstCuelightFrameNr;