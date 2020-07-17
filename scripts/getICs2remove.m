function [ sublist, EEGICs2remove, MEGICs2remove ] = getICs2remove()
%GETICS2REMOVE returns the list of subjects and the EEG/MEG ICs to remove
%for each. Used at multiple points in the analysis.

% Commented out S08 here, as this is the easiest way to drop a subject.

% Independent components to remove from EEG/MEG: (previously determined)
EEGICs2remove = {
    'S01'   [ 1 2 3       ];
    'S02'   [ 1 2 3 4     ];
    'S03'   [ 1 2 3 4     ];
    'S04'   [ 1 2 3 4 9   ];
    'S05'   [ 1 2 3 4 5 8 ];
    'S06'   [ 4 5 6 8     ];
    'S07'   [ 1 2 3       ]; 
%     'S08'   [ 1 3 4 6 14  ];
    'S09'   [ 2 3 4 5     ];
    'S10'   [1 2 3 4 5 6 7 8];
};

MEGICs2remove = {
    'S01'   [ 1 ];
    'S02'   [ 1 7 ];
    'S03'   [ 3 ];
    'S04'   [ 7 ];
    'S05'   [ 2 ];
    'S06'   [ 6 ];
    'S07'   [ 2 ];
%     'S08'   [ 2 13 ];
    'S09'   [ 2 6 ];
    'S10'   [ 6 ];
};

assert(isequal(EEGICs2remove(:,1), MEGICs2remove(:,1)));
sublist = EEGICs2remove(:,1);

end

