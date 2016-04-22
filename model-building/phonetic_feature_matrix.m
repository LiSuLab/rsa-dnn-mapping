function FEATURES = phonetic_feature_matrix()

    % Feature vectors
    % Broad categories       aa ae ah ao aw ay b  ch d  ea eh er ey f  g  hh ia ih iy jh k  l  m  n  ng oh ow oy p  r  s  sh t  th ua uh uw v  w  y  z
    FEATURES.SONORANT    = [ 1  1  1  1  1  1  0  0  0  1  1  1  1  0  0  0  1  1  1  0  0  1  1  1  1  1  1  1  0  1  0  0  0  0  1  1  1  0  1  1  0 ];
    FEATURES.VOICED      = [ 1  1  1  1  1  1  1  0  1  1  1  1  1  0  1  0  1  1  1  1  0  1  1  1  1  1  1  1  0  1  0  0  0  0  1  1  1  1  1  1  1 ];
    FEATURES.SYLLABIC    = [ 1  1  1  1  1  1  0  0  0  1  1  1  1  0  0  0  1  1  1  0  0  0  0  0  0  1  1  1  0  0  0  0  0  0  1  1  1  0  0  0  0 ];
    FEATURES.OBSTRUENT   = [ 0  0  0  0  0  0  1  1  1  0  0  0  0  1  1  1  0  0  0  1  1  1  1  1  1  0  0  0  1  1  1  1  1  1  0  0  0  1  0  0  1 ];
    % Place of articulation for obstruents
    FEATURES.LABIAL      = [ 0  0  0  0  0  0  1  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  1  0  0  0  0  0  1  0  0  0  0  0  0  0  0  1  0  0  0 ];
    FEATURES.CORONAL     = [ 0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  1  0  1  0  1  0  0  0  0  0  1  1  1  1  1  0  0  0  0  0  1  1 ];
    FEATURES.DORSAL      = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0 ];
    % Manner of articulation for obstruents
    FEATURES.STOP        = [ 0  0  0  0  0  0  1  0  1  0  0  0  0  0  1  0  0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0 ];
    FEATURES.AFFRICATE   = [ 0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ];
    FEATURES.FRICATIVE   = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  1  0  0  0  1  0  0  1 ];
    FEATURES.SIBILANT    = [ 0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0  1 ];
    FEATURES.APPROXIMANT = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  1  1  0 ];
    FEATURES.NASAL       = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ];
    % Vowel backness
    FEATURES.FRONT       = [ 0  1  0  0  0  1  0  0  0  1  1  0  1  0  0  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ];
    FEATURES.CENTRAL     = [ 0  0  0  0  1  1  0  0  0  1  0  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  1  0  0  0  0  0  0 ];
    FEATURES.BACK        = [ 1  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  0  0  0  0  0  0  1  1  1  0  0  0  0 ];
    % Vowel closeness/height
    FEATURES.CLOSE       = [ 0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  0  0  0  0 ];
    FEATURES.CLOSEMID    = [ 0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  1  0  0  0  0  0  0 ];
    FEATURES.OPENMID     = [ 0  1  1  1  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0 ];
    FEATURES.OPEN        = [ 1  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ];
    % Vowel roundedness
    FEATURES.ROUNDED     = [ 0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  0  0  0  0  0  0  1  1  1  0  0  0  0 ];

end

