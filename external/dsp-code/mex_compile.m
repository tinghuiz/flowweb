if ( ispc )
    mex -f ./mexopts.bat SPNodeMatchMex.cpp SPNodeMatch.cpp DT.cpp
    mex -f ./mexopts.bat SPBPMex.cpp SPBP.cpp SPGraph.cpp DT.cpp
    mex -f ./mexopts.bat PixelMatchMex.cpp PixelMatch.cpp 
elseif ( isunix )
    mex -f ./gccopts.sh SPNodeMatchMex.cpp SPNodeMatch.cpp DT.cpp
    mex -f ./gccopts.sh SPBPMex.cpp SPBP.cpp SPGraph.cpp DT.cpp
    mex -f ./gccopts.sh PixelMatchMex.cpp PixelMatch.cpp 
else
    disp('NOT SUPPORTED');
end
    