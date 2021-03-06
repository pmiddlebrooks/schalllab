function status = om_check_vol(vol)
%   OM_CHECK_VOL   Check meshes of volume conductor for BEM modeling
%       [STATUS] = OM_CHECK_VOL(VOL)
%
%   returns 1 if there is a problem with geometry
%   else returns 0
%
% Copyright (C) 2010, Alexandre Gramfort, INRIA

% $Id: om_check_vol.m 2212 2010-11-27 11:55:07Z roboos $
% $LastChangedBy: alegra $
% $LastChangedDate: 2010-09-06 13:58:49 +0200 (Mon, 06 Sep 2010) $
% $Revision: 2212 $

openmeeg_license
om_checkombin;

% the first compartment should be the skin, the last the source
% flip the order of the compartments if necessary
if vol.skin==length(vol.bnd) && vol.source==1
    vol.bnd    = fliplr(vol.bnd(:)');
    vol.skin   = 1;
    vol.source = length(vol.bnd);
end

assert(vol.skin == 1)
assert(vol.source == length(vol.bnd))

% Flip faces for openmeeg convention
for ii=1:length(vol.bnd)
    vol.bnd(ii).tri = fliplr(vol.bnd(ii).tri);
end

try
    % store the current path and change folder to the temporary one
    tmpfolder = cd;

    cd(tempdir)

    % write the triangulations to file
    bndfile = {};
    for ii=1:length(vol.bnd)
        [junk,tname] = fileparts(tempname);
        bndfile{ii} = [tname '.tri'];
        om_save_tri(bndfile{ii}, vol.bnd(ii).pnt, vol.bnd(ii).tri);
    end

    % these will hold the shell script and the inverted system matrix
    [junk,tname] = fileparts(tempname);
    if ~ispc
      exefile = [tname '.sh'];
    else
      exefile = [tname '.bat'];
    end

    [junk,tname] = fileparts(tempname);
    geomfile  = [tname '.geom'];

    % write conductivity and geometry files
    om_write_geom(geomfile,bndfile);

    % Exe file
    status = system(['om_check_geom -g ',geomfile]);
    cleaner(vol,bndfile,geomfile)
    cd(tmpfolder)
catch
    cleaner(vol,bndfile,geomfile)
    cd(tmpfolder)
    rethrow(lasterror)
end

function cleaner(vol,bndfile,geomfile)

% delete the temporary files
for i=1:length(vol.bnd)
    delete(bndfile{i})
end

delete(geomfile);
return
