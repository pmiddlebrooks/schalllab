function [frasters] = fat_raster( rasters, howfat )

% [frasters] = fat_raster( rasters, howfat )
%
%  Take a raster line (a single raster line), or a list of rasters where each
%  row is a single raster, and make each '1' into
%  a series of '1's, to make it fatter and easier to see with something
%  like imagesc.

frasters = rasters;
sz = size( frasters );

for rw = 1:sz(1)
    raster = frasters( rw,: );
    f = find( raster > 0 );
    newf = f;
    lf = length( f );
    lr = length( raster );
    if ~isempty( f )
        value = raster( f( 1 ) );
        fraster = raster;
        for d = 2:howfat
            newf = newf + 1;
            okf = find( newf < lr );
            fraster( newf( okf ) ) = value;
        end;
    else
        fraster = raster;
    end;
    
    frasters( rw, : ) = fraster;
end;

        