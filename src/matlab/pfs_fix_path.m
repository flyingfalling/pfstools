function [fixed_path, cygwin_home] = pfs_fix_path( path )
% Used to fix paths when running cygwin on Windows. This is an internal
% command. Do not use.

if ~ispc()
    fixed_path = path;
else
    
    pstatus = dos('set CYGWIN_HOME > nul 2>&1');
    if(pstatus == 1)
        check_dirs  = { 'c:\\cygwin', 'c:\\cygwin64', 'c:\\cygwin32' };
        cygwin_home = [];
        for kk=1:length(check_dirs)
            if exist( check_dirs{kk}, 'dir' )
                cygwin_home = check_dirs{kk};
                break;
            end
        end
        if isempty( cygwin_home )
            error( 'Cannot find Cygwin home directory. Set CYGWIN_HOME environment variable.' );
        end
    else
        [pstatus, cygwin_home] = dos('echo %CYGWIN_HOME%');
        cygwin_home = strtrim(cygwin_home); % to remove the final LF
    end
    
    [pstatus, fixed_path] = dos( [ fullfile( cygwin_home, 'bin', 'cygpath' ) ' ' path] );
    
    fixed_path = strtrim( fixed_path );
end

end