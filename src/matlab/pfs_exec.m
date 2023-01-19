function [status, cmdout] = pfs_exec( cmd_line )
% Execute commnand line in the native environment 

if ispc()
    [work_dir, cygwin_home] = pfs_fix_path( pwd() );
            
    % Replace all DOS paths with Unix paths for cygwin
    cl = strsplit( cmd_line );
    cmd_line = [];
    for kk=1:length(cl)
        if( ~isempty(strfind(cl{kk}, '\')) )
            % If a string contains \, it is assumed to be a path
            cmd_line = cat( 2, cmd_line, ' ', pfs_fix_path( cl{kk} ) );
        else
            cmd_line = cat( 2, cmd_line, ' ', cl{kk} );
        end
    end
    
    % This is put at the beginning of the shell command
    cmd = [ cygwin_home '\bin\bash -l -c "cd ''' work_dir ''';' cmd_line '"'];           
else
    work_dir = pwd();

    % It is necessary to set all ENV variables before invoking
    % pfstools commands. '/bin/bash' may need to be replaced with the
    % shell you are using.
        
    % This will remove all references to matlab libraries from the
    % LD_LIBRARY_PATH. pfstools usually do not work with matlab version
    % of the standard libraries
    set_ld_path='export LD_LIBRARY_PATH=`echo $LD_LIBRARY_PATH | sed "y/:/\n/" | grep -v "matlab" | sed ":beg;N;s/\n/:/;t beg;"`';
        
    if( strcmp( computer, 'MACI64' ) || strcmp( computer, 'MACI32' ) )
        set_ld_path=cat( 2, set_ld_path, '; export DYLD_FRAMEWORK_PATH=`echo $DYLD_FRAMEWORK_PATH | sed "y/:/\n/" | grep -v "matlab" | sed ":beg;N;s/\n/:/;t beg;"`' );
    end
        
    % Start shell, cd to the working directory and remove matlab
    % libraries
    cmd = [ '/bin/bash -l -c "cd ''' work_dir '''; ' set_ld_path '; ' cmd_line '"' ];
end

if nargout == 2
    [status, cmdout] = system( cmd );
else
    status = system( cmd );
end

end
