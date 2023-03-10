function cmd = pfs_shell( suffix )
% Returns command line that starts shell. Internal, do not use. 

if ispc()
    % -i option is needed to make sure that .bash_rc is executed and thus
    % DISPLAY and other environment variables are set
    if( exist( 'suffix', 'var' ) )
        % This is put at the end of the shell command
        cmd = '"';
    else
        [work_dir, cygwin_home] = pfs_fix_path( pwd() );
        
        % This is put at the beginning of the shell command
        cmd = sprintf('%s\\bin\\bash -l -c "cd ''%s'';', cygwin_home, work_dir);
    end
else
    if( ~exist( 'suffix', 'var' ) )
        work_dir = strrep(pwd(),'\','/');

        % It is necessary to set all ENV variables before invoking
        % pfstools commands. '/bin/bash' may need to be replaced with the
        % shell you are using.
        
        % This will remove all references to matlab libraries from the
        % LD_LIBRARY_PATH. pfstools usually do not work with matlab version
        % of the standard libraries 
        set_ld_path='export LD_LIBRARY_PATH=`echo $LD_LIBRARY_PATH | sed "y/:/\n/" | grep --ignore-case -v "matlab" | sed ":beg;N;s/\n/:/;t beg;"`';
        
        if( strcmp( computer, 'MACI64' ) || strcmp( computer, 'MACI32' ) )
            set_ld_path=cat( 2, set_ld_path, '; export DYLD_FRAMEWORK_PATH=`echo $DYLD_FRAMEWORK_PATH | sed "y/:/\n/" | grep --ignore-case -v "matlab" | sed ":beg;N;s/\n/:/;t beg;"`' );
        end
        
        % Start shell, cd to the working directory and remove matlab
        % libraries
        cmd = sprintf('/bin/bash -l -c ''cd "%s"; %s;', work_dir, set_ld_path );
    else
        cmd='''';
    end
end

end
