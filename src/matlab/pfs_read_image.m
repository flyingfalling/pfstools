function img = pfs_read_image( file_name, pfsin_command, extra_arguments )
%PFS_READ_IMAGE Read image file and return RGB, luminance or multichannel matrix.
%
% IMG = PFS_READ_IMAGE( file_name )
% IMG = PFS_READ_IMAGE( file_name, pfsin_command )
% IMG = PFS_READ_IMAGE( file_name, pfsin_command, extra_arguments )
%
% Read an HDR or SDR image using pfstools. 
%
% pfsin_command could be one of pfsin* command tools, such as 'pfsindcraw' or
% 'pfsinimgmagick'. By default 'pfsin' is used, which can automatically
% detect the image format from an extension and call the right pfsin
% commnand. This option can be useful if pfsin does not recognize the image
% format or you want to pass extra options to the command with
% "extra_arguments".
%
% extra_arguments is a string that specified extra options to be passed to
% pfsin command. For example, when used in combination with 'pfsindcraw',
% you can pass '-n' to return the image in the native camera color space
% instead of rec.709 RGB. 
%
% If input is a gray-scale or luminance image, IMG is a 2D matrix. If input is
% a color image, IMG(:,:,k) represents red, blue and green color channels for k=1,2 and 3.
% If input is a multi-channel image (channel names C1, C2, ..., Cn), IMG is a
% 3D matrix with 3rd dimension corresponding to the channels. 
%
% PFS_READ_IMAGE accepts all formats recognized by the shell "pfsin"
% command.
%
% Example: 
%   img = PFS_READ_IMAGE( 'hdr_images/memorial.exr' );
%
% See also: PFS_READ_RGB, PFS_READ_LUMINANCE, PFS_READ_XYZ,
% PFS_WRITE_IMAGE, PFSVIEW.
%
% Copyright 2009 Rafal Mantiuk


if ~exist( 'pfs_command', 'var' )
    pfsin_command = 'pfsin';
end

if ~exist( 'extra_arguments', 'var' )
    extra_arguments = '';
end

  %Check if file exists
  fid = fopen( file_name, 'rb' );
  if( fid == -1 ) 
    error( 'pfs_read_image: File "%s" does not exist', file_name );
  end
  fclose( fid );

  fid = pfspopen( sprintf( '%s%s %s ''%s''%s', pfs_shell(), pfsin_command, extra_arguments, pfs_fix_path(file_name), pfs_shell( 1 ) ), 'r' );
  pin = pfsopen( fid );
  pin = pfsget( pin );

  if( isfield( pin.channels, 'X' ) && isfield( pin.channels, 'Z' ) )
      img = pfs_transform_colorspace( 'XYZ', pin.channels.X, pin.channels.Y, pin.channels.Z, 'RGB' );
  elseif( isfield( pin.channels, 'Y' ) )
      img = pin.channels.Y;
  elseif( isfield( pin.channels, 'C1' ) )
      ch=1;
      % count the number of channels
      while( isfield( pin.channels, sprintf( 'C%d', ch ) ) )
          ch = ch+1;
      end
      ch_max = ch-1;
      img = zeros(pin.rows, pin.columns, ch_max);
      for ch=1:ch_max
          img(:,:,ch) = pin.channels.(sprintf( 'C%d', ch ));
      end
  else
      error( 'Color channels missing in the pfs frame' );
  end  
  
  pfsclose( pin );
  % TODO: Check why crashes on windows
  if ~ispc()
      pfspclose( fid );
  end
end
