%PFS_TRANSFORM_COLORSPACE Tranform between color spaces using pfs library.
%
% [C1 C2 C2] = PFS_TRANSFORM_COLORSPACE( inCSname, c1, c2, c3, outCSname );
% img_out = PFS_TRANSFORM_COLORSPACE( inCSname, img_in, outCSname );
% img_out = PFS_TRANSFORM_COLORSPACE( inCSname, c1, c2, c3, outCSname );
% [C1 C2 C2] = PFS_TRANSFORM_COLORSPACE( inCSname, img_in, outCSname );
%
%   inCSname - name of the input color space
%   c<n> - matrix with n-th channel of input color space
%   C<n> - matrix with n-th channel of output color space
%   img_in - input image given as 3D height/width/3 matrix 
%   img_out - output image given as 3D height/width/3 matrix 
%   outCSname - name of the output color space
%
% Recognized color space names:
% 'XYZ' - CIE 1931 XYZ
% 'RGB' - linear RGB with rec.709 primaries
% 'sRGB' - sRGB colorspace, the maximum value is 1
% 'YUV' - Y is luminance, uv and u'v' chromacity coordinates of the CIE L'u'v' colorspace
% 'Yxy' - CIE 1931 XYZ, but as chromatic coordinates
% 'RGB2020' - linear RGB with rec.2020 primaries
% 'PQYCbCr2020' - HDR10 MPEG color coding with the PQ transfer function. It assumes
%         10,000 cd/m^2 peak luminance. The resulting values are not quantized and
%         are not rescaled to the valid range. 
% 'YCbCr709' - colorspace used for SDR MPEG coding.
% Color space names are case insensitive.
%
% See also: PFS_READ_IMAGE, PFS_WRITE_IMAGE, PFSVIEW.
%
% Copyright 2009 Rafal Mantiuk
