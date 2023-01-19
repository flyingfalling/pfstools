 /**
 * @brief Write files in Radiance's RGBE format
 * 
 * This file is a part of PFSTOOLS package.
 * ---------------------------------------------------------------------- 
 * Copyright (C) 2003,2004 Rafal Mantiuk and Grzegorz Krawczyk
 * 
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * ---------------------------------------------------------------------- 
 * 
 * @author Grzegorz Krawczyk, <krawczyk@mpi-sb.mpg.de>
 * @author Rafal Mantiuk, <mantiuk@mpi-sb.mpg.de>
 *
 * $Id: pfsoutrgbe.cpp,v 1.3 2014/06/17 21:57:09 rafm Exp $
 */

#include <config.h>
#include <string.h>

#include <cstdlib>

#include <iostream>
#include <cmath>
#include <getopt.h>
#include <pfs.h>
#include <pfsutils.cpp>

#define PROG_NAME "pfsoutyuv"

template<typename T>
T clip3(T x, T y, T z){
  if(z<x){
    return x;
  }
  if(z>y){
    return y;
  }
  return z;
}

class MetaDataFrameFileIterator : pfs::FrameFileIteratorImpl
{
  
  unsigned int width;
  unsigned int height;
  unsigned int bitdepth;
  unsigned int fps;
  const char * colorSpace;
  const char * chromaFormat;
  //Used for YUV output - we'll override this from FrameFileIterator - because we only want one YUV file writing to - and the rest of the class is useful
  
  public:
  MetaDataFrameFileIterator( int &argc, char* argv[], const char *fopenMode,
    const char *fileNamePrefix, FILE *stdinout,
    const char *optstring, const struct option *getopt_long, unsigned int width, unsigned int height, unsigned int bitdepth, const char * colorSpace, const char * chromaFormat, int fps) : 
    pfs::FrameFileIteratorImpl(argc, argv, fopenMode, fileNamePrefix, stdinout, optstring, getopt_long), bitdepth(bitdepth), colorSpace(colorSpace), width(width), height(height), chromaFormat(chromaFormat), fps(fps) 
  {
     
  }

  pfs::FrameFile getNextFrameFile(){
    if(!currentPattern->isPattern){
      strcpy(fileName, currentPattern->pattern);
    }
    else{ 
      char metaData[100]; //reasonably safe?
      
      if(fps>0){
        sprintf(metaData, "%dx%d_%dfps_%db_%s_%s", width, height, fps, bitdepth, colorSpace, chromaFormat);
      }
      else{ 
        sprintf(metaData, "%dx%d_%db_%s_%s", width, height, bitdepth, colorSpace, chromaFormat);
      }

      sprintf(fileName,currentPattern->pattern, metaData);
    
    }
    //grab prefix 
    FILE *fh;
    if( currentPattern->stdinout != NULL )
      fh = currentPattern->stdinout;
    else
      fh = fopen( fileName, fopenMode );
        if( fh != NULL ){
         return pfs::FrameFile( fh, fileName );
        }
        else{
          throw pfs::Exception("Cound't open file for saving");
        }
        return pfs::FrameFile(NULL, NULL);
    }
  
  void closeFrameFile(pfs::FrameFile ff){
    pfs::FrameFileIteratorImpl::closeFrameFile(ff);
  }
};

class FilterableArray2D : public pfs::Array2D {
  Array2D * base;
  public:
    FilterableArray2D(pfs::Array2D * in){
      base = in;
    } 

    int getCols() const{
      return base->getCols();
    }
    int getRows() const{
      return base->getRows();
    }
    
    float& operator()( int col, int row ){
      if(col < 0 || col >= base->getCols() || row < 0 || row >= base->getRows()){
        return *(new float(0));
      }
      else{
        return (*base)(col,row);
      }
    }

    const float& operator()( int col, int row ) const{
      if(col < 0 || col >= base->getCols() || row < 0 || row >= base->getRows()){
        return *(new float(0));
      }
      else{
        return (*base)(col,row);
      }
    }
    
    float& operator()( int index ){
      return (*base)(index);
    }
    
    const float& operator()( int index ) const{
      return (*base)(index); 
    }
};

class QuietException 
{
};

enum DownsampleFilter{f0, f1};


void printHelp()
{
  std::cerr << PROG_NAME " [--verbose] [--quiet] [--bitdepth] [--colorspace] [--downsample-filter] [--chroma-format]" 
            << std::endl
            << "See man page for more information." << std::endl;
}

void downscale_444_to_420(pfs::Array2D * channel, DownsampleFilter method){
  float weightsF0[] = {0.125, 0.75, 0.125}; //weights for f0 filter method
  float weightsF1[] = {0.25, 0.5, 0.25}; //ditto for f1
  float * weights;
  switch(method){
    case f0 : weights = weightsF0; 
    case f1 : weights = weightsF1;
  }
  //downscales a single channel
  //do rows first
  
  FilterableArray2D chan(channel);
  for (int i = 0; i<chan.getCols(); i++){
    for(int j = 0; j<chan.getRows(); j+=2){
      chan(i,j) = weights[0] * chan(i, j-1) + weights[1] * chan(i,j) + weights[2]* chan(i, j+1); 
    }
  }
  //do columns second
  for( int j = 0; j<chan.getRows(); j+=2){
    for( int i = 0; i<chan.getCols(); i+=2){
      chan(i,j) = weights[0] * chan(i-1, j) + weights[1]* chan(i,j) + weights[2]* chan(i+1, j);
    }
  }
  //the downsampled pixels are now stored in the even cells of channel
}

class YUVWriter {
  
  FILE * fh;
  int bitdepth;
  bool downscale_chroma;
  DownsampleFilter downsampleFilter;
  public:
  YUVWriter(FILE * fh, int bitdepth, bool downscale_chroma, DownsampleFilter downsampleFilter) : fh(fh), bitdepth(bitdepth), downscale_chroma(downscale_chroma), downsampleFilter(downsampleFilter){  
  }
  

  template<typename T>
  void writeYuvChannel(pfs::Array2D * in, int width, int height, int stride, float scale, float offset){
  for( int j = 0; j<height; j++){ 
    T * line_buf = new T[width];

    for (int i = 0; i<width; i++){
      line_buf[i] = clip3(0, (1<<bitdepth) - 1,
                  (int) round((*in)(i<<stride % width, j<<stride % height) * scale + offset));
    }
    fwrite(line_buf, sizeof(T), width, fh);
    delete[] line_buf;
  }
  //performs floating point to 10 bit conversion
  }

  template<typename T>
  void writeYuvImage(pfs::Array2D *R, pfs::Array2D *G, pfs::Array2D *B){
  //writes the entire image to fh R - Luminance, G - Cr, B - Cb
  unsigned int width_chroma = R->getCols();
  unsigned int height_chroma = R->getRows();
  unsigned int step_chroma = 0; 
  if(downscale_chroma){
    width_chroma = R->getCols()/2;
    height_chroma = R->getRows()/2;
    step_chroma = 1;
    downscale_444_to_420(G, downsampleFilter);
    downscale_444_to_420(B, downsampleFilter);
  }
  writeYuvChannel<T>(R, R->getCols(), R->getRows(), 0, 219.0f*(1<<(bitdepth-8)), 1<<(bitdepth -4));
  writeYuvChannel<T>(G, width_chroma, height_chroma, step_chroma, 224.0f*(1<<(bitdepth-8)), 1<<(bitdepth -1));
  writeYuvChannel<T>(B, width_chroma, height_chroma, step_chroma, 224.0f*(1<<(bitdepth-8)), 1<<(bitdepth -1));
  }
};

void writeFrames( int argc, char* argv[] )
{
  pfs::DOMIO pfsio;

  bool verbose = false;
  int bit_depth = 10;
  unsigned int width; 
  unsigned int height;
  unsigned int fps = 0;
  bool quiet = false;
  bool downscale_chroma = true;
  pfs::ColorSpace opt_colorspace = pfs::CS_INVALID;
  DownsampleFilter downsampleFilter = f0;
  bool srgb_input = false;  // Input is in sRGB colour space if true, RGB otherwise
  bool opt_srgb_input = false;  
  bool opt_linear_input = false;
  
  
  // Parse command line parameters
  static struct option cmdLineOptions[] = {
    { "help", no_argument, NULL, 'h' },
    { "verbose", no_argument, NULL, 'v' },
    { "bitdepth", required_argument, NULL, 'b' },
    { "colourspace", required_argument, NULL, 'c' },
    { "quiet", no_argument, NULL, 'q' },
    { "downsample-filter", required_argument, NULL, 'a'},
    { "chroma-format", required_argument, NULL, 'f'},
    { "srgb-input", no_argument, NULL, 's'},
    { "linear-input", no_argument, NULL, 'l'},
    { NULL, 0, NULL, 0 } };
  static const char optstring[] = "hb:vqc:f:a:sl";
    
  int optionIndex = 0;
  while( 1 ) {
    int c = getopt_long (argc, argv, optstring, cmdLineOptions, &optionIndex);
    if( c == -1 ) break;
    switch( c ) {
    case 'h':
      printHelp();
      throw QuietException();
    case 'v':
      verbose = true;
      break;
    case 'b':
      bit_depth = strtol( optarg, NULL, 10 );
      break;
    case 'a':
      if( strcmp( optarg, "f0") == 0){
        downsampleFilter = f0;}
      else if( strcmp( optarg, "f1") == 0){
              downsampleFilter = f1;}
      break;
    case 'f':
       if( strcmp( optarg, "444") == 0){
         downscale_chroma = false;
       }
       else if ( strcasecmp( optarg, "420") == 0){
         downscale_chroma = true;
         
       }
       else{
        throw pfs::Exception( "Unrecognised chroma-sampling format" );
       }
       break;
  case 'c':
    if( !strcasecmp(optarg, "pq2020" )) {
      opt_colorspace = pfs::CS_PQYCbCr2020;
    } else if( !strcasecmp(optarg, "bt709" )) {
      opt_colorspace = pfs::CS_YCbCr709;
    } else if( !strcasecmp(optarg, "hlg2020" )) {
      opt_colorspace = pfs::CS_HLGYCbCr2020;
    } else{
      throw pfs::Exception( "Unrecognized colorspace name" );
    }
    break;
    case 'q':
      quiet = true;
      break;
    case 's':
      opt_srgb_input = true;
      break;
    case 'l':
      opt_linear_input = true;
      break;
      
  case '?':
      throw QuietException();
    case ':':
      throw QuietException();
    }
  }


  //we need to open this here to grab metadata for the filename :)
  pfs::Frame *frame = pfsio.readFrame( stdin );  
  pfs::Channel *fullDimensionChannel;
  fullDimensionChannel = frame->getChannel("X");    
  if( fullDimensionChannel == NULL ) {       // No color
    throw pfs::Exception( "Missing X channel in the PFS stream" );
  }
  width = fullDimensionChannel->getCols();
  height = fullDimensionChannel->getRows();
  
  const char* luminanceTag = frame->getTags()->getString("LUMINANCE");
  const char* fpsTag = frame->getTags()->getString("FPS");
  if( opt_colorspace == pfs::CS_INVALID ) {
    // Try to guess the right colour space
    if( luminanceTag!=NULL && 0==strcmp(luminanceTag,"DISPLAY") ) {      
      opt_colorspace = pfs::CS_YCbCr709;
    } else if( luminanceTag!=NULL && 0==strcmp(luminanceTag,"ABSOLUTE") )  {
        opt_colorspace = pfs::CS_PQYCbCr2020;
    } else
      throw pfs::Exception( "Color space must be specified" );
  }

  const char * colorSpaceStr;
  switch( opt_colorspace ) {
  case pfs::CS_YCbCr709:
    colorSpaceStr = "bt709";
    break;    
  case pfs::CS_PQYCbCr2020:
    colorSpaceStr = "pq2020";
    break;    
  case pfs::CS_HLGYCbCr2020:
    colorSpaceStr = "hlg2020";
    break;    
  default:
    throw pfs::Exception( "Unrecognized color space" );    
  }

  if(fpsTag != NULL){
    VERBOSE_STR << "FPS Tag: " << fpsTag << std::endl;
    fps = atoi(fpsTag);
  }

  if( opt_srgb_input ) {    
    srgb_input = true;
  } else if( opt_linear_input ) {
    srgb_input = false;    
  } else {
    if( opt_colorspace == pfs::CS_YCbCr709 ) 
      srgb_input = true;    
  }
  
    if( !srgb_input && luminanceTag!=NULL && 0==strcmp(luminanceTag,"DISPLAY") ) {      
      std::cerr << "pfsoutyuv: Warning - assuming linear (HDR) input while the stream may be gamma-corrected (LDR)." << std::endl;
    }
    if( srgb_input && luminanceTag!=NULL && 0==strcmp(luminanceTag,"ABSOLUTE") ) {      
      std::cerr << "pfsoutyuv: Warning - assuming gamma-corrected (LDR) input while the stream may be linear (HDR)." << std::endl;
    }
    if( luminanceTag!=NULL && 0==strcmp(luminanceTag,"RELATIVE") ) {      
      std::cerr << "pfsoutyuv: Warning - relative HDR pixel values detected. Absolute values should be provided." << std::endl;
    }
    
  
  if( srgb_input )
    VERBOSE_STR << "Input data in sRGB (LDR) space" << std::endl;
  else
    VERBOSE_STR << "Input data in absolute linear (HDR) space" << std::endl;
  
  MetaDataFrameFileIterator it(argc, argv, "wb", NULL, stdout, optstring, cmdLineOptions, width,
                                height, bit_depth, colorSpaceStr, downscale_chroma ? "420" : "444", fps);
  
  pfs::FrameFile ff = it.getNextFrameFile();
  if( ff.fh == NULL ) {
    throw pfs::Exception("Couldn't open file for writing");
  } else {
      VERBOSE_STR << "Saving file: " << ff.fileName << std::endl << width << "x" << height << ", Bit-depth: " << bit_depth << ", Colorspace: " << getColorspaceString(opt_colorspace)  << ", " << (downscale_chroma ? "420, " : "444, ") << std::endl;
  }

  while(true){
    YUVWriter writer( ff.fh, bit_depth, downscale_chroma, downsampleFilter);

    pfs::Channel *X, *Y, *Z;
    frame->getXYZChannels( X, Y, Z );    
    if( X == NULL || Y == NULL || Z == NULL) {       // No color
      throw pfs::Exception( "Missing X, Y, Z channels in the PFS stream" );
    }
    if( srgb_input ) {
      pfs::transformColorSpace( pfs::CS_XYZ, X, Y, Z, pfs::CS_RGB, X, Y, Z );
      pfs::transformColorSpace( pfs::CS_SRGB, X, Y, Z, opt_colorspace, X, Y, Z );      
    } else{
      pfs::transformColorSpace( pfs::CS_XYZ, X, Y, Z, opt_colorspace, X, Y, Z );
    }
    if( bit_depth <= 8){
      writer.writeYuvImage<unsigned char>( X, Y, Z );
    }
    else{
      writer.writeYuvImage<unsigned short>( X, Y, Z );
    }

    pfsio.freeFrame( frame );

    frame = pfsio.readFrame( stdin );
    if( frame == NULL ) {
      break; // No more frames
    }
  }
  it.closeFrameFile( ff );
}


int main( int argc, char* argv[] )
{
  try {
    writeFrames( argc, argv );
  }
  catch( pfs::Exception ex ) {
    std::cerr << PROG_NAME << " error: " << ex.getMessage() << std::endl;
    return EXIT_FAILURE;
  }        
  catch( QuietException  ex ) {
    return EXIT_FAILURE;
  }        
  return EXIT_SUCCESS;
}
