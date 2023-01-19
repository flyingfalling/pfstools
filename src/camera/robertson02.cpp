/**
 * @brief Robertson02 algorithm for automatic self-calibration.
 *
 * 
 * This file is a part of PFS CALIBRATION package.
 * ---------------------------------------------------------------------- 
 * Copyright (C) 2004 Grzegorz Krawczyk
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
 * @author Grzegorz Krawczyk, <gkrawczyk@users.sourceforge.net>
 *         Ivo Ihrke        , <ihrke@mmci.uni-saarland.de>
 *
 * $Id: robertson02.cpp,v 1.11 2011/02/25 13:45:14 ihrke Exp $
 */


#include <config.h>

#include <iostream>
#include <vector>
#include <cstdlib>

#include <math.h>

#include <responses.h>
#include <robertson02.h>

#include <iostream>
#include <cstring>

using namespace std;


#define PROG_NAME "robertson02"

// maximum iterations after algorithm accepts local minima
#define MAXIT 100

// maximum accepted error
#define MAX_DELTA 1e-5f

extern bool verbose; /* verbose should be declared for each standalone code */

float normalize_rcurve( float *rcurve, int M );

inline float max3( float a[3])
{
  float max = (a[0]>a[1]) ? a[0] : a[1];
  return (a[2]>max) ? a[2] : max;
}

float get_exposure_compensation( const Exposure &ex )
{
  return ex.exposure_time * ex.aperture*ex.aperture * ex.iso/3.125f / (1.0592f * 11.4f / 3.125f);
}

int robertson02_applyResponseRGB( pfs::Array2D *rgb_out[],         
  const ExposureList *imgs[], 
  const float *resp_curve[], 
  const float *weights, 
  int M,
  const noise_parameters *camera,
  const float deghosting_threshold)
{
  // number of exposures
  int N = imgs[0]->size();
  
  // frame size
  int width  = rgb_out[0] -> getCols( );
  int height = rgb_out[0] -> getRows( );

  // number of saturated pixels
  int saturated_pixels = 0;

  float mmax[3]= {1e-30, 1e-30, 1e-30};

/*
  float pw[M];
  for( int i = 0 ; i < M ; i++ ) {
    pw[i] = weights[i];
  }
  for( int i = 1; i < 10; i++ )
    pw[M-i] = 0;
  
  float ew[N]; // Per-exposure weight to account for noise charactetistic

  // For each exposure
  for( int i = 0 ; i < N ; i++ ) {
    const Exposure &ex = (*imgs[0])[i];

    // TODO: Estimate read-out noise and use for exposure weighting
    // TODO: Better blending if difference is weights is high
	
    // This is a noise-optimial reconstruction assuming that the input is affected only
    // by the photon noise (no read-out noise), which is most often the case
    ew[i] = ex.exposure_time; 
      
    VERBOSE_STR << "Ex: " << i << " exp. time: " << ex.exposure_time << "\t iso: " << ex.iso 
                << "\t aperture: " <<  ex.aperture << "\t weight: " << ew[i] << std::endl;
  }
*/

  bool photon_noise_only = false;
  if( strcmp( camera->camera_name, "Empty" ) == 0 )
  {
    photon_noise_only = true;
    VERBOSE_STR << "Merging using the Poisson Photon Noise Estimator" << endl;
  }
  else
    VERBOSE_STR << "Merging using the Variance Weighted Estimator" << endl;

  // For each pixel
  int skipped_deghost = 0;
  for( int j = 0; j < width * height; j++ )
  {
    bool saturatedRGB[3] = { false, false, false };

    // First compute scene radiances, weights and saturated exposures
    float X[3][N];
    float w[3][N];
    bool saturated_exp[N];
    fill(saturated_exp, saturated_exp + N, false);
    bool skip_exp[N];
    fill(skip_exp, skip_exp + N, false);
    float t[N];
    float g[N];

    for( int i = 0; i < N; i++ )
      for( int cc = 0; cc < 3; cc++ )
      {
        const Exposure &ex = (*imgs[cc])[i];
        t[i] = ex.exposure_time;
        g[i] = ex.iso / 100;
      
        int   m  = (int) (*ex.yi)( j );
        X[cc][i] = resp_curve[cc][ m ] / ( t[i] * g[i] * camera->color_coefficients[cc] );
        if( photon_noise_only )
          w[cc][i] = t[i];
        else
        {
          float var = resp_curve[cc][ m ] * g[i] * camera->color_coefficients[cc]
                      + powf(camera->std_readout * g[i] * camera->color_coefficients[cc], 2)
                      + powf(camera->std_adc * camera->color_coefficients[cc], 2);
            w[cc][i] = powf( t[i] * g[i] * camera->color_coefficients[cc], 2 ) / var;
        }

        // Test for overexposed pixels
        saturated_exp[i] |= resp_curve[cc][m] >= resp_curve[cc][M-1];
      }

    // Deghosting - find the exposure that could be used as a reference
    // Currently using middle exposure
    // TODO: Perform histogram analysis to pick best reference
    int ref = int(N/2);
    if( deghosting_threshold != -1 )
    {
      // Compute the number of pixels to be skipped due to deghosting
      for( int i = 0; i < N; i++ )
        if(i != ref && !saturated_exp[ref] && !saturated_exp[i])
        {
          // First compute the Pearson correlation coefficient wrt reference
          float ref_mean = 0, ref_var = 1e-5, cur_mean = 0, cur_var = 1e-5, cov = 0;
          for( int cc = 0; cc < 3; cc++ )
          {
            ref_mean += X[cc][ref];
            cur_mean += X[cc][i];
          }
          ref_mean /= 3;
          cur_mean /= 3;
          for( int cc = 0; cc < 3; cc++ )
          {
            ref_var += powf(X[cc][ref] - ref_mean, 2);
            cur_var += powf(X[cc][i] - cur_mean, 2);
            cov += (X[cc][ref] - ref_mean)*(X[cc][i] - cur_mean);
          }

          // Next check whether total deviation from the reference is acceptable
          float std_sum = 0;
          float deviation_from_ref = 0;
          for( int cc = 0; cc < 3; cc++ )
          {
            // Reciprocal of the weight is the exposure compensated variance
            // Square root of this is thus the exposure compensated standard deviation
            // Assuming a gaussian distribution, 95% of the population should lie 2 stds away from mean
            std_sum += deghosting_threshold * sqrt(1/w[cc][ref] + 1/w[cc][i]);
            deviation_from_ref += fabs( (float)(X[cc][ref] - X[cc][i]) );
          }

          // Skipping the Pearson coefficient test as reults are not satisfactory
          // Fixing a better reference might solve this
          // float rho = cov / sqrt(ref_var/2 * cur_var/2);
          // skip_exp[i] = rho < 0.9 || deviation_from_ref > std_sum;
          skip_exp[i] = deviation_from_ref > std_sum;
          if( skip_exp[i] )
            skipped_deghost++;
        }
    }
    	
/*     // Deghosting
      if( ref != -1 && !saturated_exp[ref] ) {

        const Exposure &ex = (*imgs[0])[ref];          
        float ti = get_exposure_compensation( ex );
      
        float Y_ref = (X[ref*3] + X[ref*3+1] + X[ref*3+2])/3.f;
        const float noise_level = 0.01;      
        float var_ref = (noise_level*noise_level)*Y_ref/ti;      // variance of the photon noise
      
        for( int i = 0 ; i < N ; i++ )
        {

          const Exposure &ex = (*imgs[0])[i];          
          float ti = get_exposure_compensation( ex );
        
          if( saturated_exp[i] ) {
//          skip_exp[i] = true;          
            continue;
          }
        
          //continue; // Skip deghostig for now
		
          float rho = corr_vec3( &X[i*3], &X[ref*3] );
//     float chroma = std::max( std::max( X[i*3], X[i*3+1] ), X[i*3+2] ) - std::min( std::min( X[i*3], X[i*3+1] ), X[i*3+2] );     
//     cerr << rho << " - " << chroma << ", ";
//        cerr << i << "(" << saturated_exp[i] << "): " << rho << endl;

          float Y_exp = (X[i*3] + X[i*3+1] + X[i*3+2])/3.f;
          float var_exp = (noise_level*noise_level)*Y_exp / ti;

          float ci = 1.96 * sqrt(var_ref+var_exp);
        
          // Skip a frame if color correlation too low or difference to reference too high
          skip_exp[i] = ( rho < 0.9 ) || fabs( Y_exp - Y_ref ) > ci;        
          if( skip_exp[i] ) 
            skipped_deghost++;
        }
      }
    
    
      int skip_bits = 0;
    
      for( int i = 0 ; i < N ; i++ ) {

      if( !skip_exp[i] )
      skip_bits++; // |= (1<<i);
      }
    

    
     for( int i = 0 ; i < N ; i++ ) {        
     skip_exp[i] = false;
     }
*/
      
    
//   cerr << endl;    
    
    // For each color channels
    for( int cc=0; cc < 3; cc++ )
    {

      // all exposures for each pixel
      float sum = 0.0f;
      float div = 0.0f;

      // For each exposure
      for( int i = 0 ; i < N ; i++ )
      {
        if( saturated_exp[i] || (deghosting_threshold != -1 && skip_exp[i]) )
          continue;

        sum += X[cc][i] * w[cc][i];
        div += w[cc][i];

//        if( j < 10 && cc == 1 )
//         std::cerr << "j = " << j << " ;i = " << i << " ;w = " << w[cc][i] << std::endl;
      }

      if( div >= 1e-15 )
      {
        (*rgb_out[cc])(j) = sum / div;
        mmax[cc] = ( mmax[cc] > (*rgb_out[cc])(j) ) ? mmax[cc] : (*rgb_out[cc])(j);
      }
      else
      {
        (*rgb_out[cc])(j) = -1;
        saturatedRGB[cc] = true;
      }
    }
    

    if( saturatedRGB[0] || saturatedRGB[1] || saturatedRGB[2] ) {      
      saturated_pixels++;      
    }
    
/*  for( int cc=0; cc < 3; cc++ ) {
      // TODO: Some fancy stuff to deal with distorted colors
      
      if( saturatedRGB[cc] ) {        
        // If none of the exposures can give actual pixel value, make
        // the best (clipped) estimate using the longest or the
        // shortest exposure;

        const Exposure &ex = (*imgs[cc])[0];
        // If pixel value > gray level, use the shortest exposure;        
        float short_long = ( (*ex.yi)(j) > M/2 ) ? 1.f : -1.f;        
        
        float best_ti = 1e10f * short_long;
        int best_v = short_long == 1.f ? M-1 : 0;        
        for( int i = 0 ; i < N ; i++ )
	{
          if( deghosting && skip_exp[i] )
            continue;
          
          const Exposure &ex = (*imgs[cc])[i];          
	  int   m  = (int) (*ex.yi)( j );
	  float ti = get_exposure_compensation( ex );;          
          if( ti*short_long < best_ti*short_long ) {            
            best_ti = ti;
            best_v = (int)(*ex.yi)(j);
          }          
	}        
        (*rgb_out[cc])(j) = 1/best_ti * resp_curve[cc][best_v];
      }
      
    }*/
    
    
  }

  // Fill in nan values and normalise
  for( int j = 0; j < width * height; j++ )
    for( int cc=0; cc < 3; cc++ )
    {
      if( (*rgb_out[cc])(j) == -1)
        (*rgb_out[cc])(j) = mmax[cc];
      (*rgb_out[cc])(j) /= max3(mmax);
    }

  VERBOSE_STR << "Exposure pixels skipped due to deghosting: " << 
    (float)skipped_deghost*100.f/(float)(width*height*N) << "%" << std::endl;

  return saturated_pixels;
  
}



int robertson02_applyResponse( pfs::Array2D *xj,         
  const ExposureList &imgs, 
  const float *rcurve, 
  const float *weights, 
  int M )
{

  // number of exposures
  int N = imgs.size();

  // frame size
  int width  = xj -> getCols( );
  int height = xj -> getRows( );

  // number of saturated pixels
  int saturated_pixels = 0;

//  cerr << "M = " << M << endl;
//  cerr << "W[0] = " << weights[0] << endl;
//  cerr << "W[end] = " << weights[M-1] << endl;
  
  // all pixels
  for( int j = 0; j < width * height; j++ )
  {
    // all exposures for each pixel
    float sum = 0.0f;
    float div = 0.0f;

    for( int i = 0 ; i < N ; i++ )
    {
      int   m  = (int) (*imgs[ i ].yi)( j );
      float ti = imgs[ i ].ti;
	   
      sum += weights [ m ] / ti * rcurve [ m ];
      div += weights [ m ];
    }

    
/*      if( div != 0.0f )
	(*xj)( j ) = sum / div;
        else
        (*xj)( j ) = 0.0f;*/

    if( div >= 0.01f )
      (*xj)( j ) = sum / div;
    else
      (*xj)( j ) = 0;

    /*
      {
      float best_ti = 1e10;
      int best_v = M-1;        
      for( int i = 0 ; i < N ; i++ )
      {
      int   m  = (int) (*imgs[ i ].yi)( j );
      float ti = imgs[ i ].ti;          
      if( ti < best_ti ) {            
      best_ti = ti;
      best_v = (int)(*imgs[i].yi)(j);
      }          
      }        
      (*xj)( j ) = (M-1) / best_ti * rcurve [best_v];        
        
      }*/
      
      
    
  }

  return saturated_pixels;

}


int robertson02_getResponse( pfs::Array2D       *xj, 
  const ExposureList &imgs,
  float              *rcurve, 
  const float        *weights, 
  int M )
{

// number of exposures
  int N = imgs.size();
    
  // frame size
  int width  = imgs[0].yi -> getCols( );
  int height = imgs[0].yi -> getRows( );

  // number of saturated pixels
  int saturated_pixels = 0;

  // indices
  int i, j, m;
  
  float* rcurve_prev = new float[ M ];	// previous response

  if( rcurve_prev == NULL )
  {
    std::cerr << "robertson02: could not allocate memory for camera response" << std::endl;
    exit(1);
  }

  // 0. Initialization
  
  normalize_rcurve( rcurve, M );
  
  for( m = 0 ; m < M ; m++ ) {
    //  cerr << "m = " << m << " rc = " << rcurve [ m ] << endl;    
    rcurve_prev [ m ] = rcurve [ m ];
  }
  
  robertson02_applyResponse( xj, imgs, rcurve, weights, M );

  // Optimization process
  bool   converged  = false;
  long  *cardEm     = new long [ M ];
  float *sum        = new float[ M ];

  if( sum == NULL || 
    cardEm == NULL )
  {
    std::cerr << "robertson02: could not allocate memory for optimization process" << std::endl;
    exit(1);
  }
  
  int   cur_it  = 0;
  float pdelta  = 0.0f;

  while( !converged )
  {

    // Display response curve - for debugging purposes
/*    for( m = 0 ; m < M ; m+=32 ) {
      cerr << "m = " << m << " rc = " << rcurve [ m ] << endl;
      }*/
    
    // 1. Minimize with respect to rcurve
    for( m = 0 ; m < M ; m++ )
    {
      cardEm [ m ] = 0;
      sum[ m ] = 0.0f;
    }
    
    // For each exposure
    for( i = 0 ; i < N ; i++ )
    {

      pfs::Array2D* yi = imgs[i].yi;
      float         ti = imgs[ i ].ti;

      for( j = 0 ; j < width * height ; j++ )
      {
        m = (int) (*yi)( j );

        if( m < M && m >= 0 )
        {
          sum[ m ] += ti * (*xj)( j );
          cardEm[ m ] ++;
        }
        else
          std::cerr << "robertson02: m out of range: " << m << std::endl;
      }
    }


    for( m = 0; m < M ; m++ )
    {
      if( cardEm[ m ] != 0 )
        rcurve [ m ] = sum [ m ] / cardEm [ m ];
      else
        rcurve [ m ] = 0.0f;
    }

    // 2. Normalize rcurve
    float middle_response = normalize_rcurve( rcurve, M );    
    
    // 3. Apply new response
    saturated_pixels = robertson02_applyResponse( xj, imgs, rcurve, weights, M );
    
    // 4. Check stopping condition
    float delta = 0.0f;
    int   hits  = 0;

    for( m = 0 ; m < M ; m++ )
    {
      if( rcurve[ m ] != 0.0f )
      {
        float diff = rcurve [ m ] - rcurve_prev [ m ];
	      
        delta += diff * diff;
	      
        rcurve_prev [ m ] = rcurve[ m ];
	      
        hits++;
      }
    }

    delta /= hits;

    VERBOSE_STR << " #" << cur_it
                << " delta=" << delta
                << " (coverage: " << 100*hits/M << "%)\n";
    
    if( delta < MAX_DELTA )
      converged = true;
    else if( isnan(delta) || cur_it > MAXIT )
    {
      VERBOSE_STR << "algorithm failed to converge, too noisy data in range\n";
      break;
    }

    pdelta = delta;
    cur_it++;
  }

  if( converged )
    VERBOSE_STR << " #" << cur_it
                << " delta=" << pdelta << " <- converged\n";

  delete[] rcurve_prev;
  delete[] cardEm;
  delete[] sum;
  
  return saturated_pixels;
}


//----------------------------------------------------------
// private part


int comp_floats( const void *a, const void *b )
{
  return ( (*((float*)a))< (*((float*)b)) ) ;
}

float normalize_rcurve( float* rcurve, int M )
{
  int   FILTER_SIZE =  M / 256;
  float mean;
  float rcurve_filt [ M ];
  float to_sort [ 2 * FILTER_SIZE + 1 ];

  mean = 0.f;
  for ( int i = 0; i < M; i++ )
  {
    mean += rcurve [ i ];
  }
  mean /= M;

  if( mean != 0.0f )
    for( int m = 0 ; m < M ; m++ )
    {
      rcurve [ m ] /= mean;

      /* filtered curve - initialization */
      rcurve_filt [ m ] = 0.0f;
    }

  /* median filter response curve */
  for ( int m = FILTER_SIZE ; m < M - FILTER_SIZE; m++ )
  {
    for ( int i = -FILTER_SIZE; i <= FILTER_SIZE; i++ )
      to_sort [ i + FILTER_SIZE ] = rcurve[ m + i ];
      
    qsort ( to_sort, 2 * FILTER_SIZE + 1 , sizeof(float), comp_floats );
      
    rcurve_filt [ m ]  = to_sort [ FILTER_SIZE ]; 
      
  }

  /* boundaries */
  for( int m = 0 ; m < FILTER_SIZE ; m++ )
  {
    rcurve_filt [ m ]     = rcurve_filt [ FILTER_SIZE ];
    rcurve_filt [ M - m - 1 ] = rcurve_filt [ M - FILTER_SIZE - 1 ];
  }
    
  /* copy curve */
  for( int m = 0 ; m < M ; m++ )
  {
    rcurve [ m ] = rcurve_filt [ m ];
  }

  return mean;
}


/*float normalize_rcurve_old( float* rcurve, int M )
  {

  int Mmin, Mmax;
  // find min max
  for( Mmin=0 ; Mmin<M && rcurve[Mmin]==0 ; Mmin++ );
  for( Mmax=M-1 ; Mmax>0 && rcurve[Mmax]==0 ; Mmax-- );
  
  int Mmid = Mmin+(Mmax-Mmin)/2;
  float mid = rcurve[Mmid];

//   std::cerr << "robertson02: middle response, mid=" << mid
//             << " [" << Mmid << "]"
//             << " " << Mmin << ".." << Mmax << std::endl;
  
if( mid==0.0f )
{
// find first non-zero middle response
while( Mmid<Mmax && rcurve[Mmid]==0.0f )
Mmid++;
mid = rcurve[Mmid];
}

if( mid!=0.0f )
for( int m=0 ; m<M ; m++ )
rcurve[m] /= mid;
return mid;
}
*/
