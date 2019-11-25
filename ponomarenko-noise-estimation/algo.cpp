/*
 *
 * This file is part of the Ponomarenko Noise Estimation algorithm.
 *
 * Copyright(c) 2011 Miguel Colom.
 * miguel.colom@cmla.ens-cachan.fr
 *
 * This file may be licensed under the terms of of the
 * GNU General Public License Version 2 (the ``GPL'').
 *
 * Software distributed under the License is distributed
 * on an ``AS IS'' basis, WITHOUT WARRANTY OF ANY KIND, either
 * express or implied. See the GPL for the specific language
 * governing rights and limitations.
 *
 * You should have received a copy of the GPL along with this
 * program. If not, go to http://www.gnu.org/licenses/gpl.html
 * or write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <fftw3.h>
#include <iostream>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "algo.h"
#include "curve_filter.h"
//
#include "framework/CFramework.h"
#include "framework/CImage.h"
#include "framework/libparser.h"
#include "framework/operations.cpp"
#include "framework/CHistogram.cpp"

//! Computes the delta matrix and returns the normalization factor theta
/*!
  \param *delta Delta matrix (mask for the low/high freqs in the block)
  \param w Block side
  \param T Number of low-freq coefficients, excluding DC
  \return theta Normalization factor for the matrix delta
*/
int compute_delta(float *delta, int w, int T) {
  int theta = 0;
  for (int j = 0; j < w; j++)
    for (int i = 0; i < w; i++) {    
      int value = (i + j < T && i + j != 0 ? 1 : 0);
      delta[j*w+i] = value;
      theta += value;
    }
  return theta;
}

//! Computes the set of variances computed form the low-frequency coefficients of the given blocks
/*!
  \param *VL Output set of variances
  \param M number of blocks taken into account
  \param w Block side
  \param *delta Delta matrix (mask for the low/high freqs in the block)s
  \param **blocks_ptr List of pointers to the blocks
  \param theta Normalization factor for the matrix delta
*/
void compute_VL(float *VL, int M, int w, float *delta, float **blocks_ptr,
                        int theta) {
  for (int m = 0; m < M; m++) {
    float *block = blocks_ptr[m];
    VL[m] = 0;

    for (int j = 0; j < w; j++) {
      for (int i = 0; i < w; i++)
        if (delta[j*w+i] != 0)
          VL[m] += pow(block[j*w+i], 2);
    }
    VL[m] /= theta;    
  }
}

//! Computes the set of variances computed form the high-frequency coefficients of the given blocks
/*!
  \param *VH Output set of variances
  \param **blocks_ptr List of pointers to the blocks
  \param *indices_VL Sorting indices for the blocks_ptr list (by low-freqs)
  \param w Block side
  \param T Number of low-freq coefficients, excluding DC
  \param K Number of blocks that should be used
  \return Length of the returned variances list
*/
int compute_VH(float *VH, float **blocks_ptr, int *indices_VL, int w,
                int T, int K) {
  int VH_count = 0;
  
  //#pragma omp parallel for      
  for (int q = 0; q < w*w; q++) {
    int j = q / w;
    int i = q - j*w;

    if (i + j >= T) {
      float s = 0.0;
      for (int k = 0; k < K; k++) {
        float *block = blocks_ptr[indices_VL[k]];
        s += pow(block[q], 2); // q == j*w+i
      }
      VH[VH_count++] = s / K;              
    }      
  }
  return VH_count;
}

//! Computes the optimal K parameter using Ponomarenko's original article loop
/*!
  \param M Number of variance values in VL to use
  \param *VL List of variances obtained for low-freq coefficients
  \return The optimal K
*/
int get_optimal_K_ponom_orig(int M, float *VL) {
  int K = sqrt(M);
  //
  for (int i = 0; i < 7; i++) {
    float U = 1.3 * VL[K/2];

    int m_min = arg_find<float>(U, VL, M);

    int K1 = m_min;
    if (K1 > 0)
      K = K1;
  }

  // Set K = K / 5 to provide robustness
  int K1 = int(K / 5.0);
  if (K1 > 0)
    K = K1;
  
  return K;
}

//! Return the optimal T parameter according to the given block side
/*!
  \param w Block side
  \return The optimal T parameter
*/
int get_T(int w) {
  switch (w) {
    case 3: return 3;
    case 5: return 5;
    case 7: return 8;
    case 8: return 9;
    case 11: return 13;
    case 15: return 17;
    case 21: return 24;
    default:
      PRINT_ERROR("Unknown block side: %d\n", w);
      exit(-1);
  }
}

//! Reads all valid blocks (all neighbor pixels are different when the mask
//! is active) in the image
/*!
  \param *D Output list of blocks
  \param *u Input image
  \param Nx Length of a row in the image
  \param Ny Length of a column in the image
  \param w Block side
  \param num_blocks Number of blocks
  \param *mask Mask to determine if a pixel is valid or not
  \return Number of valid block copied into the output list
*/
void read_all_valid_blocks(float *D,
                           float *u,
                           int Nx, int Ny,
                           int w, unsigned num_blocks, int *mask) {  
  if (mask == NULL) {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (unsigned q = 0; q < num_blocks; q++) {
      for (int j = 0; j < w; j++) {
        for (int i = 0; i < w; i++) {
            D[q*w*w+j*w+i] = u[j*Nx+i+q];
        }
      }
    }    
  }
  else {
    unsigned *valid_coords = new unsigned[num_blocks];
    int count_coords = 0;
    //    
    for (int i = 0; i < Nx*Ny; i++) {
      if (mask[i] == 0)
        valid_coords[count_coords++] = i;
    }
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (unsigned q = 0; q < num_blocks; q++) {
      int addr = valid_coords[q];      

      for (int j = 0; j < w; j++) {
        for (int i = 0; i < w; i++) {
            D[q*w*w+j*w+i] = u[j*Nx+i+addr];
        }
      }
    }
    delete[] valid_coords;
  }
}

//! Computes the mean of all given blocks
/*!
  \param *means Output list of means of blocks
  \param *blocks Input list of blocks to compute their means
  \param w Block side
  \param num_blocks Number of block in the input list
*/
void compute_means(float *means, float *blocks, int w, int num_blocks) {
  float ONE_DIV_w2 = 1.0 / (w*w);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int b = 0; b < num_blocks; b++) {
    float mean = 0.0;
    for (int p = 0; p < w*w; p++) {
      mean += blocks[b*w*w+p];
    }
    mean *= ONE_DIV_w2;
    means[b] = mean;
  }
}

int get_max(int *data, int N) {
  int max = data[0];
  for (int i = 1; i < N; i++)
    if (data[i] > max)
      max = data[i];
  return max;
}

void copy(float *dest, float *orig, int N) {
  for (int i = 0; i < N; i++)
    dest[i] = orig[i];
}

//! Determines if the given string corresponds to the custom percentile code
/*!
  \param *num_str Input string
  \return true if the input string corresponds to the custom percentile code or false if not.
*/
bool is_custom_percentile(const char *num_str) {
  char buffer[1024];
  float value = atof(num_str);
  sprintf(buffer, "%.4f", value);
  return strcmp(buffer, "0.0000") != 0;
}

//! Returns the mean associated to the given bin
/*!
  \param mean_method Method to compute the mean (1: mean of means, 2: median of means)
  \param K Number of elements in the bin
  \param indices_VL List of variances computed with low-freq coefficients
  \param histo Histogram object
  \return The mean of the bin
*/
float get_bin_mean(int mean_method, int K, int *indices_VL,
                   int bin, CHistogram<float*> *histo) {
  float *means_bin = histo->get_datal_bin(bin);
  float bin_mean;
  
  float *values = new float[K+1];
  for (int i = 0; i <= K; i++)
    values[i] = means_bin[indices_VL[i]];

  switch (mean_method) {
    case 1: { // mean of means
      bin_mean = 0.0;
      for (int i = 0; i <= K; i++) {
        bin_mean += values[i];
      }
      bin_mean /= (K+1);
      break;
    }
    case 2: { // median of means
      bin_mean = median<float>(values, K+1);
      break;
    }
    default: {
      PRINT_ERROR("Unknown mean method: %d\n", mean_method);
      exit(-1);
    }
  }
  delete[] values;
  return bin_mean;
}

//! In-place Normalization of the FFTW output in order to get a orthonormal 2D DCT-II
/*!
  \param *blocks Input/output list of transformed blocks
  \param w Block side
  \param num_blocks Number of blocks in the list
*/
void normalize_FFTW(float *blocks, int w, int num_blocks) {
  const float ONE_DIV_2w = 1.0 / (2.0*w);
  const float ONE_DIV_SQRT_2 = 1.0 / sqrtf(2);
  
  // Divide all coefficients by 2*w
  //#pragma omp parallel for shared(blocks)
  for (int i = 0; i < num_blocks*w*w; i++)  
    blocks[i] *= ONE_DIV_2w;
  
  #ifdef _OPENMP
  #pragma omp parallel for shared(blocks) schedule(static)
  #endif
  for (int b = 0; b < num_blocks; b++) {
    // {(i, j)} with i == 0 or j == 0
    for (int j = 0; j < w; j++) {
      int i = 0;
      blocks[b*w*w+j*w+i] *= ONE_DIV_SQRT_2;
    }
    for (int i = 0; i < w; i++) {
      int j = 0;
      blocks[b*w*w+j*w+i] *= ONE_DIV_SQRT_2;
    }
  }
}

/**
 * @brief Build a mask for valide pixel. If mask(i, j) = true, the pixels will not be used.
 *
 * @param i_im : noisy image;
 * @param o_mask : will contain the mask for all pixel in the image size;
 * @param p_imSize : size of the image;
 * @param p_sizePatch : size of a patch.
 *
 * @return number of valid blocks.
 *
 **/
unsigned buildMask(CImage &i_im, int *o_mask,
                    unsigned Nx, unsigned Ny, unsigned w,
                    unsigned num_channels) {
  unsigned count  = 0;

  for (unsigned ij = 0; ij < Nx*Ny; ij++) {
    const unsigned j = ij / Nx;
    const unsigned i = ij - j * Nx;

    //! Look if the pixel is not to close to the boundaries of the image
    if (i < Nx - w + 1 && j < Ny - w + 1) {
      for (unsigned c = 0; c < num_channels; c++) {
        float *u = i_im.get_channel(c);

        //! Look if the square 2x2 of pixels is constant
        int invalid_pixel = (c == 0 ? 1 : o_mask[ij]);

        // Try to validate pixel
        if (fabs(u[ij] - u[ij + 1]) > 0.001f) {
                invalid_pixel = 0;
        }
        else
        if (fabs(u[ij + 1] - u[ij + Nx]) > 0.001f) {
                invalid_pixel = 0;
        }
        else
        if (fabs(u[ij + Nx] - u[ij + Nx + 1]) > 0.001f) {
                invalid_pixel = 0;
        }
        o_mask[ij] = invalid_pixel;
      }
    }
    else {
            o_mask[ij] = 1; // Not valid
    }

    if (o_mask[ij] == 0)
        count++;
  }

  return count;
}

//! Ponomarenko et al. AVIRIS noise estimation algorithm.
/*!
  \param argc Number of arguments of the program
  \param **argv Arguments of the program
*/
void algorithm(int argc, char **argv) {
  vector <OptStruct *> options;
  vector <ParStruct *> parameters;
  //
  OptStruct owin   =  {"w:", 8, "8", NULL, "Block side"};
  options.push_back(&owin);
  OptStruct opercentile = {"p:", 1, "0.005", NULL, "Percentile"};
  options.push_back(&opercentile);
  OptStruct ore = {"r", 0, NULL, NULL, "Flag to remove equal pixels"};
  options.push_back(&ore);  
  OptStruct obins = {"b:", 0, "0", NULL, "Number of bins"};
  options.push_back(&obins);
  OptStruct oD = {"D:", 7, "7", NULL, "Filtering distance"};
  options.push_back(&oD);
  OptStruct ofiltercurve = {"g:", 5, "5", NULL, "Filter curve iterations"};
  options.push_back(&ofiltercurve);  
  OptStruct omeanMethod = {"m:", 2, "2", NULL, "Mean computation method"};
  options.push_back(&omeanMethod);  

  ParStruct pinput = {"input", NULL, "input file"};
  parameters.push_back(&pinput);
  //
  if (!parsecmdline("ponomarenko", "Ponomarenko SD noise estimation algorithm",
          argc, argv, options, parameters)) {
    printf("\n");
    printf("(c) 2012 Miguel Colom. Under license GNU GPL.\n");
    printf("http://mcolom.perso.math.cnrs.fr/\n");
    printf("\n");
    exit(-1);
  }

  // Read parameters
  int w = atoi(owin.value);
  int T = get_T(w);
  float p = atof(opercentile.value);
  int num_bins = atoi(obins.value);
  int D = atoi(oD.value);
  int curve_filter_iterations = atoi(ofiltercurve.value);
  int mean_method = atoi(omeanMethod.value);
  bool remove_equal_pixels_blocks = ore.flag;

  // Parallelization config
  #ifdef _OPENMP
  omp_set_num_threads(omp_get_num_procs());
  #endif

  // Load input image
  CImage input;
  input.load((char*)pinput.value);
  
  // Get image properties
  int Nx = input.get_width();
  int Ny = input.get_height();
  int num_channels = input.get_num_channels();
  
  // Set number of bins
  if (num_bins <= 0) num_bins = Nx * Ny / 42000;
  if (num_bins <= 0) num_bins = 1; // Force at least one bin  
  
  // Custom percentile or given by the user?
  bool custom_percentile = is_custom_percentile(opercentile.value);
  
  int total_blocks = (Nx-w+1) * (Ny-w+1); // Number of overlapping blocks

  // Create equal pixels mask
  int *mask_all;
  int num_blocks;

  if (remove_equal_pixels_blocks) {
    mask_all = new int[Nx*Ny];
    num_blocks = buildMask(input, mask_all, Nx, Ny, w, num_channels);
  }
  else {
    mask_all = NULL;
    num_blocks = total_blocks;
  }
  
  // Compute delta and theta
  CFramework *fw = CFramework::get_framework();
  float *delta = fw->create_array(w*w);
  int theta = compute_delta(delta, w, T);
  
  // Arrays to store the final means and noise estimations
  float *vmeans = new float[num_channels * num_bins];
  float *vstds  = new float[num_channels * num_bins];
  
  float *means = new float[num_blocks];  
  float *blocks = new float[num_blocks*w*w];
  
  // Init FFTW threads
  fftwf_init_threads();
  
  int nbTable[2] = {w,w};
  int nembed[2] = {w,w};
  
  #ifdef _OPENMP
  fftwf_plan_with_nthreads(omp_get_num_procs());
  #endif

  fftwf_r2r_kind kindTable[2] = {FFTW_REDFT10, FFTW_REDFT10};

  fftwf_plan fft_plan = fftwf_plan_many_r2r(2, nbTable, num_blocks, blocks,
                                            nembed, 1, w*w, blocks, nembed,
                                            1, w*w, kindTable, FFTW_ESTIMATE);

  // Process each channel
  for (int ch = 0; ch < num_channels; ch++) {
    float *u = input.get_channel(ch);
    
    read_all_valid_blocks(blocks, u, Nx, Ny, w, num_blocks, mask_all);

    // Compute means
    compute_means(means, blocks, w, num_blocks);

    // Compute 2D-DCT of all the blocks
    //
    // Transform blocks with FFTW
    fftwf_execute_r2r(fft_plan, blocks, blocks);    

    // Normalize FFTW output
    normalize_FFTW(blocks, w, num_blocks);

    // Create a list of pointers of the groups
    float **blocks_ptr = new float*[num_blocks];
    for (int i = 0; i < num_blocks; i++)
      blocks_ptr[i] = &blocks[i*w*w];

    // Create histogram according to the means
    CHistogram<float*> histo = CHistogram<float*>(num_bins,
                                                  blocks_ptr,
                                                  means,
                                                  num_blocks);

    // Process each bin
    #ifdef _OPENMP
    #pragma omp parallel for shared(vmeans, vstds, histo) schedule(static)
    #endif
    for (int bin = 0; bin < num_bins; bin++) {
      int elems_bin = histo.get_num_elements_bin(bin);

      float **block_ptr_bin = histo.get_data_bin(bin);

      float *VL = new float[elems_bin];

      // Compute VL
      compute_VL(VL, elems_bin, w, delta, block_ptr_bin, theta);

      // Get optimal K
      int K;
      if (custom_percentile)
        K = elems_bin * p;
      else // Using Ponomarenko's article loop
        K = get_optimal_K_ponom_orig(elems_bin, VL);

      // Compute VH
      int *indices_VL = new int[elems_bin];
      argsort(VL, indices_VL, elems_bin);

      // Array VH
      float *VH = new float[w*w];

      int VH_count = compute_VH(VH, block_ptr_bin, indices_VL, w, T, K);

      float bin_mean = get_bin_mean(mean_method, K, indices_VL,
                                    bin, &histo);
      float tilde_sigma = sqrt(median(VH, VH_count));

      // Store results
      vmeans[ch*num_bins+bin] = bin_mean;
      vstds[ch*num_bins+bin] = tilde_sigma;

      delete[] VL;
      delete[] VH;
      delete[] indices_VL;        
    }
    
    delete[] blocks_ptr;
  }
  
  // Filter noise curve
  float *new_std_control = new float[num_bins * num_channels];
  copy(new_std_control, vstds, num_channels*num_bins);
  //
  for (int ch = 0; ch < num_channels; ch++)
    for (int filt_iter = 0; filt_iter < curve_filter_iterations; filt_iter++) {
      bool allow_up = (filt_iter < 3);
      filter_curve(&vmeans[ch*num_bins], &new_std_control[ch*num_bins],
                   num_bins,
                   &new_std_control[ch*num_bins],
                   D, allow_up);
    }
  
  // Print results
  for (int bin = 0; bin < num_bins; bin++) {
    // Means
    for (int ch = 0; ch < num_channels; ch++)
        printf("%f  ", vmeans[ch*num_bins+bin]);
        
    // Standard deviations
    for (int ch = 0; ch < num_channels; ch++)
        printf("%f  ", new_std_control[ch*num_bins+bin]);
    //
    printf("\n");
  }

  // FFTW Cleanup
  fftwf_destroy_plan(fft_plan);
  fftwf_cleanup_threads();
  fftwf_cleanup();  

  // Clean up memory
  if (mask_all != NULL) delete[] mask_all;
  delete[] new_std_control;  
  delete[] vmeans;
  delete[] vstds;
  delete[] means;
  delete[] blocks;
 }
