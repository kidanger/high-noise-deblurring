/*
 *
 * This file is part of the Image Processing Framework.
 *
 * Copyright(c) 2011 Miguel Colom.
 * Website: http://mcolom.info ; miguel at mcolom.info
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

#ifndef _CHISTOGRAM_H_
#define	_CHISTOGRAM_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "CHistogram.cpp"
#include "operations.cpp"

//! Histogram class
template <class T>
class CHistogram {
public:
    //! Creates a histogram
    /*!
      \param bins Number of bins
      \param *data Data that goes inside the bins
      \param *datal Data to used to compute the histogram limits
      \param N Number of elements in data and datal
    */
    CHistogram<T>(int bins, T *data, float *datal,
               int N);
    virtual ~CHistogram<T>();
    //
    //! Gets bin begining
    /*!
      \param bin Bin number
      \return bin beginning
    */
    float get_limit_begin(int bin);

    //! Gets bin end
    /*!
      \param bin Bin number
      \return bin end
    */
    float get_limit_end(int bin);

    //! Gets the number of samples in the bin
    /*!
      \param bin Bin number
      \return number of samples in the bin
    */
    int get_num_elements_bin(int bin);

    //! Gets the elements (data) inside the bin
    /*!
      \param bin Bin number
      \return pointer to the elements (data) in the bin
    */

    T *get_data_bin(int bin);

    //! Gets the elements (datal) inside the bin
    /*!
      \param bin Bin number
      \return pointer to the elements (datal) in the bin
    */
    float *get_datal_bin(int bin);

    //! Gets the number of bins in the histogram
    /*!
      \return number of bins in the histogram
    */
    int get_num_bins();
private:
    int bins;
    float *limits_begin;
    float *limits_end;
    float *num_elements;
    T **data_bins;
    float **datal_bins;
};

template <class T>
CHistogram<T>::CHistogram(int bins, T *data, float *datal, int N) {
  this->bins = bins;

  int samples_per_bin = (int)floor(N / bins);
  if (N % bins != 0)
    samples_per_bin++;

  // Assign memory
  this->limits_begin = new float[bins];
  this->limits_end   = new float[bins];

  this->num_elements = new float[bins];

  this->data_bins = new T*[bins];
  this->datal_bins = new float*[bins];

  // Initialize this->data_bins and this->datal_bins to NULL
  for (int i = 0; i < this->bins; i++) {
    this->data_bins[i] = NULL;
    this->datal_bins[i] = NULL;
  }
  
  // Allocate samples for each bin.
  T *buffer = new T[N];
  float *bufferl = new float[N];
  
  // Sort data by datal
  int *indices = new int[N];
  argsort(datal, indices, N);

  // Number of elements in the current bin  
  int elements_count = 0;

  int bin = 0;
  for (int idx = 0; idx < N; idx++) {      
      // Keep loading...
      buffer[elements_count] = data[indices[idx]];
      bufferl[elements_count] = datal[indices[idx]];
      
      elements_count++;
      
      if (idx == N-1 || elements_count == samples_per_bin) {
          // Buffer for the bin full

          // Copy data to current bin
          this->datal_bins[bin] = new float[elements_count];

          this->data_bins[bin] = new T[elements_count];
          memcpy(this->data_bins[bin],
                  buffer,
                  sizeof(T) * elements_count);

          memcpy(this->datal_bins[bin],
                  bufferl,
                  sizeof(float) * elements_count);

          // Write bin metadata
          this->limits_begin[bin] = bufferl[0];
          this->limits_end[bin] = bufferl[elements_count-1];
          this->num_elements[bin] = elements_count;

          /* printf("[%.3f, %.3f], %d\n",
               this->limits_begin[bin], this->limits_end[bin],
               elements_count); */

          // Prepare for next bin
          elements_count = 0;
          bin++;
      }     
  }
  
  // Free memory
  delete[] indices;
  delete[] buffer;
  delete[] bufferl;
}

template <class T>
CHistogram<T>::~CHistogram() {
  // Delete data bins
  for (int i = 0; i < this->bins; i++) {
      if (this->data_bins[i] != NULL) delete[] this->data_bins[i];
      if (this->datal_bins[i] != NULL) delete[] this->datal_bins[i];
  }
  delete[] this->data_bins;
  delete[] this->datal_bins;
  
  // Delete other arrays
  delete[] this->limits_begin;
  delete[] this->limits_end;
  
  delete[] this->num_elements;
}

template <class T>
float CHistogram<T>::get_limit_begin(int bin) {
    return this->limits_begin[bin];
}

template <class T>
float CHistogram<T>::get_limit_end(int bin) {
    return this->limits_end[bin];
}

template <class T>
int CHistogram<T>::get_num_elements_bin(int bin) {
    return this->num_elements[bin];
}

template <class T>
T *CHistogram<T>::get_data_bin(int bin) {
    return this->data_bins[bin];
}

template <class T>
float *CHistogram<T>::get_datal_bin(int bin) {
    return this->datal_bins[bin];
}

template <class T>
int CHistogram<T>::get_num_bins() {
    return this->bins;
}

#endif	/* CHISTOGRAM_H */
