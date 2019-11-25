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

#ifndef _OPERATIONS_H_
#define _OPERATIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

////////////////////////////////////////////////////////////////////////
//                         Header                                     //
////////////////////////////////////////////////////////////////////////


//! Fills the array with the specied value
/*!
  \param *dest Pointer to the array
  \param value Value to fill the array with
*/
template <typename T>
void set_array(T *dest, int value, int size);

//! Copies the contents of buffer orig into buffer dest
/*!
  \param *dest Origin buffer
  \param *dest Destiny buffer
  \param N Number of elements to copy
*/
template <typename T>
void copy_buffer(T *dest, T *orig, int N);

//! Sorts the first array and switches the second array elements at the same time using the Quicksort algorithm
/*!
  \param *arr Pointer to the first array
  \param *brr Pointer to the second array
  \param ilength Length of the arrays
*/
template <typename T1, typename T2>
void quicksort(T1 *arr, T2 *brr, int ilength);

//! Returns the list of indices that would sort the array
/*!
  \param *data Pointer to the array
  \param *indices Pointer to the indices array
  \param N Length of the array
*/
template <typename T>
void argsort(T *data, int *indices, int N);

//! Computes the Summed Area Table of the array
/*!
  \param *src Pointer to the source matrix array
  \param *dst Pointer to the destiny array to save the SAT to
  \param Nx Width of the input matrix
  \param Ny Height of the input matrix
*/
template <typename TS, typename TD>
void summed_area_table(TS *src, TD *dst, int Nx, int Ny);

//! Returns the index for which the store element in the array has the lowest distance with the given value
/*!
  \param value Value to look for
  \param *array The array to look into
  \param N Number of elements in the array
*/
template <typename T>
int arg_find(T value, T *array, int N);

//! Computes the mean of the elements of the array
/*!
  \param *array Input array
  \param N Number of elements in the array
  \return The mean of the N elements of the array
*/
template <typename T>
T mean(T *array, int N);

//! Computes the median of the elements of the array
/*!
  \param *array Input array
  \param N Number of elements in the array
  \return The median of the N elements of the array
*/
template <typename T>
T median(T *array, int N);

//! Computes the subtraction of two arrays
/*!
  \param *s Output array with the subtracted elements
  \param *a Input array "a" to compute "s=a-b"
  \param *b Input array "b" to compute "s=a-b"
  \param N Number of elements in the arrays
*/
template <typename T>
void subtract(T *s, T *a, T *b, int N);

//! Computes the Median of Absolute Diferences (MAD operator)
/*!
  \param *data Input array
  \param N Number of elements in the input array
  \return The MAD of the input data
*/
template <typename T>
T compute_MAD(T *data, int N);



////////////////////////////////////////////////////////////////////////
//                         Implementation                             //
////////////////////////////////////////////////////////////////////////

template <typename T>
void set_array(T *dest, int value, int size) {
  for (int i = 0; i < size; i++)
    dest[i] = value;
}

template <typename T>
void copy_buffer(T *dest, T *orig, int N) {
  for (int i = 0; i < N; i++)
    dest[i] = orig[i];
}

template <typename T>
inline void swap(T *x, T *y) {
  T aux = *x;
  *x=*y;
  *y=aux;
}

/* Quicksort,  values in arr are set in increasing order and brr elements are
 * switched at the same time*/
template <typename T1, typename T2>
void quicksort(T1 *arr, T2 *brr, int ilength) {
  int M=7,NSTACK=50;
  int i,ir,j,k,jstack=-1,l=0;
  T1 a, b;
  int istack[50];
  
  ir=ilength-1;
  
  
  for(;;){
    if(ir-l<M){
      for(j=l+1;j<=ir;j++){
        a=arr[j];
        b=brr[j];
        for(i=j-1;i>=l;i--){
          if (arr[i]<=a) break;
          arr[i+1]=arr[i];
          brr[i+1]=brr[i];
        }
        arr[i+1]=a;
        brr[i+1]=b;
        
      }
      
      if (jstack<0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      
      k=(l+ir) >> 1;
      swap(&arr[k],&arr[l+1]);
      swap(&brr[k],&brr[l+1]);
      if (arr[l]>arr[ir]){
        swap(&arr[l],&arr[ir]);
        swap(&brr[l],&brr[ir]);
      }
      if (arr[l+1]>arr[ir]){
        swap(&arr[l+1],&arr[ir]);
        swap(&brr[l+1],&brr[ir]);
      }
      if (arr[l]>arr[l+1]){
        swap(&arr[l],&arr[l+1]);
        swap(&brr[l],&brr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for(;;){
        do i++; while (arr[i]<a);
        do j--; while (arr[j]>a);
        if (j<i) break;
        swap(&arr[i],&arr[j]);
        swap(&brr[i],&brr[j]);
      }
      
      arr[l+1]=arr[j];
      arr[j]=a;
      brr[l+1]=brr[j];
      brr[j]=b;
      jstack+=2;
      
      if (jstack>=NSTACK) {
        fprintf(stderr, "quicksort: stack too small!\n");
        exit(-1);
      }

      if (ir-i+1>=j-l){
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
      
    }
  }  
}

template <typename TS, typename TD>
TD direct_acc_sum(TS *src, int x, int y, int Nx) {
  TD sum = 0;
  for (int j = 0; j <= y; j++) {
    for (int i = 0; i <= x; i++) {
      sum += src[j*Nx + i];
    }
  }
  return sum;
}

template <typename TS, typename TD>
void summed_area_table(TS *src, TD *dst, int Nx, int Ny) {
  // y = 0
  double acc = 0;
  for (int x = 0; x < Nx; x++) {
    acc += src[x];
    dst[x] = acc;
  }

  // x = 0
  acc = 0;
  for (int y = 0; y < Ny; y++) {
    acc += src[y * Nx];
    dst[y * Nx] = acc;
  }
  
  int accumulation = 0;
  
  for (int y = 1; y < Ny; y++) {
    for (int x = 1; x < Nx; x++) {
      if (accumulation == 500000) {
        // Sync. to avoid accumulating errors
        dst[y*Nx + x] = direct_acc_sum<TS, TD>(src, x, y, Nx);
        dst[y*Nx + x-1] = direct_acc_sum<TS, TD>(src, x-1, y, Nx);
        dst[(y-1)*Nx + x] = direct_acc_sum<TS, TD>(src, x, y-1, Nx);
        dst[(y-1)*Nx + x-1] = direct_acc_sum<TS, TD>(src, x-1, y-1, Nx);

        accumulation = 0;
      }
      else {
        dst[y*Nx + x] = src[y * Nx + x] +
                          dst[y * Nx + x-1] +
                          dst[(y-1) * Nx + x] -
                          dst[(y-1) * Nx + x-1];
        accumulation++;        
      }
    }
  }
}

template <typename T>
void argsort(T *data, int *indices, int N) {
 T *data2 = new T[N];
 for (int i = 0; i < N; i++) {
   data2[i] = data[i];
   indices[i] = i;
 }
 quicksort<T>(data2, indices, N);
 delete[] data2;
}

template <typename T>
int arg_find(T value, T *array, int N) {
    float dist_min = 9E9;
    int idx_min = 0;
    for (int i = 0; i < N; i++) {
        float dist = fabs(array[i] - value);
        if (dist < dist_min) {
            dist_min = dist;
            idx_min = i;
        }
    }
    return idx_min;
}

template <typename T>
T mean(T *array, int N) {
    float sum = 0.0;
    for (int i = 0; i < N; i++)
        sum += array[i];
    return sum / N;
}

template <typename T>
T median(T *array, int N) {
    if (N <= 0) {
        fprintf(stderr, "Error: not enough samples to compute median\n");
        exit(-1);
    }
    if (N == 1) return array[0];

    int *indices = new int[N];
    argsort(array, indices, N);

    float res;
    if (N % 2 == 0)
        res =  (array[indices[N/2 - 1]] + array[indices[N/2]]) * 0.5;
    else
        res = array[indices[(N-1)/2]];

    delete[] indices;
    return res;
}

template <typename T>
void subtract(T *s, T *a, T *b, int N) {
  for (int i = 0; i < N; i++)
    s[i] = a[i] - b[i];
}

template <typename T>
T compute_MAD(T *data, int N) {
  T *buffer = new T[N];
  
  float med = median<T>(data, N);

  copy_buffer<T>(buffer, data, N);
  
  for (int i = 0; i < N; i++)
    buffer[i] = fabs(buffer[i] - med);
  
  T mad = median<T>(buffer, N);  

  delete[] buffer;
  
  return mad;
}

#endif
