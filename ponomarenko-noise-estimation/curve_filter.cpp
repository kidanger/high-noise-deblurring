/*
 *
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

#include "curve_filter.h"
#include "framework/operations.cpp"
//
#include <math.h>
#include <vector>

using namespace std;

// Get the std. corresponding to mu
float get_affine_std(float mu_c, float std_c,
                       float mu_e, float std_e,
                       float mu) {
  float slope;
  if (abs(mu_c - mu_e) < 1E-6) // Avoid dividing by zero
    slope = 0;
  else {
    slope = (std_c - std_e) / (mu_c - mu_e);
  }

  float std = (mu - mu_e) * slope + std_e;
  return std;
}

// Given the intensities and stds of the noise curve, extrapolates the
// std. of the point whose intensity is mu
float get_interpolation(float *mus, float *stds, float mu, int N) {
  // First of all, find the nearest control point
  int i = arg_find(mu, mus, N);
  float m = mus[i];
  
  float m1, m2;
  float s1, s2;
  
  if (mu < m) {
    // we're on the right of mu
    if (i == 0) {
      i = 1;
      m = mus[i];
    }
    //
    m1 = mus[i-1];
    m2 = m;
    s1 = stds[i-1];
    s2 = stds[i];
  }
  else {
    // we're on the left of mu
    if (i >= N-1) {
      i = N-2;
      m = mus[i];
    }
    //
    m1 = m;
    m2 = mus[i+1];
    s1 = stds[i];
    s2 = stds[i+1];    
  }
  
  return get_affine_std(m1, s1, m2, s2, mu);
}

// Filters the noise curve
void filter_curve(float *mus_c, float *stds_c, int num_bins,
                   float *new_std_control,
                   int D, bool allow_up) {
  int num_control_points = 0;

  for (int b = 0; b < num_bins; b++) {
    float mu_current = mus_c[b];
    float std_current = stds_c[b];
    
    float left = mu_current - D;
    float right = mu_current + D;
    
    if (left < mus_c[0]) {
      float dist = mus_c[b] - mus_c[0];
      left = mu_current - dist;
      right = mu_current + dist;      
    }
    else {
      if (right > mus_c[num_bins-1]) {
        float dist = mus_c[num_bins-1] - mus_c[b];
        left = mu_current - dist;
        right = mu_current + dist;      
      }
    }
    
    // Add the extrapolated control points inside the
    // interval [left, right]
    float sum_window = 0;
    int L = 0;
    for (float mu = left; mu <= right; mu += 0.05) {
      sum_window += get_interpolation(mus_c, stds_c, mu, num_bins);
      L++;
    }
    
    float std_new = sum_window / (float)L; // mean
    
    float std_filtered;
    //
    if (allow_up)
      std_filtered = std_new;
    else
      std_filtered = (std_new < std_current ? std_new : std_current);

    new_std_control[num_control_points] = std_filtered;
    num_control_points += 1;
  }
}
