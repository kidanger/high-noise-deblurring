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

#ifndef CURVE_FILTER_H
#define	CURVE_FILTER_H

//! Returns the standard deviation at intensity mu given the
//! ending points (mu_c, std_c),
//! (mu_e, std_e) using an affine transformation.
/*!
  \param mu_c Current point intensity
  \param std_c Current point standard deviation
  \param mu_e Ending point intensity
  \param std_e Ending point standard deviation
  \param mu Intensity of the given point
  \return Standard deviation of the point at intensity mu
*/
float get_affine_std(float mu_c, float std_c,
                       float mu_e, float std_e,
                       float mu);

//! Filters the given noise curve
/*!
  \param mus_c Intensities of the control points
  \param stds_c Standard deviations of the control points
  \param num_bins Number of control points
  \param new_std_control Output filtered standard deviations
  \param D Interval diameter
  \param allow_up Controls if the points can go up and down or only down
*/
void filter_curve(float *mus_c, float *stds_c, int num_bins,
                   float *new_std_control,
                   int D, bool allow_up);

#endif	/* CURVE_FILTER_H */

