/*
 *
 * This file is part of the Image Processing Framework.
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

#ifndef _CIMAGE_H_
#define _CIMAGE_H_

#include <vector>
#include <string>
//
#include "CFramework.h"

#define CIMAGE_FLOAT_RGB 0
#define CIMAGE_PNG 1
#define CIMAGE_LUM 2

using namespace std;

//! Image I/O class
class CImage {
  private:
    int Nx, Ny;
    int bits;
    int id;
    int *channels_id;
    int num_channels;
    
    static int unique_id;
    
    int get_file_type(char *filename);
    void load_float_RGB(char *filename);
    void load_LUM(char *filename);
    void load_PNG(char *filename);
    void load_iio(char *filename);
    void save_PNG(char *filename, int bits_per_channel);
    void save_float_RGB(char *filename);
  public:
    CImage();
    CImage(int Nx, int Ny,
           int bits, int num_channels);
    //
    ~CImage();
    
    //! Creates the image channels
    /*!
      \param num_channels Number of channels
      \param Nx Image width
      \param Ny Image height
    */
    void create_channels(int num_channels, int Nx, int Ny);

    //! Loads a image from a file
    /*!
      \param *filename File name of the image to load
    */
    void load(char *filename);

    //! Saves a image to a file
    /*!
      \param filename File name of the image to save.
      \param bits_per_channel Bits per channel (ex: 8, 16).
    */
    void save(char *filename, int bits_per_channel);

    //! Saves a image to a file
    /*!
      \param *filename File name of the image to save
    */
    void save(char *filename);

    //! Returns the number of channels in the image
    /*!
      \return Number of channels in the image
    */
    int get_num_channels();

    //! Returns the width of the loaded image
    /*!
      \return The width of the loaded image
    */
    int get_width();

    //! Returns the height of the loaded image
    /*!
      \return The height of the loaded image
    */
    int get_height();

    //! Returns the number of bits per channel in the loaded image
    /*!
      \return The number of bits per channel in the loaded image
    */
    int get_bits_per_channel();

    //! Returns the data array of the specified image channel
    /*!
      \param channel Channel index (0, 1, ...)
      \return Pointer to the channel data
    */
    float *get_channel(int channel);
};

#endif
