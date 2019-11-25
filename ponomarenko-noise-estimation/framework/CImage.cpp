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


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string.h>
#include <png.h>
extern "C" {
#include "iio.h"
}
//
#include "CImage.h"
#include "CFramework.h"

int CImage::unique_id = 0;

void CImage::create_channels(int num_channels, int Nx, int Ny) {
  CFramework *fw = CFramework::get_framework();

  this->channels_id = new int[num_channels];

  int N[2] = {Nx, Ny};  
  for (int ch = 0; ch < num_channels; ch++) {
    int id_matrix;
    fw->create_matrix(2, N, id_matrix);
    this->channels_id[ch] = id_matrix;
    PRINT_VERBOSE("I_%d creates matrix M_%d\n", this->id, id_matrix);
  }
  this->num_channels = num_channels;
}

CImage::CImage() {
  this->unique_id++;
  this->id = this->unique_id;
  this->channels_id = NULL;
  PRINT_VERBOSE("CImage created, I_%d\n", this->id);
}

CImage::CImage(int Nx, int Ny,
               int bits, int num_channels) {
  this->unique_id++;
  this->id = this->unique_id;
  //
  this->Nx = Nx;
  this->Ny = Ny;
  this->bits = bits;
  //
  PRINT_VERBOSE("CImage I_%d, %dx%d, %d bits, %d channels created\n",
    this->id, this->Nx, this->Ny, this->bits, num_channels);

  // Create data channels
  this->create_channels(num_channels, Nx, Ny);
}


CImage::~CImage() {
  PRINT_VERBOSE("CImage I_%d destroy\n", this->id);
  CFramework *fw = CFramework::get_framework();
  for (int ch = 0; ch < this->get_num_channels(); ch++) {
    int id_channel = this->channels_id[ch];
    fw->delete_matrix(id_channel);
  }
  if (this->channels_id != NULL)
    delete[] this->channels_id;
}

// Obtain format by extension of filename
int CImage::get_file_type(char *filename) {
  const char *exts[] = {".RGB", ".PNG", ".LUM"};
  const int num_exts = 3;

  int type = -1;
  for (int i = 0; i < num_exts; i++) {
    char *ext = (char*)exts[i];
    int indexCmp = strlen(filename) - strlen(ext);
    if (indexCmp >= 0) {
      char *extFile = (char*)&filename[indexCmp];
      if (strcasecmp(ext, extFile) == 0) {
        type = i;
  break;
      }
    }
  }
  return type;
}

void CImage::load(char *filename) {
  // Obtain format by extension of file
  int type = get_file_type(filename);

  //////////////////////// read format
  switch (type) {
    case CIMAGE_FLOAT_RGB:
      this->load_float_RGB(filename);
      break;
    case CIMAGE_PNG:
      this->load_PNG(filename);
      break;
    case CIMAGE_LUM:
      this->load_LUM(filename);
      break;
    default:
      this->load_iio(filename);
//    PRINT_ERROR("CImage: unknown file type: %s\n", filename);
//    exit(-1);
  }
  PRINT_VERBOSE("Loaded image \"%s\", %dx%d, %d channels, %d bits in I_%d\n",
    filename, this->Nx, this->Ny,
    this->get_num_channels(), this->get_bits_per_channel(),
    this->id);
}

void CImage::save(char *filename, int bits_per_channel) {
  // Obtain format by extension of file
  int type = get_file_type(filename);

  //////////////////////// read format
  switch (type) {
    case CIMAGE_FLOAT_RGB:
      this->save_float_RGB(filename);
      break;
    case CIMAGE_PNG:
      this->save_PNG(filename, bits_per_channel);
      break;
    case CIMAGE_LUM:
      PRINT_ERROR("CImage: can't save LUM images\n");
      exit(-1);
      break;
    default:
      PRINT_ERROR("CImage: unknown file type: %s\n", filename);
      exit(-1);
  }
  PRINT_VERBOSE("Saved image \"%s\", %dx%d, %d channels, %d bits\n",
    filename, this->Nx, this->Ny,
    this->get_num_channels(), this->get_bits_per_channel());
}

void CImage::save(char *filename) {
  this->save(filename, 8);
}

// Load float RGB image
void CImage::load_float_RGB(char *filename) {
  int num_channels;
  size_t dummy;
  
  CFramework *fw = CFramework::get_framework();

  FILE *infile = fopen(filename, "rb");
  if (!infile) {
    PRINT_ERROR("CImage: error opening float RGB file %s\n!", filename);
    exit(-1);
  }

  dummy = fread(&this->Nx, sizeof(this->Nx), 1, infile);
  dummy += fread(&this->Ny, sizeof(this->Ny), 1, infile);
  dummy += fread(&num_channels, sizeof(num_channels), 1, infile);
  dummy = dummy; // To prevent the warning

  this->bits = 16;
  PRINT_VERBOSE("CImage load RGB (%dx%d), %d channels, in I_%d\n",
    this->Nx, this->Ny, num_channels, this->id);

  int id_data_all;
  float *data_ALL = fw->create_array(num_channels*this->Nx*this->Ny, id_data_all);
  dummy = fread(data_ALL, sizeof(float), num_channels*this->Nx*this->Ny, infile);
  fclose(infile);

  // Create data channels
  create_channels(num_channels, Nx, Ny);

  int idx_data = 0;
  for (int ch = 0; ch < num_channels; ch++) {
    float *m = this->get_channel(ch);
    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
          m[y*this->Nx+x] = data_ALL[idx_data];
          idx_data++;
        }
    }
  }

  fw->delete_matrix(id_data_all);
}

void CImage::load_LUM(char *filename) {
  unsigned char lumWidthArray[4];
  unsigned char lumHeightArray[4];
  size_t dummy;

  FILE *infile = fopen(filename, "rb");
  if (!infile) {
    PRINT_ERROR("CImage: error opening LUM file %s\n!", filename);
    exit(-1);
  }

  int num_channels = 1;
  this->bits = 16;

  PRINT_VERBOSE("CImage load RGB (%dx%d), %d channels, in I_%d\n",
  this->Nx, this->Ny, num_channels, this->id);

  dummy = fread(&lumWidthArray, sizeof(char), 4, infile);
  this->Nx = lumWidthArray[0] + 256*lumWidthArray[1] +
              65536*lumWidthArray[2] + 16777216*lumWidthArray[3];
  dummy += fread(&lumHeightArray, sizeof(char), 4, infile);
  dummy = dummy; // To prevent the warning
  this->Ny = lumHeightArray[0] + 256*lumHeightArray[1] +
           65536*lumHeightArray[2] + 16777216*lumHeightArray[3];

  // Identify file format
  char lumFormat[5];
  dummy = fread(lumFormat, 1, 5, infile);
  if (strcmp(lumFormat, "12LI") != 0 &&
      strcmp(lumFormat, "16LI") != 0 &&
      strcmp(lumFormat, "16LU") != 0) {
    PRINT_ERROR("Unknown LUM format (not 12LI, 16LI, 16LU): %s\n", lumFormat);
    exit(-1);
  }
  
  // Skip offset
  fseek(infile, -2 * this->Nx * this->Ny, SEEK_END);
  
  // Create data channels
  create_channels(num_channels, Nx, Ny);

  // Read data and put it in the channels
  for (int ch = 0; ch < num_channels; ch++) {    
    float *m = this->get_channel(ch);

    unsigned char luminanceArray[2];

    int ptr = 0;
    for (int y = 0; y < this->Ny; y++)
      for (int x = 0; x < this->Nx; x++) {
        dummy = fread(&luminanceArray, 1, 2, infile);
        m[ptr++] = luminanceArray[0] + 256*luminanceArray[1];
      }
  }

  fclose(infile);
}

void CImage::load_iio(char *filename) {
  int num_channels;
  size_t dummy;

  int w, h, c;
  float *image = iio_read_image_float_split(filename, &w, &h, &c);
//w = 10; h = 10; c = 1;
//float *image = (float *)malloc(100 * sizeof(float));

  this->Nx = w;
  this->Ny = h;
  num_channels = c;

  this->bits = 16; // ??
  PRINT_VERBOSE("CImage load image (%dx%d), %d channels, in I_%d\n",
    this->Nx, this->Ny, num_channels, this->id);

  // Create data channels
  create_channels(num_channels, Nx, Ny);

  int idx_data = 0;
  for (int ch = 0; ch < num_channels; ch++) {
    float *m = this->get_channel(ch);
    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
          m[y*this->Nx+x] = image[idx_data];
          idx_data++;
        }
    }
  }

  free(image);
}

static void *read_png_abort(FILE * fp,
                            png_structp * png_ptr_p, png_infop * info_ptr_p)
{
    png_destroy_read_struct(png_ptr_p, info_ptr_p, NULL);
    if (NULL != fp && stdin != fp)
        (void) fclose(fp);
    return NULL;
}

// Load PNG image
// http://zarb.org/~gc/html/libpng.html
void CImage::load_PNG(char *filename) {
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep * row_pointers;
  const unsigned PNG_SIG_LEN = 4;

  png_byte header[8];    // 8 is the maximum size that can be checked
  size_t dummy;
  
  /* open file and test for it being a png */
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    PRINT_ERROR("error opening PNG file %s!\n", filename);
    exit(-1);
  }

  dummy = fread(header, 1, PNG_SIG_LEN, fp);
  dummy = dummy; // To prevent the warning
  if (png_sig_cmp(header, 0, PNG_SIG_LEN)) {
    PRINT_ERROR("file %s is not recognized as a PNG file!\n", filename);
    exit(-1);
  }

  /* initialize stuff */
  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr) {
    PRINT_ERROR("png_create_read_struct failed\n");
    exit(-1);
  }

  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    PRINT_ERROR("png_create_info_struct failed\n");
    exit(-1);
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    PRINT_ERROR("Error during setjmp\n");
    exit(-1);
  }

  png_init_io(png_ptr, fp);
  png_set_sig_bytes(png_ptr, PNG_SIG_LEN);

  png_read_info(png_ptr, info_ptr);
  
  //png_get_IHDR(png_ptr, info_ptr, &this->Nx, &this->Ny, &this->bits,
  //    color_type, NULL, NULL, NULL);

  this->Nx = png_get_image_width(png_ptr, info_ptr);
  this->Ny = png_get_image_height(png_ptr, info_ptr);
  //color_type = png_get_color_type(png_ptr, info_ptr);
  this->bits = png_get_bit_depth(png_ptr, info_ptr);

  png_set_interlace_handling(png_ptr);
  png_read_update_info(png_ptr, info_ptr);

  /* read file */
  if (setjmp(png_jmpbuf(png_ptr))) {
    PRINT_ERROR("Error during read_image\n");
    exit(-1);
  }

  int bytes_per_row = png_get_rowbytes(png_ptr, info_ptr);  
  png_byte *png_buffer = new png_byte[this->Ny * bytes_per_row];
  
  row_pointers = (png_bytep*)png_malloc(png_ptr, this->Ny * sizeof(png_bytep));
  for (int y = 0; y < this->Ny; y++)
    row_pointers[y] = &png_buffer[y * bytes_per_row];  

  png_read_image(png_ptr, row_pointers);
  int num_channels = png_get_channels(png_ptr, info_ptr);
  fclose(fp);
  
  // Create data channels
  create_channels(num_channels, Nx, Ny);  

  switch (this->bits) {
    case 8:
      for (int ch = 0; ch < num_channels; ch++) {
        float *m = this->get_channel(ch);

        for (int y = 0; y < this->Ny; y++) {
          for (int x = 0; x < this->Nx; x++) {
            int pos = x * num_channels  + ch;
            m[this->Nx * y + x] = row_pointers[y][pos];
          }
        }
      }
      break;
    case 16:
      for (int ch = 0; ch < num_channels; ch++) {
        float *m = this->get_channel(ch);

        for (int y = 0; y < this->Ny; y++) {
          for (int x = 0; x < this->Nx; x++) {
            int pos = (num_channels * x + ch) * 2;
            m[this->Nx * y + x] = 256 * row_pointers[y][pos] + row_pointers[y][pos + 1];
          }
        }
      }
      break;
  }

  // Free memory
  free(row_pointers);
  delete[] png_buffer;
}

void CImage::save_PNG(char *filename, int bits_per_channel) {
  png_byte color_type;
  png_byte bit_depth;
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep * row_pointers;

  /* create file */
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    PRINT_ERROR("File %s could not be opened for writing\n", filename);
    exit(-1);
  }

  /* initialize stuff */
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr) {
    PRINT_ERROR("png_create_write_struct failed\n");
    exit(-1);
  }

  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    PRINT_ERROR("png_create_info_struct failed\n");
    exit(-1);
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    PRINT_ERROR("Error during init_io\n");
    exit(-1);
  }

  png_init_io(png_ptr, fp);

  /* write header */
  if (setjmp(png_jmpbuf(png_ptr))) {
    PRINT_ERROR("Error during writing header\n");
    exit(-1);
  }

  bit_depth = bits_per_channel;
  color_type = (this->get_num_channels() == 1 || bit_depth < 8 ?
                PNG_COLOR_TYPE_GRAY : PNG_COLOR_TYPE_RGB);

  int height = this->get_height();
  int width = this->get_width();
  
  png_set_IHDR(png_ptr, info_ptr, this->get_width(), height,
         bit_depth, color_type, PNG_INTERLACE_NONE,
         PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  png_write_info(png_ptr, info_ptr);

  // Create row pointers
  int bytes_per_row = png_get_rowbytes(png_ptr, info_ptr);  
  png_byte *png_buffer = new png_byte[this->Ny * bytes_per_row];
  
  row_pointers = (png_bytep*)png_malloc(png_ptr, this->Ny * sizeof(png_bytep));
  for (int y = 0; y < this->Ny; y++)
    row_pointers[y] = &png_buffer[y * bytes_per_row];

  /* write bytes */
  if (setjmp(png_jmpbuf(png_ptr))) {
    PRINT_ERROR("error during writing bytes\n");
    exit(-1);
  }

  int channels = this->get_num_channels();
  
  this->bits = bits_per_channel;

  if (bits_per_channel == 8) {
    for (int ch = 0; ch < channels; ch++) {
      float *orig = this->get_channel(ch);
      for (int y = 0; y < height; y++) {
        for (int x = 0; x < this->get_width(); x++)
          row_pointers[y][x * channels + ch] = orig[width * y + x];
      }
    }
  }

  if (bits_per_channel == 16) {
    float multiplier = (this->bits == 8 ? 256 : 1);
    for (int ch = 0; ch < channels; ch++) {
      float *orig = this->get_channel(ch);
      for (int y = 0; y < height; y++) {
        for (int x = 0; x < this->get_width(); x++) {
          int data = orig[width * y + x] * multiplier;
          int pos = (x * channels + ch) * 2;
          row_pointers[y][pos] = (data & 0xFF00) >> 8;
          row_pointers[y][pos + 1] = (data & 0x00FF);
        }
      }
    }
  }

  png_write_image(png_ptr, row_pointers);

  /* end write */
  if (setjmp(png_jmpbuf(png_ptr))) {
    PRINT_ERROR("error during end of write\n");
    exit(-1);
  }

  png_write_end(png_ptr, NULL);
  fclose(fp);

  /* cleanup heap allocation */
  // Free memory
  free(row_pointers);
  delete[] png_buffer;
}

void CImage::save_float_RGB(char *filename) {
  int num_channels;
  size_t dummy;

  FILE *outfile = fopen(filename, "wb");
  if (!outfile) {
    PRINT_ERROR("CImage: error creating float RGB file %s\n!", filename);
    exit(-1);
  }

  num_channels = this->get_num_channels();
  
  dummy = fwrite(&this->Nx, sizeof(this->Nx), 1, outfile);
  dummy += fwrite(&this->Ny, sizeof(this->Ny), 1, outfile);
  dummy += fwrite(&num_channels, sizeof(num_channels), 1, outfile);
  dummy = dummy; // To prevent the warning
  this->bits = this->get_bits_per_channel();
  PRINT_VERBOSE("CImage save RGB (%d, %d), %d channels, id %d\n",
    this->Nx, this->Ny, num_channels, this->id);

  for (int ch = 0; ch < num_channels; ch++) {
    float *m = this->get_channel(ch);
    dummy = fwrite(m, sizeof(float), this->Nx * this->Ny, outfile);
  }
  fclose(outfile);
}


int CImage::get_num_channels() {
  return this->num_channels;
}

int CImage::get_width() {
  return this->Nx;
}

int CImage::get_height() {
  return this->Ny;
}

int CImage::get_bits_per_channel() {
  return this->bits;
}

float *CImage::get_channel(int channel) {
  if (channel < 0 || channel > this->get_num_channels()) {
    PRINT_ERROR("invalid channel %d for image id %d\n", channel, this->id);
    exit(-1);
  }
  CFramework *fw = CFramework::get_framework();

  int channel_id = this->channels_id[channel];
  float *m = fw->find_matrix(channel_id)->data;  
  return m;
}
