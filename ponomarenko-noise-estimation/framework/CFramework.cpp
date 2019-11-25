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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <execinfo.h>
#include "CFramework.h"

bool CFramework::CFramework_verbose;
CFramework *CFramework::CFramework_ref = NULL;

CFramework::CFramework() {
  PRINT_VERBOSE("CFramework created\n");
  CFramework::CFramework_verbose = false;
  this->unique_index_counter = 1;
  CFramework::CFramework_ref = this;
}

CFramework::~CFramework() {
  PRINT_VERBOSE("Cleaning up %d matrices\n",
    (int)this->mmatrices.size());
  //
  if (CFramework_verbose)
    this->print_matrix_list();
  //
  while (this->mmatrices.size() != 0) {
    map<int, MatrixInfo*>::iterator it = this->mmatrices.begin();
    MatrixInfo *matrix_info = it->second;
    this->delete_matrix(matrix_info->id);
  }
  
  // Delete singleton framework object
  if (CFramework::CFramework_ref != NULL) {
   if (CFramework_verbose)
     PRINT_VERBOSE("Deleting singleton framework object\n");
     
   CFramework::CFramework_ref = NULL;  
   delete CFramework::CFramework_ref;
  }
  
}

CFramework* CFramework::get_framework() {
  if (CFramework::CFramework_ref == NULL) {
    CFramework::CFramework_ref = new CFramework();
    if (CFramework_verbose)
      PRINT_VERBOSE("Creating singleton framework object\n");
  }
  return CFramework::CFramework_ref;
}

void CFramework::print_matrix_list() {
  printf("\nFramework matrix list (%d matrices):\n",
    (int)this->mmatrices.size());
  printf("--------------------------------------------\n");
  map<int, MatrixInfo*>::iterator end = this->mmatrices.end(); 
  int i = 1;
  for (map<int, MatrixInfo*>::iterator it = this->mmatrices.begin();
       it != end; ++it) {
    MatrixInfo *matrix_info = it->second;
    printf("%d) Matrix M_%d data at %p\n",
      i, matrix_info->id, matrix_info->data);
    i++;
  }
  printf("\n");
}

void CFramework::set_verbose(bool verbose) {
  CFramework::CFramework_verbose = verbose;
}

float *CFramework::create_matrix(int num_dims, int *N, int &out_id) {
  int product_dims = 1;
  for (int i = 0; i < num_dims; i++)
    product_dims *= N[i];

  float *data = new float[product_dims];
  int id = this->unique_index_counter++;
  MatrixInfo *matrix_info = new MatrixInfo(id, num_dims, N, data);
  this->mmatrices[id] = matrix_info;
  out_id = id;

  PRINT_VERBOSE("Created matrix M_%d at %p, data at %p\n",
    id, matrix_info, matrix_info->data);
  return data;
}

float *CFramework::create_matrix(int num_dims, int *N) {
  int dummy;
  return this->create_matrix(num_dims, N, dummy);
}


float *CFramework::create_array(int N, int &id) {
  int N_dims[1] = {N};
  return this->create_matrix(1, N_dims, id);
}

float *CFramework::create_array(int N) {
  int dummy;
  return this->create_array(N, dummy);
}

void CFramework::delete_matrix(int id) {
  MatrixInfo *matrix_info = CFramework::find_matrix(id);
  if (matrix_info == NULL) {
    PRINT_ERROR("matrix M_%d not in memory!\n", id);
    exit(-1);
  }
  PRINT_VERBOSE("Deleting matrix M_%d at %p, data at %p\n",
    id, matrix_info, matrix_info->data);
  delete[] matrix_info->data;
  this->mmatrices.erase(id);
  delete matrix_info;
}

MatrixInfo *CFramework::find_matrix(int id) {
  if (this->mmatrices.find(id) != this->mmatrices.end())
    return this->mmatrices[id];
  else
    return NULL;
}

int CFramework::find_matrix(float *data) {
  map<int, MatrixInfo*>::iterator end = this->mmatrices.end(); 
  for (map<int, MatrixInfo*>::iterator it = this->mmatrices.begin();
       it != end; ++it) {
    MatrixInfo *matrix_info = it->second;
    if (matrix_info->data == data)
      return matrix_info->id;
  }
  return -1;
}
// Used mainly by the Inspector.
MatrixInfo *CFramework::find_matrix_by_list_pos(int pos) {
  map<int, MatrixInfo*>::iterator end = this->mmatrices.end(); 
  int i = 1;
  for (map<int, MatrixInfo*>::iterator it = this->mmatrices.begin();
       it != end; ++it) {
    MatrixInfo *matrix_info = it->second;
    if (i == pos)
      return matrix_info;
    i++;
  }
  return NULL;
}

void CFramework::save_matrix(MatrixInfo *matrix, char *filename) {
  FILE *fp = fopen(filename, "wb");

  int num_elements = 1;
  for (int i = 0; i < matrix->num_dims; i++)
    num_elements *= matrix->N[i];

  fwrite(&matrix->num_dims, sizeof(matrix->num_dims), 1, fp);
  for (int i = 0; i < matrix->num_dims; i++)
    fwrite(&matrix->N[i], sizeof(matrix->N[0]), 1, fp);

  size_t written = fwrite(matrix->data,
                          sizeof(matrix->data[0]),
                          num_elements, fp);
  if (written != num_elements) {
    fprintf(stderr, "Write error!\n");
    exit(-1);
  }

  fclose(fp);
}

int CFramework::load_matrix(MatrixInfo *matrix, char *filename) {
  FILE *fp = fopen(filename, "rb");
  if (fp == NULL)
    return -1;

  // Read number of dimensions
  int num_dims; 
  size_t dummy = fread(&num_dims, sizeof(matrix->num_dims), 1, fp);

  // Check number of dimensions
  if (num_dims != matrix->num_dims) {
    fclose(fp);
    return -2;
  }

  // Read matrix size
  dummy = fread(matrix->N, sizeof(matrix->N[0]), num_dims, fp);

  // Get number of elements
  int num_elements = 1;
  for (int i = 0; i < num_dims; i++)
    num_elements *= matrix->N[i];

  // Get matrix name and delete it
  delete[] matrix->data;
  matrix->data = new float[num_elements];

  // Read data
  size_t read = fread(matrix->data,
                      sizeof(matrix->N[0]),
                      num_elements, fp);
  if (read != num_elements) {
    fprintf(stderr, "Read error!\n");
    return -3;
  }
  
  fclose(fp);
  return 0;
}

// Source:
// http://tombarta.wordpress.com/2008/08/01/c-stack-traces-with-gcc/
void CFramework::print_backtrace(FILE *out) {
    const size_t max_depth = 100;
    size_t stack_depth;
    void *stack_addrs[max_depth];
    char **stack_strings;

    stack_depth = backtrace(stack_addrs, max_depth);
    stack_strings = backtrace_symbols(stack_addrs, stack_depth);

    for (size_t i = 1; i < stack_depth; i++)
        fprintf(out, "    %s\n", stack_strings[i]);

    free(stack_strings); // malloc()ed by backtrace_symbols
    fflush(out);
}
