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

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CInspector.h"
#include "CFramework.h"

MatrixInfo *CInspector::G_using_matrix = NULL;
int *CInspector::prod_dims = NULL;

// Credit for the trim functions:
// http://snippets.dzone.com/posts/show/7964
inline string CInspector::trim_right (const string & s, const string & t)
{
    string d (s);
    string::size_type i (d.find_last_not_of (t));
    if (i == string::npos)
        return "";
    else
        return d.erase (d.find_last_not_of (t) + 1) ;
}  
//
inline string CInspector::trim_left (const string & s, const string & t)
{
    string d (s);
    return d.erase (0, s.find_first_not_of (t)) ;
}  
//
inline string CInspector::trim (const string & s, const string & t)
{
    string d (s);
    return trim_left (trim_right (d, t), t) ;
}

/////

int CInspector::read_tokens(vector<string> *tokens, string input) {
  tokens->clear();

  istringstream iss(input);
  do {
    string sub;
    iss >> sub;
    if (sub.c_str() != NULL && sub != "")
      tokens->push_back(sub);
  } while (iss);

  return tokens->size();
}

void list_of_ints(char *dest, const char *separator, int *numbers, int N) {
  char buffer[1024];
  strcpy(dest, "");
  for (int i = 0; i < N; i++) {
    sprintf(buffer, "%d", numbers[i]);
    strcat(dest, buffer);
    if (i != N-1)
      strcat(dest, separator);
  }
}


void CInspector::describe_matrix(MatrixInfo *matrix) {
  if (matrix == NULL) {
    printf("No matrix chosen yet.\n");
    return;
  }
  //
  char dims_str[1024];
  list_of_ints(dims_str, " x ", matrix->N, matrix->num_dims);
  printf("Matrix id=%d, data at %p, %d dims: (%s).\n", matrix->id, matrix->data, matrix->num_dims, dims_str);
}

bool CInspector::check_command_arguments(string command, int num_tokens) {
  if (command == "list" && num_tokens != 1) return false;
  if (command == "use" && num_tokens != 2) return false;
  return true; // OK! :-)
}


void acc_product_of_dims(MatrixInfo *matrix, int *P) {
  P[0] = 1;
  for (int i = 1; i < matrix->num_dims; i++) {
    P[i] = matrix->N[i] * P[i-1];
  }
}

int CInspector::get_pos_by_coordinates(int *x, MatrixInfo *matrix) {
  int pos = 0;
  for (int i = 0; i < matrix->num_dims; i++) {
    pos += CInspector::prod_dims[i] * x[i];
  }

  return pos;
}

int CInspector::process_input(string input, CFramework *fw) {
  vector<string> tokens;

  int num_tokens = read_tokens(&tokens, input);
  if (num_tokens == 0) return INSPECTOR_OK; // just pressed <ENTER>

  string command = tokens[0];
  if (!check_command_arguments(command, num_tokens)) {
    return INSPECTOR_ERR_BAD_ARGUMENTS;
  }

  if (command == "l" || command == "list") {
    fw->print_matrix_list();
    return INSPECTOR_OK;
  }

  if (command == "u" || command == "use") {
    if (tokens.size() != 2) {
      return INSPECTOR_ERR_NO_MATRIX_ARGUMENT;
    }
    string name = tokens[1];

    MatrixInfo *matrix = fw->find_matrix_by_list_pos(atoi(name.c_str()));
    //MatrixInfo *matrix = fw->find_matrix(name);
    if (matrix == NULL) {
      return INSPECTOR_ERR_MATRIX_DOES_NOT_EXIST;
    }
    //
    G_using_matrix = matrix;
    describe_matrix(matrix);

    CInspector::prod_dims = new int[matrix->num_dims];
    acc_product_of_dims(matrix, CInspector::prod_dims);
  }

  if (command == "bt") {
    printf("Backtrace:\n");
    fw->print_backtrace(stdout);
    return INSPECTOR_OK;
  }


  if (command == "desc") {
    describe_matrix(G_using_matrix);
    return INSPECTOR_OK;
  }

  if (command == "g" || command == "get") {
    if (G_using_matrix == NULL) {
      return INSPECTOR_ERR_NO_MATRIX_SELECTED_YET;
    }

    // Check number of coordinates
    int num_coords = tokens.size() - 1;
    if (num_coords != G_using_matrix->num_dims) {
      return INSPECTOR_ERR_COORDINATES_ARGUMENT;
    }

    int coords[num_coords];
    for (int i = 0; i < num_coords; i++) {
      int coord = atoi(tokens[i+1].c_str());
      if (coord >= G_using_matrix->N[i]) {
        return INSPECTOR_ERR_OUT_OF_BOUNDS;
      }
      coords[i] = coord;
    }
    //
    int pos = get_pos_by_coordinates(coords, G_using_matrix);

    char coords_buffer[1024];
    list_of_ints(coords_buffer, ", ", coords, num_coords);


    printf("Matrix_%d[%s] == %f (pos %d).\n", G_using_matrix->id, coords_buffer, G_using_matrix->data[pos], pos);
    return INSPECTOR_OK;
  }

  if (command == "s" || command == "set") {
    if (G_using_matrix == NULL) {
      return INSPECTOR_ERR_NO_MATRIX_SELECTED_YET;
    }

    // Check number of coordinates
    int num_coords = tokens.size() - 2;
    if (num_coords != G_using_matrix->num_dims) {
      return INSPECTOR_ERR_COORDINATES_ARGUMENT;
    }

    int coords[num_coords];
    for (int i = 0; i < num_coords; i++) {
      int coord = atoi(tokens[i+1].c_str());
      if (coord >= G_using_matrix->N[i]) {
        return INSPECTOR_ERR_OUT_OF_BOUNDS;
      }
      coords[i] = coord;
    }
    //
    int pos = get_pos_by_coordinates(coords, G_using_matrix);
    float new_value = atof(tokens[num_coords+1].c_str());

    char coords_buffer[1024];
    list_of_ints(coords_buffer, ", ", coords, num_coords);

    printf("Before: M_%d[%s] == %f (pos %d).\n", G_using_matrix->id, coords_buffer, G_using_matrix->data[pos], pos);
    G_using_matrix->data[pos] = new_value;
    printf("After: M_%d[%s] == %f (pos %d).\n", G_using_matrix->id, coords_buffer, G_using_matrix->data[pos], pos);
    return INSPECTOR_OK;
  }

  if (command == "mem") {
    rusage usage;
    int status = getrusage(RUSAGE_SELF, &usage);
    if (status == 0) {
      printf("maximum resident set size: %lu\n", usage.ru_maxrss);
      printf("integral shared memory size: %lu\n", usage.ru_ixrss);
      printf("integral unshared data size: %lu\n", usage.ru_idrss);
      printf("integral unshared stack size: %lu\n", usage.ru_isrss);
      printf("page reclaims (soft page faults): %lu\n", usage.ru_minflt);
      printf("page faults (hard page faults): %lu\n", usage.ru_majflt);
      printf("swaps: %lu\n", usage.ru_nswap);
      printf("block input operations: %lu\n", usage.ru_inblock);
      printf("block output operations: %lu\n", usage.ru_oublock);
      printf("IPC messages sent: %lu\n", usage.ru_msgsnd);
      printf("IPC messages received: %lu\n", usage.ru_msgrcv);
      printf("signals received: %lu\n", usage.ru_nsignals);
      printf("voluntary context switches: %lu\n", usage.ru_nvcsw);
      printf("involuntary context switches: %lu\n", usage.ru_nivcsw);
    }
    return INSPECTOR_OK;
  }

  if (command == "save") {
    if (G_using_matrix == NULL)
      return INSPECTOR_ERR_NO_MATRIX_SELECTED_YET;
    if (tokens.size() != 2)
      return INSPECTOR_ERR_NO_FILENAME;
    //
    char *filename = (char*)tokens[1].c_str();
    int id = G_using_matrix->id;
    fw->save_matrix(G_using_matrix, filename);
    G_using_matrix = fw->find_matrix(id);
    return INSPECTOR_OK;
  }

  if (command == "load") {
    if (G_using_matrix == NULL)
      return INSPECTOR_ERR_NO_MATRIX_SELECTED_YET;
    if (tokens.size() != 2)
      return INSPECTOR_ERR_NO_FILENAME;
    //
    char *filename = (char*)tokens[1].c_str();
    int status = fw->load_matrix(G_using_matrix, filename);
    switch (status) {
      case  0: return INSPECTOR_OK; break;
      case -1: return INSPECTOR_ERR_CAN_NOT_OPEN_FILE; break;
      case -2: return INSPECTOR_ERR_INCOMPATIBLE_DIMENSIONS; break;
      case -3: return INSPECTOR_ERR_READ_ERROR; break;
    }
    return INSPECTOR_OK;
  }


  if (command == "h" || command == "help") {
    printf("Framework CLI inpector. Help not available yet, sorry.\n");
    return INSPECTOR_OK;
  }

  return INSPECTOR_OK;
}

void CInspector::process_error(int err_code) {
  const char *msg = "OK";
  switch (err_code) {
    case INSPECTOR_ERR_BAD_ARGUMENTS:
      msg = "incorrect number of arguments";
      break;
    case INSPECTOR_ERR_NO_MATRIX_ARGUMENT:
      msg = "no matrix identifier specified";
      break;
    case INSPECTOR_ERR_MATRIX_DOES_NOT_EXIST:
      msg = "matrix doesn't exist";
      break;
    case INSPECTOR_ERR_NO_MATRIX_SELECTED_YET:
      msg = "no matrix selected yet";
      break;
    case INSPECTOR_ERR_COORDINATES_ARGUMENT:
      msg = "wrong number of coordinates for the selected matrix";
      break;
    case INSPECTOR_ERR_OUT_OF_BOUNDS:
      msg = "out of bounds";
      break;
    case INSPECTOR_ERR_NO_FILENAME:
      msg = "no filename specified";
      break;
    case INSPECTOR_ERR_MEMORY_INFO:
      msg = "can't get memory information";
      break;
    case INSPECTOR_ERR_READ_ERROR:
      msg = "read error";
      break;
    case INSPECTOR_ERR_INCOMPATIBLE_DIMENSIONS:
      msg = "incompatible dimensions";
      break;
    case INSPECTOR_ERR_CAN_NOT_OPEN_FILE:
      msg = "can't open file";
      break;
    default:
      msg = "Unknown error";
  }

  printf("Error %d: %s.\n", err_code, msg);
}

void CInspector::inspector(CFramework *fw, const char *file, int line) {
  printf("\nCLI Inpector called at file %s:%d\n", file, line);
  printf("(C) 2011 Miguel Colom.\n");

  string PROMPT = "> ";
  char buffer[4096];

  string input = "";
  while (input != "q" && input != "quit" && input != "exit") {
    cout << PROMPT;
    cin.getline(buffer, 4096);
    input = string(buffer);
    input = trim(input);
    int err_code = process_input(input, fw);
    if (err_code != INSPECTOR_OK)
      process_error(err_code);
  }
  printf("\nCLI Inpector exiting.\nResuming execution from %s:%d.\n\n", file, line);
}

