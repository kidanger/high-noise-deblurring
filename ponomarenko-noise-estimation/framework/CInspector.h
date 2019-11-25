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

#ifndef _CINSPECTOR_H_
#define _CINSPECTOR_H_

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sys/time.h>
#include <sys/resource.h>
#include "CFramework.h"

#define START_INSPECTOR CInspector::inspector(fw, __FILE__, __LINE__);

#define INSPECTOR_SPACES " \t\r\n"

#define INSPECTOR_OK 0
#define INSPECTOR_ERR_BAD_ARGUMENTS -1
#define INSPECTOR_ERR_NO_MATRIX_ARGUMENT -2
#define INSPECTOR_ERR_MATRIX_DOES_NOT_EXIST -3
#define INSPECTOR_ERR_NO_MATRIX_SELECTED_YET -4
#define INSPECTOR_ERR_COORDINATES_ARGUMENT -5
#define INSPECTOR_ERR_OUT_OF_BOUNDS -6
#define INSPECTOR_ERR_NO_FILENAME -7
#define INSPECTOR_ERR_MEMORY_INFO -8
#define INSPECTOR_ERR_READ_ERROR -9
#define INSPECTOR_ERR_CAN_NOT_OPEN_FILE -10
#define INSPECTOR_ERR_INCOMPATIBLE_DIMENSIONS -11

using namespace std;

//! Data Inspector CLI Utility class
class CInspector {
  private:
    static MatrixInfo *G_using_matrix;
    static int *prod_dims;

    static inline string trim_right (const string & s, const string & t = " \t\r\n");
    static inline string trim_left (const string & s, const string & t = INSPECTOR_SPACES);
    static inline string trim (const string & s, const string & t = INSPECTOR_SPACES);
    static int read_tokens(vector<string> *tokens, string input);

    static void describe_matrix(MatrixInfo *matrix);
    static bool check_command_arguments(string command, int num_tokens);

    static int get_pos_by_coordinates(int *x, MatrixInfo *matrix);

    static int process_input(string input, CFramework *fw);
    static void process_error(int status);


  public:
    //! Starts the Inspector CLI
    /*!
      \param *fw Pointer to the framework object
      \param *file String to the __FILE__ in which the execution was interrupted
      \param line __LINE__ in which the execution was interrupted
    */
    static void inspector(CFramework *fw, const char *file, int line);
};

#endif
