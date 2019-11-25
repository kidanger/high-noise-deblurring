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

#ifndef _CFRAMEWORK_H_
#define _CFRAMEWORK_H_

#include <map>
#include <string>

#define PRINT_VERBOSE(...) if (CFramework::CFramework_verbose) printf(__VA_ARGS__)
#define PRINT_ERROR(...) fprintf(stderr, "Error in line %d at %s: ", __LINE__, __FILE__); fprintf(stderr, __VA_ARGS__)

using namespace std;

//! Information container for matrices
class MatrixInfo {
  public:
    int id;
    int num_dims;
    int *N; // Array of dimensions
    float *data;
    MatrixInfo(int id, int num_dims, int *N, float *data) {
      this->id = id;
      this->num_dims = num_dims;
      //
      this->N = new int[num_dims];
      for (int i = 0; i < num_dims; i++) {
        this->N[i] = N[i];
      }
      //
      this->data = data;
    }

    virtual ~MatrixInfo() {
      delete[] this->N;
    }
};

//! Main framework class
class CFramework {
  private:
    map <int, MatrixInfo*> mmatrices;
    int unique_index_counter;
    static CFramework *CFramework_ref; // Singleton reference
  protected:
    CFramework(); // Singleton constructor
  public:
    static bool CFramework_verbose;
    
    virtual ~CFramework();

    //! Return a reference to the framework.
    /*!
      \return pointer to the existing or created framework
    */
    static CFramework *get_framework();

    //! Creates a matrix
    /*!
      \param num_dims Number of dimensions
      \param *N Pointer to the array that contains the dimensions
      \return pointer to the matrix data
    */
    float *create_matrix(int num_dims, int *N);
    float *create_matrix(int num_dims, int *N, int &out_id);

    //! Creates a 1D matrix
    /*!
      \param N Length of the array
      \return pointer to the matrix data
    */
    float *create_array(int N);
    float *create_array(int N, int &out_id);

    //! Deletes a matrix
    /*!
      \param name Name of the matrix
    */
    void delete_matrix(int id);

    //! Sets the verbose level
    /*!
      \param verbose Boolean to activate or deactivate the verbosity
    */
    void set_verbose(bool verbose);

    //! Prints the matrix list
    /*!
    */
    void print_matrix_list();

    //! Looks for the matrix with the given name
    /*!
      \param name Name of the matrix
      \return information object for the matrix or NULL if not found
    */
    MatrixInfo *find_matrix(int id);
    
    //! Looks for the matrix with the given data pointer
    /*!
      \param *data Data pointer
      \return matrix ID of -1 if not found
    */
    int find_matrix(float *data);


    //! Returns the matrix in the specified position in the matrix list
    /*!
      \param pos Index (1..N) in the matrix list.
      \return information object for the matrix or NULL if not found
    */
    MatrixInfo *find_matrix_by_list_pos(int pos);

    //! Saves the contents of the matrix to a file
    /*!
      \param matrix Information object for the matrix
      \param filename Name of the file to write the data to
    */
    void save_matrix(MatrixInfo *matrix, char *filename);

    //! Loads data from a file into an existing matrix
    /*!
      \param matrix Information object for the matrix
      \param filename Name of the file to read the data from
      \return Status (0=OK)
    */
    int load_matrix(MatrixInfo *matrix, char *filename);

    //! Prints a backtrace of the call stack
    /*!
      \param out File descriptor to write the backtrace to
    */
    void print_backtrace(FILE *out);
};

#endif
