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

#ifndef _LIBPARSER_H_
#define _LIBPARSER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <getopt.h>

using namespace std;
#include <vector>
#include <list>


typedef struct optstruct
{	
	const char *gp;
	int flag;
	const char *defvalue;
	const char *value;
	const char *comment;
	
} OptStruct;

typedef struct parstruct
{
	const char * name;
	const char * value;
	const char * comment;
} ParStruct;

int parsecmdline(const char *pname,
		 const char *function,
		 int argc, char **argv,
		 vector <OptStruct*> & opt,
		 vector <ParStruct*> & par);

#endif

