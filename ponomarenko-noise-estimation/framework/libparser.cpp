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

/*
 *  libparser.cpp
 *  wxSource.2010.2
 *
 *  Created by Antoni Buades on 17/06/2010.
 *
 */

#include "libparser.h"




void  printusage(const char *pname,
		 const char *gp,
		 vector<OptStruct*>  &opt,
		 vector<ParStruct*>  &par)
{
	
	//int nopt = opt.size();
	int npar = par.size();
	
	///// USAGE
	printf("\nusage: %s ", pname);
	for(int i=0; i < (int) strlen(gp); i++)
		if (gp[i] != ':')
		{
			printf("[-%c",gp[i]);
			
			if (i+1 < (int) strlen(gp) && gp[i+1] ==  ':')
			  printf(" %c] ", gp[i]);
			else
			  printf("] ");
			
		}
	
	for(int i=0; i < npar; i++)
		printf(" %s ", par[i]->name);
	
	printf("\n");
	//// PARAMETERS
	
	int j=0;
	for(int i=0; i < (int) strlen(gp); i++)
		if (gp[i] != ':')
		{
			printf("\t-%c",gp[i]);
			
			if (i+1 < (int) strlen(gp) && gp[i+1] ==  ':') {
				
				printf("  %c\t %s ", gp[i], opt[j]->comment);
				if (opt[j]->defvalue != NULL)
				  printf("(Default: %s)",opt[j]->defvalue);
				
				printf("\n");
				
			}		
			else printf("\t %s \n", opt[j]->comment);
			
			j++;
		}
	
	
	for(int i=0; i < npar; i++)
	{
		printf("\t%s",par[i]->name);
		printf("\t %s\n", par[i]->comment);
	}
	
	
}




int parsecmdline(const char *pname,
		 const char *function,
		 const int argc, char **argv,
		 vector <OptStruct*> & opt,
		 vector <ParStruct*> & par) {
	int nopt = opt.size();
	int npar = par.size();
	
	
	char gp[2*nopt+1];
	gp[0]='\0';
	
	
	for(int i=0; i < nopt; i++) {
		opt[i]->flag = 0;
		opt[i]->value=NULL;
		strcat(gp, opt[i]->gp);
	}
	
	for(int i=0; i < npar; i++)
		par[i]->value = NULL;
	
	opterr = 0;	// No messages by getopt
	
	int c;
	while ((c = getopt (argc, argv, gp)) != -1)
	{
		
		int j=0;
		for(unsigned int i=0; i < strlen(gp); i++)
			if (c == gp[i])
			{
				
				opt[j]->flag = 1;
				/*				if (optarg != NULL && optarg[0] == '-')
				 {	 
				 printf("\n%s: %s\n", pname, function);
				 fprintf (stderr,
				   "\nerror: option -%c requires an argument.\n", c);
				 printusage(pname,function, gp, argc, argv,
				            opt, nopt, par, npar);
				 return 0;
				 
				 }*/
				
				opt[j]->value = optarg;
				break;
				
			} else if (gp[i] != ':') j++;
		
		
		
		if (c == '?')
		{	
			
			unsigned int i = 0;
			for(i=0; i < strlen(gp); i++)
				if (optopt == gp[i])
				{
					printf("\n%s: %s\n", pname, function);
					fprintf (stderr,
					  "\nerror: option -%c requires an argument.\n",
					  optopt);
					break;	
				}
			
			if (i == strlen(gp)) { 	printf("\n%s: %s\n", pname, function);
				fprintf (stderr, "\nerror: unknown option `-%c'.\n",
				  optopt);
			}
			
			printusage(pname, gp,  opt,  par);
			return 0;
			
		}
		
	}
	
	
	//// Setting default values for non selected options
	for (int j=0; j < nopt; j++)
		if (opt[j]->flag == 0 && opt[j]->defvalue != NULL)
		  opt[j]->value =  opt[j]->defvalue;
	
	if (argc - optind != npar) {
		printf("\n%s: %s\n", pname, function);
		fprintf (stderr, "\nerror: incorrect number of parameters\n");
		printusage(pname, gp,  opt,par);
		return 0;
	}
	
	int i=0;
	for (int index = optind; index < argc ; index++, i++)
		par[i]->value = argv[index];
	
	return 1;
}

