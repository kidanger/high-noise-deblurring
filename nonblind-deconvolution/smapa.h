#pragma once

#include <stdlib.h>

static const char* int_fmt = "%d";
static const char* float_fmt = "%f";
static const char* double_fmt = "%lf";

#define SMART_PARAMETER(t,n,v) static t& n(void)\
{\
	static bool smapa_known_ ## n = false;\
	static t smapa_value_ ## n = v;\
	if (!smapa_known_ ## n)\
	{\
		char *sv = getenv(#n);\
		t y;\
		if (sv && sscanf(sv, t##_fmt, &y) == 1)\
		{\
			fprintf(stderr, "%s = ", #n);\
			fprintf(stderr, t##_fmt, y); \
			fprintf(stderr, "\n"); \
			smapa_value_ ## n = y;\
		}\
		smapa_known_ ## n = true;\
	}\
	return smapa_value_ ## n;\
}

#define SMART_PARAMETER_STR(n) static const char*& n(void)\
{\
	static bool smapa_known_ ## n = false;\
	static const char* smapa_value_ ## n = 0;\
	if (!smapa_known_ ## n)\
	{\
		char *sv = getenv(#n);\
		if (sv && *sv) {\
			fprintf(stderr, "%s = %s\n", #n, sv);\
			smapa_value_ ## n = sv;\
		}\
		smapa_known_ ## n = true;\
	}\
	return smapa_value_ ## n;\
}


