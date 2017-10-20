/*
 * IO.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: niber
 */

#include "IO.h"
#include <sstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>

namespace io {

int do_mkdir(char * path, mode_t mode)
{
	struct stat st;
	int status = 0;

	if (stat(path, &st) != 0) {
		if (::mkdir(path, mode) != 0 && errno != EEXIST)
			status = -1;
	} else if (!S_ISDIR(st.st_mode)) {
		errno = ENOTDIR;
		status = -1;
	}

	return status;
}

void mkdir(const char * folder)
{
	/* Create output folder structure */
	std::string folder_tmp(folder);
	char * foldercopy = const_cast<char *>(folder_tmp.c_str());
	char * start, * end;

	int status = 0;
	start = foldercopy;

	while(status == 0 && (end = strchr(start, '/')) != 0) {
		if(start != end) {
			*end = '\0';
			status = do_mkdir(foldercopy, 0777);
			*end = '/';
		}

		start = end + 1;
	}

	if(status == 0)
		do_mkdir(foldercopy, 0777);
}

} /* namespace io */
