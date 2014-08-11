/*
 * readcols.h
 *
 *  Created on: Aug 11, 2014
 *      Author: kaiyin
 */

#ifndef READCOLS_H_
#define READCOLS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "ncols.h"
#include "countlines.h"
#include "armaheader.h"

Rcpp::CharacterMatrix readcols(std::string fn,
		std::vector<unsigned int> colsel, size_t nFirstSkipLines,
		size_t nSkipUnit);

#endif /* READCOLS_H_ */
