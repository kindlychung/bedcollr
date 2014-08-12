/*
 * readcols.cpp
 *
 *  Created on: Aug 11, 2014
 *      Author: kaiyin
 */

#include "readcols.h"


// [[Rcpp::export]]
Rcpp::CharacterMatrix readcols(std::string fn,
		std::vector<int> colsel, size_t nFirstSkipLines,
		size_t nSkipUnit) {

	if (colsel.empty()) {
		throw std::string("You didn't select any column!");
	}


	size_t nc_file = ncols(fn);
	std::cout << "File has " << nc_file << " columns \n";
	size_t nr_file = countlines(fn);
	std::cout << "File has " << nr_file << " rows \n";
	std::cout << "You want to skip the first " <<  nFirstSkipLines << " lines \n";
	size_t nr_left = nr_file - nFirstSkipLines;
	std::cout << "Lines left: " << nr_left << "\n";

	std::cout << "Of the rest you want to select one line out of every " << nSkipUnit << " lines \n";
	size_t nr = (size_t) (nr_left / nSkipUnit);
	std::cout << "You selected " << nr << " rows \n";
	size_t nc = colsel.size();
	std::cout << "You selected " << nc << " columns \n";


	{
		size_t remainder = (size_t) ((nr_file - nFirstSkipLines) % nSkipUnit);
		if (remainder != 0) {
			std::cerr << "Number of lines to read is not a multiple of nSkipUnit! \n";
			nr++;
		}
		else {
			std::cout << "Number of lines to read is a multiple of nSkipUnit \n";
		}
	}

	std::cout << "Calculating max column number...\n";
	unsigned int colsel_max = *std::max_element(colsel.begin(), colsel.end());
	std::cout << "Max column number is " << colsel_max << "\n";

	if (colsel_max > nc_file) {
		throw std::string("Some col number(s) are out of range!");
	}

	// c++ is 0-based, adjust for it
	for (unsigned int i = 0; i < colsel.size(); i++) {
		colsel[i]--;
	}

	// initialize a 2d vector (matrix) with fixed size
//	std::vector< std::vector<std::string> > res(nc, std::vector<std::string>(nr));
    std::cout << "Allocating a matrix of " << nr << " rows and " << nc << " columns.\n";
	Rcpp::CharacterMatrix res(nr, nc);
	std::ifstream infile(fn.c_str());
	std::string tmpline;

	// skip lines in the beginning
	for (int lineIter = 0; lineIter < nFirstSkipLines; lineIter++) {
		std::cout << "Skipping line " << lineIter + 1 << "\n";
		getline(infile, tmpline);
	}


	size_t rowIter = 0;
	for (size_t lineIter = 0; lineIter < nr_left; lineIter++) {
		std::string line;
		getline(infile, line);
		if (lineIter % nSkipUnit == 0) {
			std::istringstream lineStream(line);
			size_t colIter = 0;
			for (size_t wordIter = 0; wordIter <= colsel_max; wordIter++) {
				std::string tmpword;
				lineStream >> tmpword;
				if (std::find(colsel.begin(), colsel.end(), wordIter)
						!= colsel.end()) {
					res(rowIter, colIter) = tmpword;
					colIter++;
				}
			}
			rowIter++;
		}
	}

    return res;

}
