/*
 * ncols.cpp
 *
 *  Created on: Jul 26, 2014
 *      Author: kaiyin
 */



#include "ncols.h"

// Counts the number of columns of tab / space delimited file.
// Quoted fields not supported, assumes that no field contains any whitespace!

size_t ncols(std::string fn) {
    try {
        fileExists(fn);
        std::ifstream in_file(fn.c_str());
        std::string tmpline;
        std::getline(in_file, tmpline);
        std::vector<std::string> strs;
        size_t numCols = 0;
        std::istringstream lineStream(tmpline);
        std::string tmpword = "";
        while(lineStream >> tmpword) {
        	if(not tmpword.empty()) {
//        		Rcpp::Rcout << "I got word: " << tmpword << "\n";
        		std::cout << "I got word: " << tmpword << "\n";
        		numCols++;
        	}
        }

        return numCols;
    } catch (const std::string& e) {
        std::cerr << "\n" << e << "\n";
        exit(EXIT_FAILURE);
    }
}
