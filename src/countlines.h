#pragma once



#include <iostream>
#include <fstream> 
#include "armaheader.h"


using namespace std;

inline size_t countlines(std::string fn);


// count lines of a file
// [[Rcpp::export]]
size_t countlines(std::string fn) {
    std::ifstream in_file(fn);
    size_t nlines =  (size_t) std::count(std::istreambuf_iterator<char>(in_file), std::istreambuf_iterator<char>(), '\n');

    // nlines needs be corrected when the last line is not ended with a newline
    in_file.seekg(-1, in_file.end);
    unsigned char lastchar = in_file.get();
    if(lastchar != '\n') {
        ++nlines;
    }
    return nlines;
}







