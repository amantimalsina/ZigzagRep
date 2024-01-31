#ifndef _ZZREP_H_
#define _ZZREP_H_

#include <vector>
#include <tuple>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <memory>

namespace ZZREP { 
/* HELPER DATA STRUCTURES: */

//typedef int Integer;
//typedef std::vector<Integer> Simplex;

class ZigzagRep {
public:
    /*
      Description: 'filt_simp' and 'filt_op' should have the same length which altogether
      specify the input zigzag filtration. 'filt_simp' specifies the simplices
      being added or deleted (following the order of the filtration) and 
      'filt_op' specifies whether it's an addition (true) or deletion (false).
      'persistence' returns the barcode, with the first element of the tuple
      being the birth, the second element being the death, the third
      being the dimension, along with the representative. Here, m denotes the dimension of the highest dimensional simplex from the filtration.
    */
    void compute(
        const std::vector<std::vector<int> > &filt_simp, 
        const std::vector<bool> &filt_op,
        std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > *persistence,
        std::vector< std::vector<int>> *id_to_i,
        const int m);
};

}
#endif