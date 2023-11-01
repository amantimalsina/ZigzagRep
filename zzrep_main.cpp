#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <tuple>

#include "zzrep.h"

void getFilePurename(const std::string& filename, std::string *purename) {
    auto pos = filename.find_last_of("/\\");
    if (pos != std::string::npos) {
        *purename = filename.substr(pos + 1);
    } else {
        *purename = filename;
    }
    
    pos = purename->find_last_of(".");
    *purename = purename->substr(0, pos);
}

void parseSimplex(const std::string& str, char &op, std::vector<int> &simp) {
    std::istringstream iss(str);
    iss >> op;

    int index;
    while (iss >> index) { simp.push_back(index); }
}

int main(const int argc, const char *argv[]) {
    if (argc < 2) 
    { std::cerr << "Err: input not large enough" << std::endl; return -1; }

    const std::string infilename(argv[1]);
    std::ifstream filt_fin(infilename);

    if (!filt_fin)
    { std::cerr << "Err: input file open failed" << std::endl; return -1; }

    std::string line;
    char op;
    std::vector<std::vector<int> > filt_simp;
    std::vector<bool> filt_op;
    int m = 0;

    while (filt_fin) {
        std::getline(filt_fin, line);
        if (line.size() == 0) { continue; }

        filt_simp.emplace_back();
        parseSimplex(line, op, filt_simp.back());
        m = std::max(m, (int)filt_simp.back().size() - 1);

        if (op == 'i') {
            filt_op.push_back(true);
        } else {
            //assert(op == 'd');
            filt_op.push_back(false);
        }
    }


    filt_fin.close();
    filt_simp.pop_back();
    filt_op.pop_back();
    
    // Let's measure the time it takes to compute the zigzag rep:
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > persistence;
    std::vector <std::map<int, int>> i_to_id;

    ZZREP::ZigzagRep zzr;
    zzr.compute(
        filt_simp, 
        filt_op,
        &persistence,
        &i_to_id,
        m);
    std::string purename;
    getFilePurename(infilename, &purename);
    std::ofstream pers_fout(purename + "_pers");
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
        << "[ms]" << std::endl;

    // Change this to add the representatives to the file.
    for (const auto& e : persistence) {
        pers_fout << std::get<2>(e) << "-dimensional bar [" << std::get<0>(e) << ", " << std::get<1>(e) << "]" << std::endl;   
        pers_fout << "Representatives: " << std::endl;
        for (auto i : std::get<3>(e)) {
            pers_fout << "From " << std::get<0>(i) << ": ";
            for (size_t k = 0; k < std::get<1>(i).size(); ++k) {
                pers_fout << std::get<1>(i)[k];
                if (k != std::get<1>(i).size() - 1) {
                    pers_fout << " + ";
                }
            }
            pers_fout << std::endl;
        } 
        pers_fout << "-----------------------" << std::endl;    
    }

    pers_fout.close();

    // Now, produce a file with the map from the simplices to the unique_id that they are first assigned. For this, we can simply use the i_to_id map.
    std::ofstream i_to_id_fout(purename + "_id_to_simp");
    // Iterate over filt_simp and for i in i_to_id, add the simp -> id mapping to the file.
    for (size_t p = 0; p <= m; ++p) 
    {
        i_to_id_fout << "Dimension " << p << ": " << std::endl;
        for (auto i : i_to_id[p]) 
        {
            i_to_id_fout << i.second << " -> ";
            for (auto j : filt_simp[i.first]) 
            {
                i_to_id_fout << j << " ";
            }
            i_to_id_fout << std::endl;
        }
        i_to_id_fout << "-----------------------" << std::endl;  
    }
    i_to_id_fout.close();

    return 0;
}