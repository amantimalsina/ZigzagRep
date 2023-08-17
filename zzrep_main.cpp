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
            assert(op == 'd');
            filt_op.push_back(false);
        }
    }


    filt_fin.close();

    std::vector< std::tuple<int, int, int, std::vector<std::vector<int> >>> persistence;
    ZZREP::ZigzagRep zzr;
    zzr.compute(
        filt_simp, 
        filt_op,
        &persistence,
        m);
    std::string purename;
    getFilePurename(infilename, &purename);
    std::ofstream pers_fout(purename + "_pers");

    // Change this to also display representatives.
    for (const auto& e : persistence) {
        pers_fout << std::get<2>(e) << " " << std::get<0>(e) 
            << " " << std::get<1>(e) << std::endl;   
        pers_fout << "rep: " << std::endl;
        for (auto simp : std::get<3>(e)) {
            for (auto i : simp) {pers_fout << i << " ";}
            pers_fout << std::endl;
        } 
        pers_fout << std::endl << "-----------------------" << std::endl;    
    }

    return 0;
}