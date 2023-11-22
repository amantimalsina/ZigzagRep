#include <cassert>
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


void run_persistence(
    const std::vector<std::vector<int> > &filt_simp, 
    const std::vector<bool> &filt_op,
    const int m,
    const std::string &infilename
)
{
        ZZREP::ZigzagRep zzr;
        // Let's measure the time it takes to compute the zigzag rep:
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > persistence;
        std::vector <std::vector<int> > id_to_i(m+1, std::vector<int>());
        zzr.compute(
            filt_simp, 
            filt_op,
            &persistence,
            &id_to_i,
            m);
        std::string purename;
        getFilePurename(infilename, &purename);
        std::ofstream pers_fout("../outputs/" + purename + "_pers_");
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        std::cout << "Runtime = " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
        << "[ms]" << std::endl;

        // Change this to add the representatives to the file.
        pers_fout << "n: " << filt_simp.size() << std::endl;
        pers_fout << "m: " << m << std::endl;
        pers_fout << "Runtime (in seconds): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000.0 << std::endl;

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
        std::ofstream i_to_id_fout("../outputs/" + purename + "_id_to_simp_");
        // Iterate over filt_simp and for i in i_to_id, add the simp -> id mapping to the file.
        for (size_t p = 0; p <= m; ++p) 
        {
            i_to_id_fout << "Dimension " << p << ": " << std::endl;
            for (size_t i = 0; i < id_to_i[p].size(); ++i) 
            {
                i_to_id_fout << i << " -> ";
                std::vector<int> simp_i = filt_simp[id_to_i[p][i]];
                for (auto j: simp_i) 
                {
                    i_to_id_fout << j << " ";
                }
                i_to_id_fout << std::endl;
            }
        }
        i_to_id_fout.close();
}

void get_runtime(
    const std::vector<std::vector<int> > &filt_simp, 
    const std::vector<bool> &filt_op,
    const int m,
    const std::string &infilename
)
{
        std::string purename;
        getFilePurename(infilename, &purename);
        std::ofstream runtime_fout("../outputs/" + purename + "_runtime_");

        for (size_t i = 0; i < filt_simp.size(); i += 50000) {
                ZZREP::ZigzagRep zzr;
                // Let's measure the time it takes to compute the zigzag rep:
                std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > persistence;
                std::vector <std::vector<int> > id_to_i(m+1, std::vector<int>());
                std::vector<std::vector<int>> filt_simp_trunc(filt_simp.begin(), filt_simp.begin() + i);
                std::vector<bool> filt_op_trunc(filt_op.begin(), filt_op.begin() + i);
                zzr.compute(
                    filt_simp_trunc, 
                    filt_op_trunc,
                    &persistence,
                    &id_to_i,
                    m);
                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                double runtime =  std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()/1000.0;
                runtime_fout << i << " " << runtime << std::endl;
        }
        runtime_fout.close();
}

void representative_test(
    const std::vector<std::vector<int> > &filt_simp, 
    const std::vector<bool> &filt_op,
    const std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > &persistence,
    const std::vector <std::vector<int> > &id_to_i
)
{
    /*
    Test for checking the correctness of representatives. Here are the criterias:
    1. If filtop[b - 1] = 1 (the map is injective), the simplex being inserted at b - 1 should be in the representative at b.
    2. If filtop[b - 1] = 0 (the map is surjective), the representative's class is the non-zero element in the kernel of this map.
    3. For i < n, if filtop[i] = 0 (injective), the simplex being inserted at i should be in the representative at i.
    4. For i < n, if filtop[i] = 1 (surjective), the representative's class is the non-zero element in the kernel of this map.
    5. Each representative at index b <= j <= i should be present in the complex K_j.
    6. The map \psi_j: H(K_{j-1}) \lefrightarrow H(K_{j}) takes the representative at index j-1 to the representative at index j.
    */
    size_t n = filt_op.size();
    for (auto pers: persistence) {
        int birth = std::get<0>(pers);
        int death = std::get<1>(pers);
        int dim = std::get<2>(pers);
        std::vector<std::tuple<int, std::vector<int>>> reps = std::get<3>(pers);
        if (filt_op[birth - 1]) {
            bool found = false;
            std::vector<int> birth_simp = filt_simp[birth-1];
            std::vector<int> reps_at_birth = std::get<1>(reps[0]);
            for (auto simp_id: reps_at_birth) {
                int first_i = id_to_i[dim][simp_id];
                std::vector<int> first_simp = filt_simp[first_i];
                // Check if first_simp is the same as death_simp.
                if (first_simp.size() == birth_simp.size()) {
                    bool same = true;
                    for (size_t i = 0; i < first_simp.size(); ++i) {
                        if (first_simp[i] != birth_simp[i]) {
                            same = false;
                            break;
                        }
                    }
                    if (same) {
                        found = true;
                        break;
                    }
                }
            }  
        }
        if (death < n) 
        {
            if (!filt_op[death])
            {
                bool found = false;
                std::vector<int> death_simp = filt_simp[death];
                // reps_at_death is the last element of reps.
                std::vector<int> reps_at_death = std::get<1>(reps[reps.size() - 1]);
                for (auto simp_id: reps_at_death) {
                    int first_i = id_to_i[dim][simp_id];
                    std::vector<int> first_simp = filt_simp[first_i];
                    // Check if first_simp is the same as death_simp.
                    if (first_simp.size() == death_simp.size()) {
                        bool same = true;
                        for (size_t i = 0; i < first_simp.size(); ++i) {
                            if (first_simp[i] != death_simp[i]) {
                                same = false;
                                break;
                            }
                        }
                        if (same) {
                            found = true;
                            break;
                        }
                    }   
                }
                if (!found) {
                    std::cout << "death: " << death << std::endl;
                    std::cout << "dim: " << dim << std::endl;
                    std::cout << "reps_at_death: ";
                    for (auto simp_id: reps_at_death) {
                        std::cout << simp_id << " ";
                    }
                    std::cout << std::endl;
                    std::cout << "id_to_i[dim]: ";
                    for (auto simp_id: id_to_i[dim]) {
                        std::cout << simp_id << " ";
                    }
                    std::cout << std::endl;
                }
            }
        }
    }
}


int main(const int argc, const char *argv[]) {
    if (argc < 3) 
    { std::cerr << "Err: input not large enough" << std::endl; return -1; }

    // Diplay the arguments:
    for (int i = 0; i < argc; ++i) {
        std::cout << "i" << i << ": " << argv[i] << std::endl;
    }


    if (atoi(argv[2]) == 1)
    {
        std::cout << "Getting runtimes." << std::endl;
    }
    else {
        std::cout << "Running persistence." << std::endl;
    }

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

    // Run the compute algorithm by truncating the input in the increment of 50,000 simplices: 
    if (atoi(argv[2]) == 1) {
        get_runtime(filt_simp, filt_op, m,infilename);
    }
    else {
        run_persistence(filt_simp, filt_op, m, infilename);
    }

    return 0;
}