#include "zzrep.h"

#include <algorithm>
#include <unordered_map>
#include <tuple>
#include <boost/container_hash/hash.hpp> // Changed functional to container_hash

using namespace std;

namespace ZZREP { 
 
template <class ElemType>
class VecHash { 
public:
    size_t operator()(const vector<ElemType>& v) const; 
};

template <class ElemType>
size_t VecHash<ElemType>
    ::operator()(const vector<ElemType>& v) const {

    size_t seed = 0;

    for (auto e : v) { boost::hash_combine(seed, e); }

    return seed;
}

template <class ElemType>
class VecEqual { 
public:
    bool operator()(const vector<ElemType>& v1, 
        const vector<ElemType>& v2) const; 
};

template <class ElemType>
bool VecEqual<ElemType>
    ::operator()(const vector<ElemType>& v1, 
        const vector<ElemType>& v2) const {

    if (v1.size() != v2.size()) { return false; }

    for (auto i = 0; i < v1.size(); i ++) {
        if (v1[i] != v2[i]) {
            return false;
        }
    }

    return true;
}

typedef unordered_map< vector<int>, int,
    VecHash<int>, VecEqual<int> > SimplexIdMap;
/*
        Description: 
            Adds two binary vectors using bitwise XOR
        Input: 
             binary vectors a,b.
        Output: 
            pointer to c = a+b.
*/

typedef vector<int> column;


void add_columns(std::vector<int> *a, std::vector<int> *b, std::vector<int> *c);

int pivot (std::vector<int> *a);

/*
        Description: 
            Explores if there is a pivot conflict among the columns in a matrix . 
        Input: 
            A vector of vectors of ints representing a 2D matrix.
        Output: 
            (true, i, j) if there is a pivot conflict between columns i and j,
            (false, -1, -1) otherwise.
*/
tuple<int,  int, bool> pivot_conflict (std::vector<std::vector<int> > *matrix);

int pivot (std::vector<int> *a);

void reduce(vector<int> *a, vector<vector<int> > *M, vector<int> *indices);

void find_indices(vector<vector<int>> *I, vector<int> *indices);

void ZigzagRep::compute(
        const std::vector<vector<int> > &filt_simp, 
        const std::vector<bool> &filt_op,
        std::vector <tuple <int, int, int, vector<int> > > *persistence, int m) {
    
    persistence -> clear();
    int n = filt_simp.size();
    /* 
    Define a cycle matrix Z[p] and a chain matrix C[p]. 
    1. The number of columns of Z[p] and C[p-1] equals rank Z_p(K_i) and the number of rows of Z[p] and C[p-1] equals n.
    2. Each column of Z[p] and C[p] represents a p-chain in K_i such that for each simplex s in K_i, s belongs to the p-chain if and only if the bit with index id[s] in the column equals 1.
    */
    vector<vector<vector<int> > > Z;
    vector<vector<vector<int> > > C;
    /*
    Define a unique integral id less than n to each simplex.
    Moreover, for each column Z_p[j], a birth timestamp b_p[j] is maintained.
    */
    vector<SimplexIdMap> id;
    vector<vector<int> > birth_timestamp;

    for (int p = 0; p <= m; p++) {
        vector<column> Z_p;
        vector<column> C_p;
        SimplexIdMap id_p;
        column birth_timestamp_p;
        Z.push_back(Z_p);
        C.push_back(C_p);
        id.push_back(id_p);
        birth_timestamp.push_back(birth_timestamp_p);
    }
    for (int i = 0; i < n; ++i) {
        const vector<int> &simp = filt_simp[i];
        int p = simp.size() - 1; // p denotes the dimension of the simplex.
        
        if (filt_op[i]) {
            // INSERTION:
            // A p-simplex is inserted into id[p] and a (p-1)-simplex is inserted into pm1_id. Next, we add a row to Z[p] and C[p], according to the dimension of simp.
            int j = id[p].size();
            (id[p]).emplace(simp, j);
            // Add a row with all zeros at the end to Z[p] and C[p].
            int k = C[p].size();
            if (k != 0) {
                column last_column_C;
                // Add all zeros to this last column:
                for (int j = 0; j < k; j++) {
                    last_column_C.push_back(0);
                }
                C[p].push_back(last_column_C);
                for (int j = 0; j < C[p].size(); j++) {
                    if (j == C[p].size() - 1) {
                        C[p][j].push_back(1);
                    }
                    else {
                        Z[p][j].push_back(0);
                        C[p][j].push_back(0);
                    }
                }
            }
            else {
                C[p].push_back({1});
                if (p == 0)  Z[p].push_back({1});
            }
            bool all_boundary = true;
            // Compute boundary of the simplex:
            column bd_simp;
            // Represent the boundary of simp as a sum of columns of Z_{p-1} by a reduction algorithm; I is such set of columns.
            column I;
            if (p != 0) {
                int b = C[p-1].size();
                for (int j = 0; j < b; j++) {
                    bd_simp.push_back(0);
                }
                for (int i = 0; i <= p; i++) {
                    // Remove the i-th vertex.
                    vector<int> boundary_simplex = simp;
                    boundary_simplex.erase(boundary_simplex.begin() + i);
                    bd_simp[id[p-1][boundary_simplex]] = 1;
                }
                reduce(&bd_simp, &Z[p-1], &I);
                // Check the birth timestamp of a to check whether all of them are boundaries.
                for (auto a: I) {
                    if (birth_timestamp[p][a] >= 0) {
                        all_boundary = false;
                        break;
                    }
                }
            }
            if (all_boundary) // simp is a boundary
            {
                /*
                An interval in dimension p gets born: 
                Append a new column simp + \sum_{a \in I} C^{p}[a] with birth timestamp i+1 to Z^p.
                */
                vector<int> new_column;
                if (p != 0)  {
                    for (int j = 0; j < p; j++) {
                        bd_simp.push_back(0);
                    }
                    new_column[id[p][simp]] = 1;
                    for (auto a: I) {
                        add_columns(&new_column, &C[p][a], &new_column);
                    }
                    Z[p].push_back(new_column);
                }
                else if (p == 0 && k != 0) {
                    int k = C[p].back().size();
                    for (int j = 0; j < k; j++) {
                        new_column.push_back(0);
                    }
                    new_column[id[p][simp]] = 1;
                    Z[p].push_back(new_column);
                }
                birth_timestamp[p].push_back(i+1);
            }
            else {
                // An interval in dimension p-1 dies.
                /*
                Let J consist of indices in I whose corresponding columns in Z_{p−1} have non-negative birth timestamps. 
                */ 
                vector<int> J;
                for (auto a: I) {
                    if (birth_timestamp[p-1][a] >= 0) { // Gather all the non-boundary cycles.
                        J.push_back(a);
                    }
                }
                /*
                If φ_{b^{p−1}[c]−1} (the arrow at b^{p−1}[c]−1) points backward for all c in J, let k be the smallest index in J; 
                otherwise, let k be the largest c in J such that φ_{b^{p−1}[c]−1} points forward. 
                */
                int k;
                bool arrow_backward = true;
                // Check if arrow at b^{p−1}[c]−1 points backward for all c in J
                for (auto c: J) {
                    if (filt_op[birth_timestamp[p-1][c]-1]) {
                        arrow_backward = false;
                        break;
                    }
                }
                if (arrow_backward) {
                    k = *J.begin();
                }
                else {
                    k = *J.rbegin();
                }
                // Output the (p − 1)-th interval [b^{p−1}[k], i].
                persistence->push_back(std::make_tuple(birth_timestamp[p-1][k], i, p-1, Z[p-1][k]));
                // Set Z_{p−1}[k] = bd.(simp), C[p][k] = simp, and b^{p−1}[k] = −1.
                vector<int> bd_simp;
                for (int j = 0; j < p; j++) {
                    bd_simp.push_back(0);
                }
                for (int i = 0; i <= p; i++) {
                    // Remove the i-th vertex.
                    vector<int> boundary_simplex = simp;
                    boundary_simplex.erase(boundary_simplex.begin() + i);
                    bd_simp[id[p-1][boundary_simplex]] = 1;
                }
                Z[p-1][k] = bd_simp;
                C[p][k] = simp;
                birth_timestamp[p-1][k] = -1;
                /*
                Avoiding pivot conflicts (column of bd.(simp) with that of another column in Z_{p-1}) as follows:
                    
                */
                tuple<int, int, bool> pivot_conflict_tuple =  pivot_conflict(&Z[p-1]);
                int a = get<0>(pivot_conflict_tuple);
                int b = get<1>(pivot_conflict_tuple);
                bool exists_pivot_conflict = get<2>(pivot_conflict_tuple);
                while (exists_pivot_conflict) {
                    vector<int> Z_pm1_aplusb;
                    add_columns(&Z[p-1][a], &Z[p-1][b], &Z_pm1_aplusb);
                    if (birth_timestamp[p-1][a] < 0 && birth_timestamp[p-1][b] < 0) {
                        Z[p-1][a] = Z_pm1_aplusb;
                        add_columns(&C[p][a], &C[p][b], &C[p][a]);
                    }
                    else if (birth_timestamp[p-1][a] < 0 && birth_timestamp[p-1][b] >= 0) {
                        Z[p-1][b] = Z_pm1_aplusb;
                    }
                    else if (birth_timestamp[p-1][a] >= 0 && birth_timestamp[p-1][b] < 0) {
                        Z[p-1][a] = Z_pm1_aplusb;
                    }
                    else {
                        bool a_lessthan_b = ((a == b) || ((a < b) && (filt_op[b-1])) || ((a > b) && (!(filt_op[a-1]))));   
                        if (a_lessthan_b) {
                            Z[p-1][b] = Z_pm1_aplusb;
                        }
                        else {

                            Z[p-1][a] = Z_pm1_aplusb;
                        }
                    }
                    pivot_conflict_tuple =  pivot_conflict(&Z[p-1]);
                    a = get<0>(pivot_conflict_tuple);
                    b = get<1>(pivot_conflict_tuple);
                    exists_pivot_conflict = get<2>(pivot_conflict_tuple);
                    }
            }
        } else { // Case: Arrow points backward.
            bool existence = false;
            for (auto column: Z[p]) {
                if (column[id[p][simp]]==1) {
                    existence = true;
                    break;
                }
            } 
            if (existence)// Case: Birth.
            {   
                int a, b;
                bool boundary_pairs_bool;
                for (int i = 0; i < Z[p-1].size(); i++) 
                {
                    for (int j = 0; j < Z[p-1].size(); j++) 
                    {
                        if (birth_timestamp[p-1][i] < 0 && birth_timestamp[p-1][j] < 0 && (C[p][i][id[p][simp]]==1) && (C[p][j][id[p][simp]]==1)) 
                        {
                            a = i; b = j; boundary_pairs_bool = true;
                        }
                    }
                }
                while (boundary_pairs_bool) {
                    vector<int> Z_pm1_aplusb;
                    vector<int> C_p_aplusb;
                    add_columns(&Z[p-1][a], &Z[p-1][b], &Z_pm1_aplusb);
                    add_columns(&C[p][a], &C[p][b], &C_p_aplusb);
                    if (pivot(&Z[p-1][a]) > pivot(&Z[p-1][b])) 
                    {
                        Z[p-1][a] = Z_pm1_aplusb;
                        C[p][a] = C_p_aplusb;
                    }
                    else 
                    {
                        Z[p-1][b] = Z_pm1_aplusb;
                        C[p][b] = C_p_aplusb;
                    }
                    for (int i = 0; i < Z[p-1].size(); i++) 
                    {
                        for (int j = 0; j < Z[p-1].size(); j++) 
                        {
                            if (birth_timestamp[p-1][i] < 0 && birth_timestamp[p-1][j] < 0 && (C[p][i][id[p][simp]]==1) && (C[p][j][id[p][simp]]==1)) 
                            {
                                a = i; b = j; boundary_pairs_bool = true;
                            }
                        }
                    }
                }
                // find the only column Z_{p-1}[a] with negative birth timestamp such that C[p][a] contains the simplex.
                vector<int> column;
                int x = 0;
                for (x = 0; x < Z[p-1].size(); x++) {
                    if ((birth_timestamp[p-1][x] < 0) && (C[p][x][id[p][simp]]==1)) {
                        column = Z[p-1][x];
                        break;
                    }
                }
                birth_timestamp[p-1][x] = i+1;
            }
            else // Case: Death.
            {
                // Update C[p] so that no columns contain the simplex.
                for (int a = 0; a < Z[p-1].size(); a++)
                {
                    if (Z[p-1][a][id[p][simp]]==1) {
                        for (int b = 0; b < C[p].size(); b++) 
                        {
                            if ((C[p][b][id[p][simp]]==1) && (birth_timestamp[p-1][b] < 0)) // each column in C[p][b] containing the simplex for which Z_{p-1}[b] is a boundary).
                            {    
                                add_columns(&Z[p][a], &C[p][b], &C[p][b]);
                                add_columns(&C[p][b], &Z[p][a], &C[p][b]);
                            }
                        }
                    }
                }
                // Remove the simplex from Z[p]
                // I contains the columns that contain the simplex.
                vector<int> I;
                for (int a = 0; a < Z[p-1].size(); a++) // each column Z[p][a] containing the simplex
                {
                    if (Z[p-1][a][id[p][simp]]==1) I.push_back(a);
                }
                // sort I in the order of the birth timestamps where the order is the total order as above.
                sort(I.begin(), I.end(), [&](int &a, int &b){ return ((I[a] == b) || ((I[a] < b) && (filt_op[I[b]-1])) || ((I[a] > I[b]) && (!(filt_op[I[a]-1]))));});
                vector<int> z = Z[p][I[0]];
                for (auto a: I)
                {
                    if (pivot(&Z[p][a]) > pivot(&z)) {
                        add_columns(&Z[p][a], &z, &Z[p][a]);
                    }   
                    else 
                    {
                        vector<int> temp = Z[p][a];
                        add_columns(&Z[p][a], &z, &Z[p][a]);
                        z = temp;
                    }
                }
                // Output the p-th interval [birth_timestamp_p[I[0]], i].
                persistence->push_back(std::make_tuple(birth_timestamp[p][I[0]], i, p, Z[p][I[0]]));
                // Delete the column Z[p][I[0]] from Z[p] and delete birth_timestamp_p[I[0]] from birth_timestamp_p.
                Z[p].erase(Z[p].begin() + I[0]);
                birth_timestamp[p].erase(birth_timestamp[p].begin() + I[0]);
            }
        }
    }
    // For each p and each column Z[p][a] of Z[p] with non-negative birth timestamp, output the p-th interval [birth_timestamp_p[a], n].
    for (int p = 0; p <= m; p++) {
        for (int a = 0; a < Z[p].size(); a++) {
            if (birth_timestamp[p][a] >= 0) {
                persistence->push_back(std::make_tuple(birth_timestamp[p][a], n, p, Z[p][a]));
            }
        }
    }
}

void add_columns(vector<int> *a, vector<int> *b, vector<int> *c)
{
    int length = a->size();
    // Iterate over a and b and add them bitwise to c.
    for (int i = 0; i < length; i++) {
        (*c)[i] = (*a)[i] ^ (*b)[i];
    }
}

// Goes over a vector of ints and finds the position of the last non-zero element.
int pivot (vector<int> *a)
{
    int length = a->size();
    for (int i = length - 1; i >= 0; i--) {
        if ((*a)[i] != 0) {
            return i;
        }
    }
    return -1;
}

// Goes over the columns of the matrix and returns the first pair with the same pivots.
std::tuple<int,  int, bool> pivot_conflict (vector<vector<int> > *matrix)
{
    int prev_pivot = -1;
    int prev_pivot_column = -1;
    int current_pivot = -1;
    int current_pivot_column = -1;
    bool exists_pivot_conflict = false;
    for (int i = 0; i < matrix->size(); i++) {
        current_pivot = pivot(&(*matrix)[i]);
        current_pivot_column = i;
        if (current_pivot == prev_pivot) {
            exists_pivot_conflict = true;
            break;
        }
        prev_pivot = current_pivot;
        prev_pivot_column = current_pivot_column;
    }
    return std::make_tuple(prev_pivot_column, current_pivot_column, exists_pivot_conflict);
}


// Reduction algorithm: given a column a, find the indices of the columns in M such that the columns sum to a.
void reduce(vector<int> *a, vector<vector<int>> *M, vector<int> *indices)
{
    // Append a to M.
    M -> push_back(*a);
    int k = M -> size();
    std::vector<vector<int>> I;
    std::vector<int> pivot_entries;
    for (int i = 0; i < M->size(); i++) {
        pivot_entries.push_back(pivot(&((*M)[i])));
        std::vector<int> I_i;
        I.push_back(I_i);
        bool pivot_conflict = false;
        int conflicting_column;
        do
        {
            // Check for pivot conflict
            for (int j = 0; j < i; j++) {
                if (pivot_entries[i] == pivot_entries[j]) {
                    pivot_conflict = true;
                    conflicting_column = j;
                    break;
                }
                else {
                    pivot_conflict = false;
                }
            }
            // Resolve pivot conflict and keep track of the column that was added in I.
            if (pivot_conflict) {
                add_columns(&((*M)[i]), &((*M)[conflicting_column]), &((*M)[i]));
                pivot_entries[i] = pivot(&((*M)[i]));
                I[i].push_back(conflicting_column);
            }
        }
        while (pivot_conflict);
    }
    // Find all the indices that add to the zeroing of the boundary.
    find_indices(&I, indices);
}

void find_indices(vector<vector<int>> *I, vector<int> *indices)
{
    std::unordered_map<int, bool> visited;
    // Add all the indices that got added to the last column.
    indices -> push_back(I -> size() - 1);
    for (auto i = (indices -> begin()); i != indices -> end(); i++) {
        // Add all the indices that got added to this column.
        for (int j = 0; j < I[*i].size(); j++) {
            if (!visited[(*I)[*i][j]]) {
                indices -> push_back((*I)[*i][j]);
                visited[(*I)[*i][j]] = true;
            }
        }
    }
    // Remove the last column.
    indices -> pop_back();
}

}// namespace ZZREP