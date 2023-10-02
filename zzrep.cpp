#include "zzrep.h"

#include <algorithm>
#include <unordered_map>
#include <tuple>
#include <boost/bimap.hpp>
#include <boost/dynamic_bitset.hpp>

using namespace std;

namespace ZZREP { 

typedef boost::bimaps::bimap<vector<int>, int> SimplexIdMap;
typedef boost::bimaps::bimap<int, int> CycleToChainMap;
typedef std::unordered_map<int, int> RepToCycleMap;
typedef SimplexIdMap::value_type SimplexIdPair;
typedef CycleToChainMap::value_type CycleToChainPair;
typedef RepToCycleMap::value_type RepToCyclePair;

typedef boost::dynamic_bitset<> column;

void add_columns(column *a, column *b, column *c);

int pivot (column *a);

/*
        Description: 
            Explores if there is a pivot conflict among the columns in a matrix . 
        Input: 
            A vector of vectors of ints representing a 2D matrix.
        Output: 
            (true, i, j) if there is a pivot conflict between columns i and j,
            (false, -1, -1) otherwise.
*/
tuple<int,  int> pivot_conflict (std::vector<column> *matrix);

int pivot (column *a);

void reduce(column *a, vector<column> *M, vector<int> *indices);

int find_last(column *a);

int index_binary_search(vector<int> &a, int x);

void ZigzagRep::compute(
        const std::vector<vector<int> > &filt_simp, 
        const std::vector<bool> &filt_op,
        std::vector <tuple <int, int, int, std::vector<std::vector<int> > > > *persistence, int m) {
    
    persistence -> clear();
    int n = filt_simp.size();
    /* 
    Define a cycle matrix Z[p] and a chain matrix C[p]. 
    1. The number of columns of Z[p] and C[p-1] equals rank Z_p(K_i) and the number of rows of Z[p] and C[p-1] equals n.
    2. Each column of Z[p] and C[p] represents a p-chain in K_i such that for each simplex s in K_i, s belongs to the p-chain if and only if the bit with index id[s] in the column equals 1.
    */
    vector<vector<column> > Z;
    vector<vector<column> > C;
    /*
    Define a unique integral id less than n to each simplex.
    Moreover, for each column Z_p[j], a birth timestamp b_p[j] is maintained.
    */
    vector<SimplexIdMap> id;
    vector<vector<int>> birth_timestamp;
    vector<CycleToChainMap> cycle_to_chain;
    bool next_op = filt_op[1];

    /*
    Representative Data Structures:
    */
    vector<vector<vector<column>>> representatives; // representatives[p][i][j] is a p-cycle representative corresponding to the interval i at point j.
    vector<vector<vector<int>>> index_rep; // a sorted list for each p and the filtration index j corresponding to index of representatives[p][i][j] for i in Z[p].
    vector<RepToCycleMap> rep_to_cycle; // a map from the index of representatives[p][i][j] to the index of the cycle in Z[p] that it represents.

    for (int p = 0; p <= m; p++) {
        vector<column> Z_p;
        vector<column> C_p;
        SimplexIdMap id_p;
        CycleToChainMap cycle_to_chain_p;
        vector<int> birth_timestamp_p;
        Z.push_back(Z_p);
        C.push_back(C_p);
        id.push_back(id_p);
        cycle_to_chain.push_back(cycle_to_chain_p);
        birth_timestamp.push_back(birth_timestamp_p);
    }
    for (int i = 0; i < n; ++i) {
        const vector<int> &simp = filt_simp[i];
        int p = simp.size() - 1; // p denotes the dimension of the simplex.
        if (filt_op[i]) {
            // INSERTION:
            // A p-simplex is inserted into id[p] and a (p-1)-simplex is inserted into pm1_id. Next, we add a row to Z[p] and C[p], according to the dimension of simp.
            int k = id[p].size();
            (id[p]).insert(SimplexIdPair(simp, k));
            // Add a row with all zeros at the end to Z[p] and C[p].
            if (k != 0) {
                for (int j = 0; j < C[p].size(); j++) {
                    C[p][j].push_back(0);
                }
                for (int j = 0; j < Z[p].size(); ++j) {
                    Z[p][j].push_back(0);
                }
            }
            else {
                if (p == 0) {
                    Z[p].push_back(column(1,1));
                }
            }
            bool all_boundary = true;
            // Compute boundary of the simplex:
            column bd_simp;
            // Represent the boundary of simp as a sum of columns of Z_{p-1} by a reduction algorithm; I is such set of columns.
            vector<int> I;
            if (p != 0) {
                column temp(id[p-1].size(), 0);
                bd_simp = temp;
                for (int i = 0; i <= p; i++) {
                    // Remove the i-th vertex.
                    vector<int> boundary_simplex = simp;
                    boundary_simplex.erase(boundary_simplex.begin() + i);
                    int idx = id[p-1].left.at(boundary_simplex);
                    bd_simp[idx] = 1;
                }
                reduce(&bd_simp, &Z[p-1], &I);
                // Check the birth timestamp of a to check whether all of them are boundaries.
                for (auto a: I) {
                    if (birth_timestamp[p-1][a] >= 0) {
                        all_boundary = false;
                        break;
                    }
                }
            }
            if (all_boundary) // simp completes a cycle
            {
                /*
                An interval in dimension p gets born: 
                Append a new column simp + \sum_{a \in I} C^{p}[a] with birth timestamp i+1 to Z^p.
                */
                column new_column = column(id[p].size(), 0); // New column with of size id[p].size() with 1 appended at the end.
                new_column[id[p].left.at(simp)] = 1;
                if (p != 0)  {
                    for (auto a: I) {
                        int chain_a = cycle_to_chain[p-1].left.at(a);
                        new_column ^= C[p][chain_a];
                    }
                    Z[p].push_back(new_column);
                }
                else if (p == 0 && k != 0) {
                    Z[p].push_back(new_column);
                }
                birth_timestamp[p].push_back(i+1);
                /*
                // Record this cycle for storing representatives at this stage. This is getting born and we simply store the cycle at this point.
                representatives[p][i].push_back(new_column);
                index_rep[p][i].push_back(k);
                rep_to_cycle[p].insert(RepToCyclePair(k, k));
                */
            }
            else {
                // An interval in dimension p-1 dies.
                /*
                Let J consist of indices in I whose corresponding columns in Z_{p−1} have non-negative birth timestamps. 
                */ 
                vector<int> J;
                for (int a: I) {
                    if (birth_timestamp[p-1][a] >= 0) { // Gather all the non-boundary cycles.
                        J.push_back(a);
                    }
                }
                // Check if arrow at b^{p−1}[c]−1 points backward for all c in J
                int l;
                bool arrow_backward = true;
                for (auto c: J) {
                    if (filt_op[birth_timestamp[p-1][c]-1]) {
                        arrow_backward = false;
                        break;
                    }
                }
                sort(J.begin(), J.end());
                if (arrow_backward) {
                    // the arrow at b^{p−1}[c]−1) points backward for all c in J, let k be the smallest index in J.
                    l = *J.begin();
                }
                else {
                    // let k be the largest c in J such that φ_{b^{p−1}[c]−1} points forward.
                    for (auto i = J.size()-1; i >= 0; i--) {
                        if (filt_op[birth_timestamp[p-1][J[i]]-1]) {
                            l = J[i];
                            break;
                        }
                    }
                }
                // Output the (p − 1)-th interval [b^{p−1}[l], i].
                vector<vector<int>> representative;
                for (int i = 0; i < Z[p-1][l].size(); i++) {
                    if (Z[p-1][l][i] == 1) {
                        representative.push_back(id[p-1].right.at(i));
                    }
                }
                /* Representative updates for the forward case: 
                // Iterate over J and update the representatives.
                for (auto j: J) {
                    if (j != l) {
                        // Get the points where the representatives are stored for j:
                        vector<int> index_rep_j = index_rep[p-1][birth_timestamp[p-1][j]-1];  
                        // Iterate over the representatives and update them:
                        for (auto r: index_rep_j) {
                            // Find the smallest index closes to r in index_rep[p-1][l][r] and get the reprentative at that index:
                            int index_rep_l = index_binary_search(index_rep[p-1][j], r);
                            column representative_l = representatives[p-1][birth_timestamp[p-1][l]-1][index_rep_l];
                            column representative_j = representatives[p-1][birth_timestamp[p-1][j]-1][r];
                            // Update the representative by adding a column pointing to r:
                            representatives[p-1][birth_timestamp[p-1][l]-1].push_back(representative_l ^ representative_j);
                            index_rep[p-1][birth_timestamp[p-1][l]-1].push_back(r);
                            // Update the map from representatives to cycles:
                            rep_to_cycle[p-1].insert(RepToCyclePair(r, l));
                        }
                    }
                }
                */
                persistence->push_back(std::make_tuple(birth_timestamp[p-1][l], i, p-1, representative));
                // Set Z_{p−1}[l] = bd.(simp), C[p][l] = simp, and b^{p−1}[l] = −1.
                Z[p-1][l] = bd_simp;
                column current_chain(id[p].size(), 0);
                current_chain[id[p].left.at(simp)] = 1;
                C[p].push_back(current_chain);
                (cycle_to_chain[p-1]).insert(CycleToChainPair(l, C[p].size()-1));
                birth_timestamp[p-1][l] = -1;
                /*
                Avoiding pivot conflicts (column of bd.(simp) with that of another column in Z_{p-1}) as follows:
                FIXME: Check if this is correct or not.
                */
                tuple<int, int> pivot_conflict_tuple =  pivot_conflict(&Z[p-1]);
                int a = get<0>(pivot_conflict_tuple);
                int b = get<1>(pivot_conflict_tuple);
                bool exists_pivot_conflict = (a != b);
                while (exists_pivot_conflict) {
                    column Z_pm1_aplusb;
                    //add_columns(&Z[p-1][a], &Z[p-1][b], &Z_pm1_aplusb);
                    Z_pm1_aplusb = Z[p-1][a] ^ Z[p-1][b];
                    if (birth_timestamp[p-1][a] < 0 && birth_timestamp[p-1][b] < 0) {
                        Z[p-1][a] = Z_pm1_aplusb;
                        // Update chain matrices to reflect this addition (use the cycletochainmap):
                        int chain_a = cycle_to_chain[p-1].left.at(a);
                        int chain_b = cycle_to_chain[p-1].left.at(b);
                        C[p][chain_a] = C[p][chain_a] ^ C[p][chain_b];
                        // FIXME: Check if we need to update the cycle to chain map:
                        // cycle_to_chain[p-1].erase(b);
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
                    exists_pivot_conflict = (a != b);
                    }
            }
        } 
        else { // Case: Arrow points backward.
            bool existence = false;
            int idx = id[p].left.at(simp);
            for (auto column: Z[p]) {
                if (column[idx]==1) {
                    existence = true;
                    break;
                }
            } 
            if (!existence)// Case: Birth.
            {   
                int a, b;
                bool boundary_pairs_bool = false;
                for (int x = 0; x < Z[p-1].size(); x++) 
                {
                    for (int y = x + 1; y < Z[p-1].size(); y++) 
                    {
                        if (birth_timestamp[p-1][x] < 0)
                        {
                            if (birth_timestamp[p-1][y] < 0)
                            {
                                int chain_x = (cycle_to_chain[p-1]).left.at(x);
                                int chain_y = (cycle_to_chain[p-1]).left.at(y);
                                if (C[p][chain_x][idx]==1)
                                {
                                    if (C[p][chain_y][idx]==1)
                                    {
                                        a = x; b = y; boundary_pairs_bool = true;
                                    }
                                }
                            }
                        }
                    }
                }
                while (boundary_pairs_bool) {
                    column Z_pm1_aplusb;
                    column C_p_aplusb;
                    Z_pm1_aplusb = Z[p-1][a] ^ Z[p-1][b];
                    int chain_a = (cycle_to_chain[p-1]).left.at(a);
                    int chain_b = (cycle_to_chain[p-1]).left.at(b);
                    C_p_aplusb = C[p][chain_a] ^ C[p][chain_b];
                    if (pivot(&Z[p-1][a]) > pivot(&Z[p-1][b])) 
                    {
                        Z[p-1][a] = Z_pm1_aplusb;
                        C[p][chain_a] = C_p_aplusb;
                    }
                    else 
                    {
                        Z[p-1][b] = Z_pm1_aplusb;
                        C[p][chain_b] = C_p_aplusb;
                    }
                    boundary_pairs_bool = false;
                    for (int x = 0; x < Z[p-1].size(); x++) 
                    {
                        for (int y = x + 1; y < Z[p-1].size(); y++) 
                        {
                            if (birth_timestamp[p-1][x] < 0)
                            {
                                if (birth_timestamp[p-1][y] < 0)
                                {
                                    int chain_x = (cycle_to_chain[p-1]).left.at(x);
                                    int chain_y = (cycle_to_chain[p-1]).left.at(y);
                                    if (C[p][chain_x][idx]==1)
                                    {
                                        if (C[p][chain_y][idx]==1)
                                        {
                                            a = x; b = y; boundary_pairs_bool = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                // find the only column Z_{p-1}[a] with negative birth timestamp such that C[p][a] contains the simplex.
                int only_idx = 0;
                for (only_idx = 0; only_idx < Z[p-1].size(); only_idx++) {
                    int chain_x = cycle_to_chain[p-1].left.at(only_idx);
                    if ((birth_timestamp[p-1][only_idx] < 0) && (C[p][chain_x][only_idx]==1)) {
                        // Update the pth chain matrix so that the rows do not contain the simplex.
                        C[p][chain_x][only_idx] = 0;
                        birth_timestamp[p-1][only_idx] = i+1;
                        // Remove the cycle to chain record.
                        cycle_to_chain[p-1].left.erase(only_idx);
                        break;
                    }
                }
                // Record this cycle for storing representatives at this stage. This is getting born and we simply store the cycle at this point.
                /* 
                representatives[p-1][i].push_back(Z[p-1][only_idx]);
                index_rep[p-1][i].push_back(only_idx);
                rep_to_cycle[p-1].insert(RepToCyclePair(only_idx, only_idx));
                */
            }
            else // Case: Death.
            {
                // Update C[p] so that no columns contain the simplex.
                for (int a = 0; a < Z[p].size(); a++)
                {
                    if (Z[p][a][idx]==1) {
                        for (int b = 0; b < C[p].size(); b++) 
                        {
                            int cycle_b = cycle_to_chain[p].right.at(b);
                            if (birth_timestamp[p-1][cycle_b] < 0){
                                if (C[p][b][idx]==1) // each column in C[p][chain_b] containing the simplex s.t. Z[p-1][b] is a boundary.
                                {    
                                    //add_columns(&C[p][b], &Z[p][a], &C[p][b]);
                                    C[p][b] ^= Z[p][a];
                                }
                            }
                        }
                    }
                }
                // Remove the simplex from Z[p]
                // I contains the columns that contain the simplex.
                vector<int> I;
                for (int a = 0; a < Z[p].size(); a++) // each column Z[p][a] containing the simplex
                {
                    if (Z[p][a][idx]==1) I.push_back(a);
                }
                // sort I in the order of the birth timestamps where the order is the total order as above.
                sort(I.begin(), I.end(), [&](int &a, int &b){ return ((I[a] == b) || ((I[a] < b) && (filt_op[I[b]-1])) || ((I[a] > I[b]) && (!(filt_op[I[a]-1]))));});
                int alpha = I[0];
                I.erase(I.begin());
                column z = Z[p][alpha];
                for (auto a: I)
                {
                    // FIXME: Don't we need to update the cycle to chain map?
                    if (pivot(&Z[p][a]) > pivot(&z)) {
                        // add_columns(&Z[p][a], &z, &Z[p][a]);
                        Z[p][a] = Z[p][a] ^ z;
                    }   
                    else 
                    {
                        column temp = Z[p][a];
                        // add_columns(&Z[p][a], &z, &Z[p][a]);
                        Z[p][a] = Z[p][a] ^ z;
                        z = temp;
                    }
                }
                // Output the p-th interval [birth_timestamp_p[I[0]], i].
                vector<vector<int>> representative;
                for (int i = 0; i < Z[p][alpha].size(); i++) {
                    if (Z[p][alpha][i] == 1) {
                        representative.push_back(id[p].right.at(i));
                    }
                }
                persistence->push_back(std::make_tuple(birth_timestamp[p][alpha], i, p, representative));
                /* Representative updates for the backward case:
                // Iterate over I and update the representatives.
                for (auto a: I) {
                    if (a != alpha) {
                        // Get the points where the representatives are stored for j:
                        vector<int> index_rep_alpha = index_rep[p][birth_timestamp[p][alpha]-1];
                        // Iterate over the representatives and update them:
                        for (auto r: index_rep_alpha) {
                            // Find the smallest index closes to r in index_rep[p-1][l][r] and get the reprentative at that index:
                            int index_rep_a = index_binary_search(index_rep[p][a], r);
                            column representative_a = representatives[p-1][birth_timestamp[p-1][a]-1][index_rep_a];
                            column representative_alpha = representatives[p-1][birth_timestamp[p][alpha]-1][r];
                            // Update the representative by adding a column pointing to r:
                            representatives[p-1][birth_timestamp[p-1][a]-1].push_back(representative_a ^ representative_alpha);
                            index_rep[p-1][birth_timestamp[p-1][a]-1].push_back(r);
                            rep_to_cycle[p-1].insert(RepToCyclePair(r, a));
                        }
                    }
                }
                */
                // Delete the column Z[p][I[0]] from Z[p] and delete birth_timestamp_p[I[0]] from birth_timestamp_p.
                Z[p].erase(Z[p].begin() + alpha);
                birth_timestamp[p].erase(birth_timestamp[p].begin() + alpha);
            }
        }
    }
    // For each p and each column Z[p][a] of Z[p] with non-negative birth timestamp, output the p-th interval [birth_timestamp_p[a], n].
    for (int p = 0; p <= m; p++) {
        for (int a = 0; a < Z[p].size(); a++) {
            if (birth_timestamp[p][a] >= 0) {
                vector<vector<int>> representative;
                for (int i = 0; i < Z[p][a].size(); i++) {
                    if (Z[p][a][i] == 1) {
                        representative.push_back(id[p].right.at(i));
                    }
                }
                persistence->push_back(std::make_tuple(birth_timestamp[p][a], n, p, representative));
            }
        }
    }
}

void add_columns(column *a, column *b, column *c)
{
    *c = (*a) ^ (*b);
}

// Goes over a vector of ints and finds the position of the last non-zero element.
int pivot (column *a)
{
    int pivot = -1;
    for (int i = 0; i < a->size(); ++i) {
        if ((*a)[i] == 1) {
            pivot = i;
        }
    }
    return pivot;
}

// Goes over the columns of the matrix and returns the first pair with the same pivots.
std::tuple<int,  int> pivot_conflict (vector<column > *matrix)
{
    vector<int> pivot_entries;
    // Recompute pivots:
    for (int i = 0; i < matrix -> size(); ++i) {
        pivot_entries.push_back(pivot(&(*matrix)[i]));
    }
    for (int i = 0; i < matrix->size()-1; ++i) {
        if (pivot_entries[i] != -1) {
            int j = i + 1;
            for (; j < matrix->size(); ++j) {
                if (pivot_entries[i] == pivot_entries[j]) {
                    return std::make_tuple(i, j);
                }
            }
        }
    }
    return std::make_tuple(0, 0);
}

// Reduction algorithm: given a column a, find the indices of the columns in M such that the columns sum to a.
void reduce(column *a, vector<column> *M, vector<int> *indices)
{
    // Append a to M.
    M -> push_back(*a);
    int k = M -> size();
    vector<vector<int>> I;
    vector<int> pivot_entries;
    for (int i = 0; i < M->size(); i++) {
        pivot_entries.push_back(pivot(&((*M)[i])));
        vector<int> I_i;
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
                pivot_conflict = false;
            }
            // Resolve pivot conflict and keep track of the column that was added in I.
            if (pivot_conflict) {
                (*M)[i] ^= (*M)[conflicting_column];
                pivot_entries[i] = pivot(&((*M)[i]));
                I[i].push_back(conflicting_column);
            }
        }
        while (pivot_conflict);
    }
    M -> pop_back();
    // Find all the indices that add to the zeroing of the boundary.
    std::unordered_map<int, bool> visited;
    // Iterate over the last column of I and add all the indices that got added to this column.
    for (int i = 0; i < I[k-1].size(); i++) {
        int c = I[k-1][i];
        indices -> push_back(c);
        visited.insert({c, true});
        if (I[c].size() != 0){
            for (int j = 0; j < I[c].size(); j++) {
                if (!visited[I[c][j]]) {
                    indices -> push_back(I[c][j]);
                    visited.insert({I[c][j], true});
                }
            }
        }
    }
}

/* 
Implementation of index binary search:
*/
int index_binary_search(vector<int> &a, int x){
    int l = 0;
    int r = a.size() - 1;
    int m;
    while (l <= r) {
        m = l + (r - l) / 2;
        if (a[m] == x) {
            return m;
        }
        else if (a[m] < x) {
            l = m + 1;
        }
        else {
            r = m - 1;
        }
    }
    return m;
}


}// namespace ZZREP