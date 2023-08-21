#include "zzrep.h"

#include <algorithm>
#include <unordered_map>
#include <tuple>
#include <boost/bimap.hpp>
#include <boost/dynamic_bitset.hpp>

using namespace std;

namespace ZZREP { 

typedef boost::bimaps::bimap<vector<int>, int> SimplexIdMap;
typedef SimplexIdMap::value_type SimplexIdPair;

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
tuple<int,  int, bool> pivot_conflict (std::vector<column> *matrix);

int pivot (column *a);

void reduce(column *a, vector<column> *M, vector<int> *indices);

int find_last(column *a);

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

    for (int p = 0; p <= m; p++) {
        vector<column> Z_p;
        vector<column> C_p;
        SimplexIdMap id_p;
        vector<int> birth_timestamp_p;
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
            (id[p]).insert(SimplexIdPair(simp, j));
            // Add a row with all zeros at the end to Z[p] and C[p].
            int k = id[p].size() - 1;
            if (k != 0) {
                column last_column_C(k, 0);
                // Add all zeros to this last column:
                C[p].push_back(last_column_C);
                for (int j = 0; j < C[p].size(); j++) {
                    if (j == C[p].size() - 1) {
                        C[p][j].push_back(1);
                    }
                    else {
                        // Z[p][j].push_back(0);
                        C[p][j].push_back(0);
                        if (p == 0)  Z[p][j].push_back(0);
                    }
                }
            }
            else {
                C[p].push_back(column(1,1));
                if (p == 0)  Z[p].push_back(column(1,1));
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
                    bd_simp[id[p-1].left.at(boundary_simplex)] = 1;
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
            if (all_boundary) // simp is a boundary
            {
                /*
                An interval in dimension p gets born: 
                Append a new column simp + \sum_{a \in I} C^{p}[a] with birth timestamp i+1 to Z^p.
                */
                column new_column;
                if (p != 0)  {
                    new_column = C[p].back();
                    for (auto a: I) {
                        add_columns(&new_column, &C[p][a], &new_column);
                    }
                    Z[p].push_back(new_column);
                }
                else if (p == 0 && k != 0) {
                    column temp(C[p].back().size(), 0);
                    new_column = temp;
                    new_column[id[p].left.at(simp)] = 1;
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
                persistence->push_back(std::make_tuple(birth_timestamp[p-1][l], i, p-1, representative));
                // Set Z_{p−1}[l] = bd.(simp), C[p][l] = simp, and b^{p−1}[l] = −1.
                Z[p-1][l] = bd_simp;
                if (C[p][l].size() != 0)
                {
                    C[p][l] = C[p].back();
                }
                else {
                    C[p].push_back(C[p].back());
                }
                birth_timestamp[p-1][l] = -1;
                /*
                Avoiding pivot conflicts (column of bd.(simp) with that of another column in Z_{p-1}) as follows:
                    
                */
                tuple<int, int, bool> pivot_conflict_tuple =  pivot_conflict(&Z[p-1]);
                int a = get<0>(pivot_conflict_tuple);
                int b = get<1>(pivot_conflict_tuple);
                bool exists_pivot_conflict = get<2>(pivot_conflict_tuple);
                while (exists_pivot_conflict) {
                    column Z_pm1_aplusb;
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
        } 
        else { // Case: Arrow points backward.
            bool existence = false;
            for (auto column: Z[p]) {
                if (column[id[p].left.at(simp)]==1) {
                    existence = true;
                    break;
                }
            } 
            if (!existence)// Case: Birth.
            {   
                int a, b;
                bool boundary_pairs_bool = false;
                for (int i = 0; i < Z[p-1].size(); i++) 
                {
                    for (int j = 0; j < Z[p-1].size(); j++) 
                    {
                        if (i != j && birth_timestamp[p-1][i] < 0 && birth_timestamp[p-1][j] < 0 && (C[p][i][id[p].left.at(simp)]==1) && (C[p][j][id[p].left.at(simp)]==1)) 
                        {
                            a = i; b = j; boundary_pairs_bool = true;
                        }
                    }
                }
                while (boundary_pairs_bool) {
                    column Z_pm1_aplusb;
                    column C_p_aplusb;
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
                            if (birth_timestamp[p-1][i] < 0 && birth_timestamp[p-1][j] < 0 && (C[p][i][id[p].left.at(simp)]==1) && (C[p][j][id[p].left.at(simp)]==1)) 
                            {
                                a = i; b = j; boundary_pairs_bool = true;
                            }
                        }
                    }
                }
                // find the only column Z_{p-1}[a] with negative birth timestamp such that C[p][a] contains the simplex.
                column column;
                int x = 0;
                for (x = 0; x < Z[p-1].size(); x++) {
                    if ((birth_timestamp[p-1][x] < 0) && (C[p][x][id[p].left.at(simp)]==1)) {
                        column = Z[p-1][x];
                        break;
                    }
                }
                birth_timestamp[p-1][x] = i+1;
            }
            else // Case: Death.
            {
                // Update C[p] so that no columns contain the simplex.
                for (int a = 0; a < Z[p].size(); a++)
                {
                    if (Z[p][a][id[p].left.at(simp)]==1) {
                        for (int b = 0; b < C[p].size(); b++) 
                        {
                            if (C[p][b][id[p].left.at(simp)]==1) // each column in C[p][b] containing the simplex.
                            {    
                                add_columns(&C[p][b], &Z[p][a], &C[p][b]);
                            }
                        }
                    }
                }
                // Remove the simplex from Z[p]
                // I contains the columns that contain the simplex.
                vector<int> I;
                for (int a = 0; a < Z[p].size(); a++) // each column Z[p][a] containing the simplex
                {
                    if (Z[p][a][id[p].left.at(simp)]==1) I.push_back(a);
                }
                // sort I in the order of the birth timestamps where the order is the total order as above.
                sort(I.begin(), I.end(), [&](int &a, int &b){ return ((I[a] == b) || ((I[a] < b) && (filt_op[I[b]-1])) || ((I[a] > I[b]) && (!(filt_op[I[a]-1]))));});
                int alpha = I[0];
                I.erase(I.begin());
                column z = Z[p][alpha];
                for (auto a: I)
                {
                    if (pivot(&Z[p][a]) > pivot(&z)) {
                        add_columns(&Z[p][a], &z, &Z[p][a]);
                    }   
                    else 
                    {
                        column temp = Z[p][a];
                        add_columns(&Z[p][a], &z, &Z[p][a]);
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
                // Delete the column Z[p][I[0]] from Z[p] and delete birth_timestamp_p[I[0]] from birth_timestamp_p.
                Z[p].erase(Z[p].begin() + alpha);
                birth_timestamp[p].erase(birth_timestamp[p].begin() + alpha);
                /* 
                UPDATE: There is no point in clearing this. This may cause memory overload but allows for bitwise operations.
                // Remove the zero row and update id as well.
                for (unsigned i = 0; i < C[p].size(); ++i)
                {
                if (C[p][i].size() > id[p].left.at(simp))
                {
                    C[p][i].erase(C[p][i].begin() + id[p][simp]);
                }
                }
                for (unsigned i = 0; i < Z[p].size(); ++i)
                {
                if (Z[p][i].size() > id[p].left.at(simp))
                {
                    Z[p][i].erase(Z[p][i].begin() + id[p][simp]);
                }
                }
                // Iterate over id[p] and decrease the indices of the simplices that have index greater than id[simp].
                SimplexIdMap::right_iterator it = id[p].right.lower_bound(id[p].left.at(simp));
                for (; it != id[p].right.end(); ++it)
                {
                    id[p].right.modify_key( it, it->first - 1 );
                }
                id[p].erase(simp);
                */
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
    /*
    typedef boost::dynamic_bitset<>::size_type size_type;
    const size_type npos = boost::dynamic_bitset<>::npos;
    size_type first_idx = (*a).find_first();
    size_type current_idx = first_idx;
    if (first_idx != npos)
    {
      do {
         current_idx = (*a).find_next(current_idx);
      } while ((*a).find_next(current_idx) != boost::dynamic_bitset<>::npos);
    }
    return current_idx;
    */
   return ((*a).size() - 1) - (*a).find_first();
}

// Goes over the columns of the matrix and returns the first pair with the same pivots.
std::tuple<int,  int, bool> pivot_conflict (vector<column > *matrix)
{
    vector<int> pivot_entries;
    // Recompute pivots:
    for (int i = 0; i < matrix -> size(); i++) {
        pivot_entries.push_back(pivot(&(*matrix)[i]));
    }
    for (int i = 0; i < matrix->size()-1; i++) {
        if (pivot_entries[i] != -1) {
            for (int j = i+1; j < matrix->size(); j++) {
                if (pivot_entries[i] == pivot_entries[j]) {
                    return std::make_tuple(i, j, true);
                }
            }
        }
    }
    return std::make_tuple(-1, -1, false);
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

}// namespace ZZREP