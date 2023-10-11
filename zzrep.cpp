#include "zzrep.h"

#include <algorithm>
#include <unordered_map>
#include <tuple>
#include <boost/bimap.hpp>
#include <boost/dynamic_bitset.hpp>
#include <vector>

using namespace std;

namespace ZZREP { 

/* HELPER DATA STRUCTURES: */
typedef boost::bimaps::bimap<vector<int>, int> SimplexIdMap;
typedef boost::bimaps::bimap<int, int> BdToChainMap;
typedef SimplexIdMap::value_type SimplexIdPair;
typedef BdToChainMap::value_type BdToChainPair;
typedef boost::dynamic_bitset<> column;

/* DECLARATION OF HELPER FUNCTIONS: */
int pivot (column *a);
tuple<int,  int> pivot_conflict (std::vector<column> *matrix);
void reduce(column *a, vector<column> *M, vector<int> *indices);

void ZigzagRep::compute(
        const std::vector<vector<int> > &filt_simp, 
        const std::vector<bool> &filt_op,
        std::vector <tuple <int, int, int, std::vector<std::vector<int> > > > *persistence, int m) {
    
    persistence -> clear();
    int n = filt_simp.size();
    /*
    DECLARATION of the data structures:
    */
    vector<SimplexIdMap> id(m); //  A unique integral id for each p-simplex.
    vector<vector<column> > Z(m); // A cycle matrix Z[p] for each dimension p that is maintained such that each column of Z[p] represents a p-cycle in K_i.
    vector<vector<column> > C(m); // A chain matrix C[p] for each dimension p so that the chains in C[p] are associated with the the boundary cycles in Z[p-1].
    vector<vector<int> > birth_timestamp(m); // For each column Z[p][j], a birth timestamp is maintained, where a negative value implies that the column pertains to a boundary.
    vector<BdToChainMap> bd_to_chain(m); // Maps between indices of Z[p-1] and C[p] such that \partial C[p+1][j] = Z[p][bd_to_chain[j]] for all j such that birth_timestamp[p-1][j] < 0.
    /*
    COMPUTATION of zigzag persistence:
    */
    for (int i = 0; i < n; ++i) {
        const vector<int> &simp = filt_simp[i];
        int p = simp.size() - 1; // p denotes the dimension of the simplex.
        if (filt_op[i]) { // INSERTION
            // A p-simplex is inserted into id[p].
            int k = id[p].size();
            (id[p]).insert(SimplexIdPair(simp, k));
            // Next, we add a row of all zeros to Z[p] and C[p] for indicating the insertion of a new p-simplex which increases the respective ranks by 1.
            if (k != 0) {
                for (int j = 0; j < C[p].size(); j++) {
                    C[p][j].push_back(0);
                }
                for (int j = 0; j < Z[p].size(); ++j) {
                    Z[p][j].push_back(0);
                }
            }
            // Represent the boundary of simp as a sum of columns of Z_{p-1} by a reduction algorithm; I is such set of columns.
            bool all_boundary = true;
            column bd_simp;
            vector<int> I;
            if (p != 0) {
                // Compute boundary of the simplex:
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
                // Check the birth timestamps to check whether all of them are boundaries.
                for (auto a: I) {
                    if (birth_timestamp[p-1][a] >= 0) {
                        all_boundary = false;
                        break;
                    }
                }
            }
            if (all_boundary) // The inserted simplex covers a (p-1)-boundary (an interval in dimension p gets born) FIXME: Do we need to update the cycle to chain map here?
            {
                /*
                Append a new column simp + \sum_{a \in I} C^{p}[a] with birth timestamp i+1 to Z^p.
                */
                column new_column = column(id[p].size(), 0); // New column with of size id[p].size() with 1 appended at the end.
                new_column[id[p].left.at(simp)] = 1;
                if (p != 0)  {
                    for (auto a: I) {
                        int chain_a = bd_to_chain[p-1].left.at(a);
                        new_column ^= C[p][chain_a];
                    }
                }
                Z[p].push_back(new_column);
                birth_timestamp[p].push_back(i+1);
            }
            else { // The inserted simplex kills a (p-1)-cycle (an interval in dimension (p-1) dies).
                /*
                Let J consist of indices in I whose corresponding columns in Z[p−1] have non-negative birth timestamps. 
                */ 
                vector<int> J;
                for (int a: I) {
                    if (birth_timestamp[p-1][a] >= 0) { // Gather all the non-boundary cycles.
                        J.push_back(a);
                    }
                }
                sort(J.begin(), J.end());
                // Check if arrow at b^{p−1}[c]−1 points backward for all c in J
                bool arrow_backward = true;
                for (auto c: J) {
                    if (filt_op[birth_timestamp[p-1][c]-1]) {
                        arrow_backward = false;
                        break;
                    }
                }
                // If the arrow at {b^{p−1}[c]−1} points backward for all c in J, let l be the smallest index in J. Otherwise, let l be the largest c in J such that the arrow at {b^{p−1}[c]−1} points forward.
                int l;
                if (arrow_backward) {
                    l = *J.begin();
                }
                else {
                    for (auto i = J.size()-1; i >= 0; i--) {
                        if (filt_op[birth_timestamp[p-1][J[i]]-1]) {
                            l = J[i];
                            break;
                        }
                    }
                }
                // Output the (p − 1)-th interval [b^{p−1}[l], i] after gathering the relevant representative.
                vector<vector<int>> representative;
                for (int i = 0; i < Z[p-1][l].size(); i++) {
                    if (Z[p-1][l][i] == 1) {
                        representative.push_back(id[p-1].right.at(i));
                    }
                }
                persistence->push_back(std::make_tuple(birth_timestamp[p-1][l], i, p-1, representative));
                /* 
                UPDATE:
                */
                // Set Z[p−1][l] = bd_simp.
                Z[p-1][l] = bd_simp;
                // Set C[p][l] = simp and update the boundary-to-chain map.
                column current_chain(id[p].size(), 0);
                current_chain[id[p].left.at(simp)] = 1;
                C[p].push_back(current_chain);
                (bd_to_chain[p-1]).insert(BdToChainPair(l, C[p].size()-1));
                // Set b^{p−1}[l] = −1 and record l as a boundary cycle index.
                birth_timestamp[p-1][l] = -1;
                /*
                Avoiding pivot conflicts (column of bd_simp with that of another column in Z[p-1]) as follows:
                */
                // TODO: Intially the only pivot conflict can happen with the boundary simplex so we can just check for that.
                tuple<int, int> pivot_conflict_tuple =  pivot_conflict(&Z[p-1]);
                int a = get<0>(pivot_conflict_tuple);
                int b = get<1>(pivot_conflict_tuple);
                bool exists_pivot_conflict = (a != b);
                while (exists_pivot_conflict) {
                    column Z_pm1_aplusb;
                    Z_pm1_aplusb = Z[p-1][a] ^ Z[p-1][b];
                    if (birth_timestamp[p-1][a] < 0 && birth_timestamp[p-1][b] < 0) {
                        Z[p-1][a] = Z_pm1_aplusb;
                        // Update chain matrices to reflect this addition (use the BdToChainMap):
                        int chain_a = bd_to_chain[p-1].left.at(a);
                        int chain_b = bd_to_chain[p-1].left.at(b);
                        C[p][chain_a] = C[p][chain_a] ^ C[p][chain_b];
                        // FIXME: Check if we need to update the cycle to chain map:
                        // bd_to_chain[p-1].erase(b);
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
        else { // DELETION:
            // Find the index of the simplex in id[p].
            bool existence = false;
            int idx = id[p].left.at(simp);
            // Check if the simplex exists in Z[p].
            for (auto column: Z[p]) {
                if (column[idx]==1) {
                    existence = true;
                    break;
                }
            } 
            if (!existence)// The deleted simplex does not constitute a cycle and its boundary now becomes a (p-1)-cycle (An interval in dimension (p-1) gets born).
            {   
                // FIXME: There's obviously a bug here; make sure that the cyle to chain indices are getting retrieved correctly.
                int a, b;
                bool boundary_pairs_bool = false;
                for (int x = 0; x < bd_to_chain[p-1].size(); x++) 
                {
                    for (int y = x + 1; y < bd_to_chain[p-1].size(); y++) 
                    {
                        int cycle_x = bd_to_chain[p-1].left[x];
                        int cycle_y = bd_to_chain[p-1].left[y];
                        int chain_x = (bd_to_chain[p-1]).left.at(cycle_x);
                        int chain_y = (bd_to_chain[p-1]).left.at(cycle_y);
                        if (C[p][chain_x][idx]==1)
                        {
                            if (C[p][chain_y][idx]==1)
                            {
                                a = cycle_x; b = cycle_y; boundary_pairs_bool = true;
                            }
                        }
                    }
                }
                while (boundary_pairs_bool) {
                    column Z_pm1_aplusb;
                    column C_p_aplusb;
                    Z_pm1_aplusb = Z[p-1][a] ^ Z[p-1][b];
                    int chain_a = (bd_to_chain[p-1]).left.at(a);
                    int chain_b = (bd_to_chain[p-1]).left.at(b);
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
                    for (int x = 0; x < bd_to_chain[p-1].size(); x++) 
                    {
                        for (int y = x + 1; y < bd_to_chain[p-1].size(); y++) 
                        {
                            int cycle_x = bd_to_chain[p-1].left[x];
                            int cycle_y = bd_to_chain[p-1].left[y];
                            int chain_x = (bd_to_chain[p-1]).left.at(cycle_x);
                            int chain_y = (bd_to_chain[p-1]).left.at(cycle_y);
                            if (C[p][chain_x][idx]==1)
                            {
                                if (C[p][chain_y][idx]==1)
                                {
                                    a = cycle_x; b = cycle_y; boundary_pairs_bool = true;
                                }
                            }
                        }
                    }
                }
                // find the only column Z[p-1][a] with negative birth timestamp such that C[p][a] contains the simplex.
                int alpha;
                for (int b_idx = 0; b_idx < bd_to_chain[p-1].size(); b_idx++) {
                    alpha = bd_to_chain[p-1].left[b_idx];
                    int chain_only_idx = bd_to_chain[p-1].left.at(alpha);
                    if (C[p][chain_only_idx][idx] == 1) {
                        C[p][chain_only_idx][idx] = 0;
                        break;
                    }
                }
                birth_timestamp[p-1][alpha] = i+1;
                // Remove the boundary to chain record.
                bd_to_chain[p-1].left.erase(alpha);
            }
            else // // The deleted simplex is part of some p-cycle (An interval in dimension p dies).
            {
                // Update C[p] so that no columns contain the simplex.
                for (int a = 0; a < Z[p].size(); a++)
                {
                    if (Z[p][a][idx]==1) {
                        for (int b = 0; b < C[p].size(); b++) 
                        {
                            if (C[p][b][idx] == 1) // each column in C[p][chain_b] containing the simplex s.t. Z[p-1][b] is a boundary. // FIXME: Why are we only iterating over these?
                            {  
                                int cycle_b = bd_to_chain[p-1].right.at(b); // FIXME: Doesn't this assume that we only iterate over the chains whose boundary is the simplex? So, we can just iterate over the boundary cycle indices?
                                if (birth_timestamp[p-1][cycle_b] < 0){  
                                        C[p][b] ^= Z[p][a];
                                }
                            }
                        }
                    }
                }
                // Remove the simplex from Z[p]
                vector<int> I; // Gather indices of columns that contain the simplex.
                for (int a = 0; a < Z[p].size(); a++)
                {
                    if (Z[p][a][idx]==1) I.push_back(a);
                }
                // sort I in the order of the birth timestamps where the order is the total order as above.
                sort(I.begin(), I.end(), [&](int &a, int &b){ return ((I[a] == b) || ((I[a] < b) && (filt_op[I[b]-1])) || ((I[a] > I[b]) && (!(filt_op[I[a]-1]))));});
                // The column to be deleted is the first column in I.
                int alpha = I[0];
                I.erase(I.begin());
                // Iterate over the rest of I and ensure that the pivots remain distinct:
                column z = Z[p][alpha];
                for (auto a: I)
                {
                    if (pivot(&Z[p][a]) > pivot(&z)) 
                    {
                        Z[p][a] = Z[p][a] ^ z; // FIXME: Don't we need to update the cycle to chain map?
                    }   
                    else 
                    {
                        column temp = Z[p][a];
                        Z[p][a] = Z[p][a] ^ z;
                        z = temp;
                    }
                }
                // Output the p-th interval [birth_timestamp_p[alpha], i] after gathering the representatives.
                vector<vector<int>> representative;
                for (int i = 0; i < Z[p][alpha].size(); i++) {
                    if (Z[p][alpha][i] == 1) {
                        representative.push_back(id[p].right.at(i));
                    }
                }
                persistence->push_back(std::make_tuple(birth_timestamp[p][alpha], i, p, representative));
                /*
                UPDATE:
                */ 
                Z[p].erase(Z[p].begin() + alpha); // Delete the column Z[p][alpha] from Z[p]
                birth_timestamp[p].erase(birth_timestamp[p].begin() + alpha); // Delete the birth_timestamp[p][alpha] from birth_timestamp[p].
            }
        }
    }
    /*
     POST-PROCESSING:
    */
    // For each p and each column Z[p][a] of Z[p] with non-negative birth timestamp, output the p-th interval [birth_timestamp[p][a], n].
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

/*
IMPLEMENTATION OF HELPER FUNCTIONS:
*/
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
}// namespace ZZREP