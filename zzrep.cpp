#include "zzrep.h"

#include <algorithm>
#include <map>
#include <unordered_map>
#include <tuple>
#include <boost/bimap.hpp>
#include <boost/dynamic_bitset.hpp>
#include <vector>

using namespace std;

namespace ZZREP { 

/* HELPER DATA STRUCTURES: */
typedef boost::dynamic_bitset<> column;
typedef boost::bimaps::bimap<vector<int>, int> SimplexIdMap;
typedef std::map<int, int> BdToChainMap;
typedef std::unordered_map<int, int> CycleToBundleMap;
typedef SimplexIdMap::value_type SimplexIdPair;
typedef BdToChainMap::value_type BdToChainPair;
typedef CycleToBundleMap::value_type CycleToBundlePair;

/* DECLARATION OF HELPER FUNCTIONS: */
int pivot (column *a);
void reduce(column *a, vector<column> *M, vector<int> *indices, vector<int> *pivots);

void ZigzagRep::compute(
        const std::vector<vector<int> > &filt_simp, 
        const std::vector<bool> &filt_op,
        std::vector <tuple <int, int, int, std::vector<std::vector<int> > > > *persistence, int m) {
    
    persistence -> clear();
    int n = filt_simp.size();
    /*
    DECLARATION of the data structures:
    */
    vector<SimplexIdMap> id; //  A unique integral id for each p-simplex.
    vector<vector<column> > Z; // A cycle matrix Z[p] for each dimension p that is maintained such that each column of Z[p] represents a p-cycle in K_i.
    vector<vector<column> > C; // A chain matrix C[p] for each dimension p so that the chains in C[p] are associated with the the boundary cycles in Z[p-1].
    vector<vector<int> > birth_timestamp; // For each column Z[p][j], a birth timestamp is maintained, where a negative value implies that the column pertains to a boundary.
    vector<BdToChainMap> bd_to_chain; // Maps between indices of Z[p-1] and C[p] such that \partial C[p+1][j] = Z[p][bd_to_chain[j]] for all j such that birth_timestamp[p-1][j] < 0.
    vector<vector<int>> pivots; // For each column Z[p][j] we have a unique pivot entry; we then store pivots[p][i] = j where j = -1 implies that no column exists in Z[p-1] with pivot i.
    // Data structures for storing the intervals and cycles in order to compute representatives: 
    /*
    vector<vector<column> > bundle; // A collection of wires: bundle[p] stores p-dimensional wires.
    vector<vector<int> > timestamps; // Wire timestamps.
    vector<CycleToBundleMap> links; // Map from the cycle matrix to the wires: links[p][j] is the collection of wires associated with Z[p][j].
    */
    /*
    INITIALIZATION:
    */
    for (int d = 0; d <= m; ++d)
    {
        SimplexIdMap id_p; 
        vector<column> Z_p; 
        vector<column> C_p; 
        vector<int> birth_timestamp_p;
        BdToChainMap bd_to_chain_p; 
        vector<int> pivots_p;
        id.push_back(id_p);
        Z.push_back(Z_p);
        C.push_back(C_p);
        birth_timestamp.push_back(birth_timestamp_p);
        bd_to_chain.push_back(bd_to_chain_p);
        pivots.push_back(pivots_p);
        /*
        vector<column> bundle_p;
        vector<int> timestamps_p;
        CycleToBundleMap links_p;
        bundle.push_back(bundle_p);
        timestamps.push_back(timestamps_p);
        links.push_back(links_p);
        */
    }
    /*
    COMPUTATION of zigzag persistence:
    */
    for (int i = 0; i < n; ++i) {
        const vector<int> &simp = filt_simp[i];
        int p = simp.size() - 1; // p denotes the dimension of the simplex.
        if (i == 428) {
            cout << "i = " << i << endl;
            cout << "p = " << p << endl;
        }
        if (filt_op[i]) { // INSERTION
            // A p-simplex is inserted into id[p].
            int k = id[p].size();
            (id[p]).insert(SimplexIdPair(simp, k));
            // Next, we add a row of all zeros to Z[p] and C[p] for indicating the insertion of a new p-simplex which increases the respective ranks by 1.
            if (k != 0) 
            {
                for (int j = 0; j < C[p].size(); j++) {
                    C[p][j].push_back(0);
                }
                for (int j = 0; j < Z[p].size(); ++j) {
                    Z[p][j].push_back(0);
                }
            }
            pivots[p].push_back(-1);
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
                // Find the columns in Z[p-1] that sum to bd_simp:
                column bd_simp_temp = bd_simp;  // Copy bd_simp to a temporary column.
                int pivot_bd = pivot(&bd_simp_temp);
                bool zeroed = (pivot_bd == -1);
                // FIXME: This routine does not work somehow, it gives duplicate indices at i = 428.
                while (!zeroed) {
                    // Find the column in M that has the same pivot as a.
                    int conflict_idx = pivots[p-1][pivot_bd];
                    // Add this column to indices.
                    bd_simp_temp ^= Z[p-1][conflict_idx];
                    I.push_back(conflict_idx);
                    // Update the pivot of a.
                    pivot_bd = pivot(&bd_simp_temp);
                    zeroed = (pivot_bd == -1);
                } 
                // Check the birth timestamps to check whether all of them are boundaries.
                for (auto a: I) {
                    if (birth_timestamp[p-1][a] >= 0) {
                        all_boundary = false;
                        break;
                    }
                }
            }
            if (all_boundary) // The inserted simplex covers a (p-1)-boundary (an interval in dimension p gets born)
            {
                /*
                Append a new column simp + \sum_{a \in I} C^{p}[a] with birth timestamp i+1 to Z^p.
                */
                column new_column = column(id[p].size(), 0); // New column with of size id[p].size() with 1 appended at the end.
                new_column[id[p].left.at(simp)] = 1;
                if (p != 0)  {
                    for (auto a: I) {
                        int chain_a = bd_to_chain[p-1].at(a);
                        new_column ^= C[p][chain_a];
                    }
                }
                Z[p].push_back(new_column);
                birth_timestamp[p].push_back(i+1);
                int pivot_new = pivot(&new_column);
                pivots[p][pivot_new] = Z[p].size() - 1;
                // Add a wire to the bundle with timestamp 1 and update the link from the cycle to the bundle.
                /*
                bundle[p].push_back(new_column);
                timestamps[p].push_back(i+1);
                links[p].insert(CycleToBundlePair(Z[p].size()-1, bundle[p].size()-1));
                */
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
                // TODO: Output the (p − 1)-th interval [b^{p−1}[l], i] after gathering the relevant representative.
                // Here, we need to gather all the wires pertaining to the cycle l and reconstruct the representatives at each index.
                /* 
                UPDATE:
                */
                // Set Z[p−1][l] = bd_simp.
                Z[p-1][l] = bd_simp;
                birth_timestamp[p-1][l] = -1; // Set b^{p−1}[l] = −1 and record l as a boundary cycle index.
                // Set C[p][l] = simp and update the boundary-to-chain map.
                column current_chain(id[p].size(), 0);
                current_chain[id[p].left.at(simp)] = 1;
                C[p].push_back(current_chain);
                (bd_to_chain[p-1]).insert(BdToChainPair(l, C[p].size()-1));
                /*
                Avoiding pivot conflicts (column of bd_simp with that of another column in Z[p-1]).
                We know that the pivots of Z[p-1] were unique before the insertion of bd_simp. 
                Thus, we need to ensure that they remain unique after the insertion of bd_simp. 
                */
                // TODO: Update the links as we add columns in Z[p-1].
                int pivot_l = pivot(&Z[p-1][l]);
                int current_idx = pivots[p-1][pivot_l];
                bool pivot_conflict;
                int a, b;
                if (current_idx  == l or current_idx == -1)
                {   
                    pivot_conflict = false;
                    pivots[p-1][pivot_l] = l;
                }
                else {
                    pivot_conflict = true;
                    a = l;
                    b = current_idx;
                }
                while (pivot_conflict)
                {
                    if (birth_timestamp[p-1][a] < 0 && birth_timestamp[p-1][b] < 0) {
                        Z[p-1][a] =Z[p-1][a] ^ Z[p-1][b];
                        // Update chain matrices to reflect this addition (use the BdToChainMap):
                        int chain_a = bd_to_chain[p-1].at(a);
                        int chain_b = bd_to_chain[p-1].at(b);
                        C[p][chain_a] = C[p][chain_a] ^ C[p][chain_b];
                        // Check for pivot conflicts again (pivot of b remains the same, but the pivot of a might have changed):
                        int pivot_a = pivot(&Z[p-1][a]);
                        int current_idx = pivots[p-1][pivot_a];
                        if (current_idx == a or current_idx == -1)
                        {
                            pivot_conflict = false;
                            pivots[p-1][pivot_a] = a;
                        }
                        else {
                            pivot_conflict = true;
                            b = current_idx;
                        }

                    }
                    else if (birth_timestamp[p-1][a] < 0 && birth_timestamp[p-1][b] >= 0) {
                        Z[p-1][b] = Z[p-1][a] ^ Z[p-1][b];
                        // No need to update the chain matrices as a boundary is added to a cycle. Instead, check for pivot conflicts again:
                        int pivot_b = pivot(&Z[p-1][b]);
                        int current_idx = pivots[p-1][pivot_b];
                        if (current_idx == b or current_idx == -1)
                        {
                            pivot_conflict = false;
                            pivots[p-1][pivot_b] = b;
                        }
                        else {
                            pivot_conflict = true;
                            a = current_idx;
                        }
                    }
                    else if (birth_timestamp[p-1][a] >= 0 && birth_timestamp[p-1][b] < 0) {
                        Z[p-1][a] = Z[p-1][a] ^ Z[p-1][b];
                        // Again, no need to update the chain matrices as a boundary is added to a cycle. Instead, check for pivot conflicts again:
                        int pivot_a = pivot(&Z[p-1][a]);
                        int current_idx = pivots[p-1][pivot_a];
                        if (current_idx == a or current_idx == -1)
                        {
                            pivot_conflict = false;
                            pivots[p-1][pivot_a] = a;
                        }
                        else {
                            pivot_conflict = true;
                            b = current_idx;
                        }
                    }
                    else {
                        bool a_lessthan_b = ((a == b) || ((a < b) && (filt_op[b-1])) || ((a > b) && (!(filt_op[a-1]))));   
                        if (a_lessthan_b) {
                            Z[p-1][b] = Z[p-1][a] ^ Z[p-1][b];
                            // Since both are cycles, no need to update the chain matrices. Instead, check for pivot conflicts again:
                            int pivot_b = pivot(&Z[p-1][b]);
                            int current_idx = pivots[p-1][pivot_b];
                            if (current_idx == b or current_idx == -1)
                            {
                                pivot_conflict = false;
                                pivots[p-1][pivot_b] = b;
                            }
                            else {
                                pivot_conflict = true;
                                a = current_idx;
                            }
                        }
                        else {
                            Z[p-1][a] = Z[p-1][a] ^ Z[p-1][b];
                            // Since both are cycles, no need to update the chain matrices. Instead, check for pivot conflicts again:
                            int pivot_a = pivot(&Z[p-1][a]);
                            int current_idx = pivots[p-1][pivot_a];
                            if (current_idx == a or current_idx == -1)
                            {
                                pivot_conflict = false;
                                pivots[p-1][pivot_a] = a;
                            }
                            else {
                                pivot_conflict = true;
                                b = current_idx;
                            }
                        }
                    }
                    /*
                    bundle[p-1].push_back(bd_simp);
                    timestamps[p-1].push_back(i+1);
                    links[p-1].insert(CycleToBundlePair(l, bundle[p-1].size()-1));
                    */
                }
            }
        } 
        else { // DELETION:
            // Find the index of the simplex in id[p].
            bool existence = false;
            int idx = id[p].left.at(simp);
            column simp_column = column(id[p].size(), 0);
            // Check if the simplex exists in Z[p].
            for (auto p_column: Z[p]) {
                if (p_column[idx]==1) {
                    existence = true;
                    simp_column = p_column;
                    break;
                }
            } 
            if (!existence)// The deleted simplex does not constitute a cycle and its boundary now becomes a (p-1)-cycle (An interval in dimension (p-1) gets born).
            {   
                // TODO: Update the links as we add columns in Z[p-1].
                int a , b, chain_a, chain_b;
                // Handling the edge case where there is only one boundary pair.
                a = bd_to_chain[p-1].begin()->first;
                chain_a = bd_to_chain[p-1].begin()->second;
                bool boundary_pairs_bool = false;
                typedef BdToChainMap::const_iterator bd_const_iterator;
                for(bd_const_iterator bd_iter = bd_to_chain[p-1].begin();
                    bd_iter != bd_to_chain[p-1].end(); ++bd_iter)
                {
                    for(bd_const_iterator bd_iter_b = bd_to_chain[p-1].begin();
                    bd_iter_b != bd_to_chain[p-1].end(); ++bd_iter_b)
                    {
                        if (bd_iter != bd_iter_b) {
                            int cycle_x = bd_iter->first;
                            int cycle_y = bd_iter_b->first;
                            int chain_x = bd_iter->second;
                            int chain_y = bd_iter_b->second;
                            if (C[p][chain_x][idx]==1)
                            {
                                a = cycle_x;
                                chain_a = chain_x;
                                if (C[p][chain_y][idx]==1)
                                {
                                    b = cycle_y;
                                    chain_b = chain_y;
                                    boundary_pairs_bool = true;
                                    break;
                                }
                            }
                        }
                    }                    
                }
                while (boundary_pairs_bool) {
                    int pivot_a = pivot(&Z[p-1][a]);
                    int pivot_b = pivot(&Z[p-1][b]);
                    if (pivot_a > pivot_b) 
                    {
                        Z[p-1][a] = Z[p-1][a] ^ Z[p-1][b];
                        C[p][chain_a] = C[p][chain_a] ^ C[p][chain_b];
                    }
                    else 
                    {
                        Z[p-1][b] = Z[p-1][a] ^ Z[p-1][b];
                        C[p][chain_b] = C[p][chain_a] ^ C[p][chain_b];
                    }
                    boundary_pairs_bool = false;
                    for(bd_const_iterator bd_iter = bd_to_chain[p-1].begin();
                    bd_iter != bd_to_chain[p-1].end(); ++bd_iter)
                    {
                        for(bd_const_iterator bd_iter_b = bd_to_chain[p-1].begin();
                        bd_iter_b != bd_to_chain[p-1].end(); ++bd_iter_b)
                        {
                            if (bd_iter != bd_iter_b) {
                                int cycle_x = bd_iter->first;
                                int cycle_y = bd_iter_b->first;
                                int chain_x = bd_iter->second;
                                int chain_y = bd_iter_b->second;
                                if (C[p][chain_x][idx]==1)
                                {
                                    a = cycle_x;
                                    chain_a = chain_x;
                                    if (C[p][chain_y][idx]==1)
                                    {
                                        b = cycle_y;
                                        chain_b = chain_y;
                                        boundary_pairs_bool = true;
                                        break;
                                    }
                                }
                            }
                        }                    
                    }
                }
                // find the only column Z[p-1][a] with negative birth timestamp such that C[p][chain_a] contains the simplex.
                int alpha = a;
                // We do not remove the column as that means changing the entire bd_to_chain map. Instead, we can zero out the column chain_a from C[p].
                for (int i = 0; i < C[p][chain_a].size(); i++) {
                    if (C[p][chain_a][i] == 1) {
                        C[p][chain_a][i] = 0;
                    }
                }
                birth_timestamp[p-1][alpha] = i+1;
                // Remove the boundary to chain record by finding the iterator and removing it.
                bd_const_iterator bd_iter = bd_to_chain[p-1].find(alpha);
                bd_to_chain[p-1].erase(bd_iter);
                // Add a wire to the bundle with timestamp i+1 and update the link from the cycle to the bundle.
                /*
                bundle[p-1].push_back(Z[p-1][alpha]);
                timestamps[p-1].push_back(i+1);
                // Update the links from a to the bundle.
                links[p-1][alpha] = bundle[p-1].size()-1;
                */
            }
            else // // The deleted simplex is part of some p-cycle (An interval in dimension p dies).
            {
                // Update C[p] so that no columns contain the simplex.
                for (int b = 0; b < C[p].size(); b++) 
                {
                    if (C[p][b][idx] == 1) // each column in C[p][chain_b] containing the simplex s.t. Z[p-1][b] is a boundary.
                    {  
                        C[p][b] ^= simp_column; // Here, we do not need to check whether Z[p-1][b] is a boundary since we only store the chains whose boundaries are non-zero.
                    }
                }
                // Remove the simplex from Z[p]
                vector<int> I; // Gather indices of columns that contain the simplex. 
                for (int a = 0; a < Z[p].size(); a++)
                {
                    if (Z[p][a][idx]==1) I.push_back(a);
                }
                // sort I in the order of the birth timestamps where the order is the total order as above.
                sort(I.begin(), I.end(), [&](int &a, int &b){ return (((a == b) || ((a < b) && (filt_op[b-1])) || ((a > b) && (!(filt_op[a-1])))));});
                // The column to be deleted is the first column in I.
                int alpha = I[0];
                I.erase(I.begin());
                // Iterate over the rest of I and ensure that the pivots remain distinct:
                column z = Z[p][alpha];
                int z_pivot = pivot(&z);
                int alpha_pivot = z_pivot;
                // TODO: Update the links as we add columns in Z[p-1].
                for (auto a: I)
                {
                    int a_pivot = pivot(&Z[p][a]);
                    if (a_pivot > z_pivot) 
                    {
                        Z[p][a] = Z[p][a] ^ z; // This does not change the pivot of Z[p][a].
                    }   
                    else 
                    {
                        column temp = Z[p][a];
                        int temp_pivot = a_pivot;
                        Z[p][a] = Z[p][a] ^ z;
                        a_pivot = z_pivot;
                        z = temp;
                        z_pivot = temp_pivot;
                        // Update the pivot of Z[p][a] in pivots[p].
                        pivots[p][a_pivot] = a;
                    }
                }
                // Output the p-th interval [birth_timestamp_p[alpha], i] after gathering the representatives.
                // TODO: Gather the representatives at each index using the links from alpha and remove the record.
                /*
                UPDATE:
                */ 
                Z[p].erase(Z[p].begin() + alpha); // Delete the column Z[p][alpha] from Z[p]
                pivots[p][alpha_pivot] = -1; // Remove the pivot of Z[p][alpha] from pivots[p].
                birth_timestamp[p].erase(birth_timestamp[p].begin() + alpha); // Delete the birth_timestamp[p][alpha] from birth_timestamp[p].
                /*
                 This will have a cascading effect on the index records of pivot and bd_to_chain.
                */
                // Update the pivot map:
                for (int i = 0; i < pivots[p].size(); i++) {
                    if (pivots[p][i] > alpha) {
                        pivots[p][i] -= 1;
                    }
                }
                // Update the bd_to_chain map by traversing the map and decrementing keys that are greater than alpha:
                for (auto it = bd_to_chain[p].begin(); it != bd_to_chain[p].end(); ++it) {
                    if (it->first > alpha) {
                        int new_key = it->first - 1;
                        int new_value = it->second;
                        bd_to_chain[p].erase(it);
                        bd_to_chain[p].insert(BdToChainPair(new_key, new_value));
                    }
                }
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
}// namespace ZZREP