#include "zzrep.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <unordered_map>
#include <tuple>
#include <boost/bimap.hpp>
#include <boost/dynamic_bitset.hpp>
#include <utility>
#include <vector>

using namespace std;

namespace ZZREP { 

/* HELPER DATA STRUCTURES: */
template <class ElemType>
class VecHash { 
public:
    size_t operator()(const std::vector<ElemType>& v) const; 
};

template <class ElemType>
size_t VecHash<ElemType>
    ::operator()(const std::vector<ElemType>& v) const {

    std::size_t seed = 0;

    for (auto e : v) { boost::hash_combine(seed, e); }

    return seed;
}

template <class ElemType>
class VecEqual { 
public:
    bool operator()(const std::vector<ElemType>& v1, 
        const std::vector<ElemType>& v2) const; 
};

template <class ElemType>
bool VecEqual<ElemType>
    ::operator()(const std::vector<ElemType>& v1, 
        const std::vector<ElemType>& v2) const {

    if (v1.size() != v2.size()) { return false; }

    for (auto i = 0; i < v1.size(); i ++) {
        if (v1[i] != v2[i]) {
            return false;
        }
    }

    return true;
}

class birth_timestamp {
public:
    // If bd is true if the index is a boundary and val_idx records to the index of the corresponding chain, while if bd is false, then val_idx records the birth timestamp of the cycle.
    bool bd;
    int val_idx;
    birth_timestamp(bool bd, int val_idx) : bd(bd), val_idx(val_idx) {} 
};

typedef boost::dynamic_bitset<> bitset;
typedef shared_ptr<bitset> column;
typedef std::map< vector<int>, int> SimplexIdMap;
typedef std::map<int, int> PivotMap;
typedef SimplexIdMap::value_type SimplexIdPair;
typedef PivotMap::value_type PivotPair;

/* DECLARATION OF HELPER FUNCTION: */
int pivot(column a);

void dynamic_xor(column a, column b);

void ZigzagRep::compute(
        const std::vector<vector<int> > &filt_simp, 
        const std::vector<bool> &filt_op,
        std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > *persistence, int m) {
    
    persistence -> clear();
    int n = filt_simp.size();
    /*
    DECLARATION:
    */
    vector<SimplexIdMap> id; //  The current integral id for each p-simplex.
    vector<vector<int>> unique_id; // The unique integral id for each p-simplex.
    vector<vector<column> > Z; // A cycle matrix Z[p] for each dimension p that is maintained such that each column of Z[p] represents a p-cycle in K_i.
    vector<vector<int>> available_columns_Z; // A list of available columns in Z[p] for each dimension p.
    vector<vector<int>> used_columns_Z; // A list of used columns in Z[p] for each dimension p.
    vector<vector<column> > C; // A chain matrix C[p] for each dimension p so that the chains in C[p] are associated with the the boundary cycles in Z[p-1].
    vector<vector<int>> available_columns_C; // A list of available columns in Z[p] for each dimension p.
    vector<vector<int>> used_columns_C; // A list of used columns in Z[p] for each dimension p.
    vector<vector<pair <int, int> > > birth_timestamp; // For each column Z[p][j], a birth timestamp is maintained, where a negative value implies that the column pertains to a boundary.
    vector<PivotMap> pivots; // For each column Z[p][j] we have a unique pivot entry; we then store pivots[p][i] = j where j = -1 implies that no column exists in Z[p-1] with pivot i.
    // Data structures for storing the intervals and cycles in order to compute representatives: 
    vector<vector<column> > bundle; // A collection of wires: bundle[p] stores p-dimensional wires.
    vector<vector<int> > timestamp; // Map from the cycle matrix to the wires: links[p][j] is the collection of wires associated with Z[p][j].
    vector<vector<column> > links; // Represents a map from the cycle matrix to the wires: links[p][j] is the collection of wires associated with Z[p][j]; we use the same used/available columns as in Z.
    /*
    INITIALIZATION:
    */
    for (int d = 0; d <= m; ++d)
    {
        SimplexIdMap id_p; 
        vector<int> unique_id_p;
        vector<column> Z_p;
        vector<int> available_columns_Z_p;
        vector<column> C_p; 
        vector<int> used_columns_Z_p;
        vector<pair <int, int>> birth_timestamp_p;
        PivotMap pivots_p;
        id.push_back(id_p);
        unique_id.push_back(unique_id_p);
        Z.push_back(Z_p);
        available_columns_Z.push_back(available_columns_Z_p);
        used_columns_Z.push_back(used_columns_Z_p);
        C.push_back(C_p);
        available_columns_C.push_back(available_columns_Z_p);
        used_columns_C.push_back(used_columns_Z_p);
        birth_timestamp.push_back(birth_timestamp_p);
        pivots.push_back(pivots_p);
        vector<column> bundle_p;
        vector<int> timestamp_p;
        vector<column> links_p;
        bundle.push_back(bundle_p);
        timestamp.push_back(timestamp_p);
        links.push_back(links_p);
    }
    /*
    COMPUTATION of zigzag persistence and the associated representatives:
    */
    for (int i = 0; i < n; ++i) {
        const vector<int> &simp = filt_simp[i];
        int p = simp.size() - 1; // p denotes the dimension of the simplex.
        if (filt_op[i]) { // INSERTION
            // Find the unique id of the simplex.
            int unique_id_simp;
            if (id[p].find(simp) != id[p].end()) // If key does exist already:
            {
                unique_id_simp = unique_id[p][id[p].at(simp)];
            }
            else {
                unique_id_simp = unique_id[p].size();
            }
            // Change the id of the simplex to the new id.
            id[p][simp] = unique_id[p].size();
            unique_id[p].push_back(unique_id_simp);
            // Represent the boundary of simp as a sum of columns of Z_{p-1} by a reduction algorithm; I is such set of columns.
            bool all_boundary = true;
            column bd_simp;
            vector<int> I;
            if (p != 0) {
                // Compute boundary of the simplex:
                bd_simp = make_shared<bitset>(unique_id[p-1].size(), 0);
                for (int i = 0; i <= p; i++) {
                    // Remove the i-th vertex.
                    vector<int> boundary_simplex = simp;
                    boundary_simplex.erase(boundary_simplex.begin() + i);
                    int idx = id[p-1].at(boundary_simplex);
                    bd_simp -> set(idx);
                }
                /* 
                Find the columns in Z[p-1] that sum to bd_simp:
                */
                // Copy bd_simp to a temporary column.
                column bd_simp_temp;
                bd_simp_temp = column(new bitset ( *bd_simp)); 
                int pivot_bd = pivot(bd_simp_temp);
                bool zeroed = (pivot_bd == -1);
                while (!zeroed) {
                    // Find the column in M that has the same pivot as a.
                    // assert(pivots[p-1].find(pivot_bd) != pivots[p-1].end());
                    int conflict_idx = pivots[p-1].at(pivot_bd);
                    // Add this column to indices.
                    dynamic_xor(bd_simp_temp, Z[p-1][conflict_idx]);
                    I.push_back(conflict_idx);
                    // Update the pivot of a.
                    pivot_bd = pivot(bd_simp_temp);
                    zeroed = (pivot_bd == -1);
                } 
                // Check the birth timestamps to check whether all of them are boundaries.
                for (auto a: I) {
                    if (birth_timestamp[p-1][a].first >= 0) {
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
                // New column with of size unique_id[p].size() with 1 appended at the end.
                column new_column; 
                new_column = make_shared<bitset>(unique_id[p].size(), 0);
                new_column -> set(id[p].at(simp));
                if (p != 0)  {
                    for (auto a: I) {
                        int chain_a = birth_timestamp[p-1][a].second; // Since a is a boundary, we know that the birth timestamp is non-negative and this is safe!
                        dynamic_xor(new_column, C[p][chain_a]);
                    }
                }
                uint pivot_new = pivot(new_column);
                // Add a wire to the bundle with timestamp i+1 and update the link from the cycle to the bundle.
                bundle[p].push_back(new_column);
                timestamp[p].push_back(i+1);
                column new_link;
                new_link = make_shared<bitset>(bundle[p].size(), 0);
                new_link -> set(bundle[p].size()-1);
                // Now, we add the new column to the cycle matrix and the new link to the link matrix.
                if (available_columns_Z[p].size() != 0) 
                {
                    uint av_idx = available_columns_Z[p].back();
                    Z[p][av_idx] = new_column;
                    links[p][av_idx] = new_link;    
                    birth_timestamp[p][av_idx] = make_pair(i+1, -1);
                    pivots[p][pivot_new] = av_idx;
                    // Update the used and available columns:
                    used_columns_Z[p].push_back(av_idx);
                    available_columns_Z[p].pop_back();
                }
                else 
                {
                    Z[p].push_back(new_column);
                    links[p].push_back(new_link);
                    birth_timestamp[p].push_back(make_pair(i+1, -1));
                    pivots[p][pivot_new] = Z[p].size() - 1;
                    // Update the used and available columns:
                    used_columns_Z[p].push_back(Z[p].size() - 1);
                }
            }
            else { // The inserted simplex kills a (p-1)-cycle (an interval in dimension (p-1) dies).
                /*
                Let J consist of indices in I whose corresponding columns in Z[p−1] have non-negative birth timestamps. 
                */ 
                vector<int> J;
                for (int a: I) {
                    if (birth_timestamp[p-1][a].first >= 0) { // Gather all the non-boundary cycles.
                        J.push_back(a);
                    }
                }
                sort(J.begin(), J.end());
                // Check if arrow at b^{p−1}[c]−1 points backward for all c in J
                bool arrow_backward = true;
                for (auto c: J) {
                    if (filt_op[birth_timestamp[p-1][c].first - 1]) {
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
                        if (filt_op[birth_timestamp[p-1][J[i]].first - 1]) {
                            l = J[i];
                            break;
                        }
                    }
                }
                /* 
                OUTPUT the (p − 1)-th interval [b^{p−1}[l], i] after gathering the relevant representative.
                */
                column wire_column = links[p-1][l];
                // Here, we need to gather all the wires pertaining to the cycle l and reconstruct the representatives at each index.
                vector<tuple<int, column>> wire_representatives;
                vector<tuple<int, vector<int>>> wire_representatives_ids;
                for (int col_idx = 0; col_idx < wire_column -> size(); ++col_idx) {
                    if ((*wire_column)[col_idx] == 1) {
                        wire_representatives.push_back(make_tuple(timestamp[p-1][col_idx], bundle[p-1][col_idx]));
                    }
                } 
                // Sort the representatives in the order of the timestamps.
                sort(wire_representatives.begin(), wire_representatives.end(), [&](tuple<int, column> &a, tuple<int, column> &b){ return (get<0>(a) < get<0>(b));});
                // Initialize a bitset column of size representative_max_size with all zeros.
                column current_representative = make_shared<bitset>(1, 0);
                uint interval_birth = birth_timestamp[p-1][l].first;
                // Add all the representatives with timestamps less than or equal to interval_birth.
                uint rep_idx;
                for (rep_idx = 0; rep_idx < wire_representatives.size(); ++rep_idx) {
                    if (get<0>(wire_representatives[rep_idx]) <= interval_birth) {
                        dynamic_xor(current_representative, get<1>(wire_representatives[rep_idx]));
                    }
                }
                vector<int> indices;
                for (int rep_idx = 0; rep_idx < current_representative -> size(); ++rep_idx) {
                    if ((*current_representative)[rep_idx] == 1) {
                        indices.push_back(unique_id[p-1][rep_idx]);
                    }
                }
                wire_representatives_ids.push_back(make_tuple(interval_birth, indices));
                for (; rep_idx < wire_representatives.size(); ++rep_idx) {
                    int time = get<0>(wire_representatives[rep_idx]);
                    vector<int> indices;
                    dynamic_xor(current_representative, get<1>(wire_representatives[rep_idx]));
                    for (int rep_idx = 0; rep_idx < current_representative -> size(); ++rep_idx) {
                        if ((*current_representative)[rep_idx] == 1) {
                            indices.push_back(unique_id[p-1][rep_idx]);
                        }
                    }
                    wire_representatives_ids.push_back(make_tuple(time, indices));
                }
                persistence -> push_back(make_tuple(interval_birth, i + 1, p - 1, wire_representatives_ids));
                /*
                Bundles and Links UPDATE:
                */
                // Add a new wire to the (p-1)-th bundle:
                bundle[p-1].push_back(bd_simp);
                timestamp[p-1].push_back(i+1);
                column new_link = make_shared<bitset>(bundle[p-1].size(), 0);
                new_link -> set(bundle[p-1].size()-1);
                links[p-1][l] = new_link;
                /* 
                UPDATE:
                */
                // Set Z[p−1][l] = bd_simp.
                uint prev_pivot = pivot(Z[p-1][l]);
                Z[p-1][l] = bd_simp;
                // Set C[p][l] = simp and update the boundary-to-chain map.
                column current_chain;
                current_chain = make_shared<bitset>(unique_id[p].size(), 0);
                current_chain -> set(id[p].at(simp)); // Set last entry to 1.
                if (available_columns_C[p].size() != 0) 
                {
                    uint av_idx = available_columns_C[p].back();
                    C[p][av_idx] = current_chain;
                    birth_timestamp[p-1][l] = make_tuple(-1, av_idx);
                    // Update the used and available columns:
                    used_columns_C[p].push_back(av_idx);
                    available_columns_C[p].pop_back();
                }
                else 
                {
                    C[p].push_back(current_chain);
                    uint av_idx = C[p].size() - 1;
                    birth_timestamp[p-1][l] = make_tuple(-1, av_idx);
                    // Update the used and available columns:
                    used_columns_C[p].push_back(av_idx);
                }
                /*
                Avoiding PIVOT CONFLICTS: (column of bd_simp with that of another column in Z[p-1]). We know that the pivots of Z[p-1] were unique before the insertion of bd_simp. Thus, we need to ensure that they remain unique after the insertion of bd_simp. 
                */
                bool pivot_conflict;
                int a, b, pivot_a, pivot_b, current_idx, current_pivot;
                pivots[p-1].erase(prev_pivot); // Remove the previous pivot
                current_pivot = pivot(Z[p-1][l]);
                if (pivots[p-1].find(current_pivot) == pivots[p-1].end())
                {   
                    pivot_conflict = false;
                    pivots[p-1][current_pivot] = l;
                }
                else {
                    current_idx = pivots[p-1][current_pivot];
                    pivot_conflict = true;
                    a = l;
                    b = current_idx;
                    pivot_a = current_pivot;
                    pivot_b = current_pivot;
                }
                while (pivot_conflict)
                {
                    // Check that a and b are valid indices in Z[p-1]:
                    if (birth_timestamp[p-1][a].first == -1 && birth_timestamp[p-1][b].first == -1) {
                        dynamic_xor(Z[p-1][a], Z[p-1][b]);
                        dynamic_xor(links[p-1][a], links[p-1][b]);
                        // Update chain matrices to reflect this addition (use the BdToChainMap):
                        dynamic_xor(C[p][birth_timestamp[p-1][a].second], C[p][birth_timestamp[p-1][b].second]);
                        // Check for pivot conflicts again (pivot of b remains the same, but the pivot of a might have changed):
                        pivots[p-1][pivot_b] = b;
                        pivot_a = pivot(Z[p-1][a]);
                        if (pivots[p-1].find(pivot_a) == pivots[p-1].end())
                        {
                            pivot_conflict = false;
                            pivots[p-1][pivot_a] = a;
                        }
                        else 
                        {
                            if (pivots[p-1][pivot_a] == a) {
                                pivot_conflict = false;
                            }
                            else 
                            {
                                current_idx = pivots[p-1][pivot_a];
                                pivot_conflict = true;
                                b = current_idx;
                                pivot_b = pivot_a;
                            }
                        }

                    }
                    else if (birth_timestamp[p-1][a].first == -1 && birth_timestamp[p-1][b].first >= 0) {
                        dynamic_xor(Z[p-1][b], Z[p-1][a]);
                        dynamic_xor(links[p-1][b], links[p-1][a]);
                        // No need to update the chain matrices as a boundary is added to a cycle. Instead, check for pivot conflicts again:
                        pivots[p-1][pivot_a] = a;
                        pivot_b = pivot(Z[p-1][b]);
                        if (pivots[p-1].find(pivot_b) == pivots[p-1].end())
                        {
                            pivot_conflict = false;
                            pivots[p-1][pivot_b] = b;
                        }
                        else 
                        {
                            if (pivots[p-1][pivot_b] == b) {
                                pivot_conflict = false;
                            }
                            else 
                            {
                                current_idx = pivots[p-1][pivot_b];
                                pivot_conflict = true;
                                a = current_idx;
                                pivot_a = pivot_b;
                            }
                        }
                    }
                    else if (birth_timestamp[p-1][a].first >= 0 && birth_timestamp[p-1][b].first == -1) {
                        dynamic_xor(Z[p-1][a], Z[p-1][b]);
                        dynamic_xor(links[p-1][a], links[p-1][b]);
                        // Again, no need to update the chain matrices as a boundary is added to a cycle. Instead, check for pivot conflicts again:
                        pivots[p-1][pivot_b] = b;
                        pivot_a = pivot(Z[p-1][a]);
                        if (pivots[p-1].find(pivot_a) == pivots[p-1].end())
                        {
                            pivot_conflict = false;
                            pivots[p-1][pivot_a] = a;
                        }
                        else {
                            if (pivots[p-1][pivot_a] == a) {
                                pivot_conflict = false;
                            }
                            else 
                            {
                                current_idx = pivots[p-1][pivot_a];
                                pivot_conflict = true;
                                b = current_idx;
                                pivot_b = pivot_a;
                            }
                        }
                    }
                    else {
                        int birth_a = birth_timestamp[p-1][a].first;
                        int birth_b = birth_timestamp[p-1][b].first;
                        bool a_lessthan_b = ((birth_a == birth_b) || ((birth_a < birth_b) && (filt_op[birth_b-1])) || ((birth_a > birth_b) && (!(filt_op[birth_a-1]))));   
                        if (a_lessthan_b) {
                            dynamic_xor(Z[p-1][b], Z[p-1][a]);
                            dynamic_xor(links[p-1][b], links[p-1][a]);
                            // Since both are cycles, no need to update the chain matrices. Instead, check for pivot conflicts again:
                            pivots[p-1][pivot_a] = a;
                            pivot_b = pivot(Z[p-1][b]);
                            if (pivots[p-1].find(pivot_b) == pivots[p-1].end())
                            {
                                pivot_conflict = false;
                                pivots[p-1][pivot_b] = b;
                            }
                            else 
                            {
                                if (pivots[p-1][pivot_b] == b) {
                                    pivot_conflict = false;
                                }
                                else {
                                    current_idx = pivots[p-1][pivot_b];
                                    pivot_conflict = true;
                                    a = current_idx;
                                    pivot_a = pivot_b;
                                }
                            }
                        }
                        else {
                            dynamic_xor(Z[p-1][a], Z[p-1][b]);
                            dynamic_xor(links[p-1][a], links[p-1][b]);
                            // Since both are cycles, no need to update the chain matrices. Instead, check for pivot conflicts again:
                            pivots[p-1][pivot_b] = b;
                            pivot_a = pivot(Z[p-1][a]);
                            if (pivots[p-1].find(pivot_a) == pivots[p-1].end())
                            {
                                pivot_conflict = false;
                                pivots[p-1][pivot_a] = a;
                            }
                            else 
                            {
                                if (pivots[p-1][pivot_a] == a) {
                                    pivot_conflict = false;
                                }
                                else {
                                    current_idx = pivots[p-1][pivot_a];
                                    pivot_conflict = true;
                                    b = current_idx;
                                    pivot_b = pivot_a;
                                }
                            }
                        }
                    }
                }
            }
        } 
        else { // DELETION:
            // Find the index of the simplex in id[p].
            int idx = id[p].at(simp);
            int exist_col_idx = -1;
            // Check if the simplex exists in Z[p].
            for (auto col_idx: used_columns_Z[p]) 
            {
                if (idx < Z[p][col_idx] -> size() && (*Z[p][col_idx])[idx] == 1) {
                    exist_col_idx = col_idx;
                    break;
                }
            }
            if (exist_col_idx == -1)// The deleted simplex does not constitute a cycle and its boundary now becomes a (p-1)-cycle (An interval in dimension (p-1) gets born).
            {   
                // Find all the columns in C[p] that contain the simplex and get the indices of their corresponding boundaries.
                vector<int> I;
                int alpha, chain_alpha; // alpha will be the only column with negative birth timestamp such that C[p][chain_a] contains the simplex.
                int smallest_pivot = unique_id[p-1].size();
                for (auto cyc_idx: used_columns_Z[p-1])
                {
                    if (birth_timestamp[p-1][cyc_idx].first == -1) 
                    {
                        int chain_idx = birth_timestamp[p-1][cyc_idx].second;
                        if (idx < C[p][chain_idx] -> size() && (*C[p][chain_idx])[idx]==1) {
                            I.push_back(cyc_idx);
                            int pivot_cyc_idx = pivot(Z[p-1][cyc_idx]);
                            if (pivot_cyc_idx < smallest_pivot) {
                                smallest_pivot = pivot_cyc_idx;
                                alpha = cyc_idx;
                                chain_alpha = chain_idx;
                            }
                        }
                    }
                }
                // Add alpha to all the other columns in I.
                for (auto cyc_idx: I) {
                    int chain_idx = birth_timestamp[p-1][cyc_idx].second;
                    if (cyc_idx != alpha) {
                        dynamic_xor(Z[p-1][cyc_idx], Z[p-1][alpha]);
                        dynamic_xor(links[p-1][cyc_idx], links[p-1][alpha]);
                        dynamic_xor(C[p][chain_idx], C[p][chain_alpha]);
                    }
                }
                C[p][chain_alpha] = nullptr; // Zero out the chain C[p][chain_alpha] from C[p].
                // Update the used and available indices for C[p]:
                available_columns_C[p].push_back(chain_alpha);
                used_columns_C[p].erase(remove(used_columns_C[p].begin(), used_columns_C[p].end(), chain_alpha), used_columns_C[p].end());
                birth_timestamp[p-1][alpha] = make_pair(i+1, -1);
                // Add a new wire to the (p-1)-th bundle:
                bundle[p-1].push_back(Z[p-1][alpha]);
                timestamp[p-1].push_back(i+1);
                column new_link = make_shared<bitset>(bundle[p-1].size(), 0);
                new_link -> set(bundle[p-1].size()-1);
                links[p-1][alpha] = new_link;
            }
            else // The deleted simplex is part of some p-cycle (An interval in dimension p dies).
            {
                // Update C[p] so that no columns contain the simplex.
                for (auto chain_idx : used_columns_C[p]) 
                {
                    if (idx < C[p][chain_idx] -> size() && (*C[p][chain_idx])[idx] == 1) // each column in C[p][chain_b] containing the simplex s.t. Z[p-1][b] is a boundary.
                    {  
                        // Here, we do not need to check whether Z[p-1][b] is a boundary since we only store the chains whose boundaries are non-zero.
                        dynamic_xor(C[p][chain_idx], Z[p][exist_col_idx]);
                    }
                }
                // Remove the simplex from Z[p]
                vector<int> I; // Gather indices of columns that contain the simplex. 
                for (auto col_idx: used_columns_Z[p]) 
                {
                    if (idx < (Z[p][col_idx] -> size()) && (*Z[p][col_idx])[idx] == 1)
                    {
                        I.push_back(col_idx);
                    }
                }
                // sort I in the order of the birth timestamps where the order is the total order as above without using the sort function.
                sort(I.begin(), I.end(), [&](int &a, int &b){
                    int birth_a = birth_timestamp[p][a].first;
                    int birth_b = birth_timestamp[p][b].first;
                     return ((birth_a == birth_b) || ((birth_a < birth_b) && (filt_op[birth_b-1])) || ((birth_a > birth_b) && (!(filt_op[birth_a-1]))));
                     });
                // The column to be deleted is the first column in I.
                column z = Z[p][I[0]];
                int alpha = I[0];
                int current_alpha = I[0]; // Keeps track of the column index being added in order to update links.
                I.erase(I.begin());
                int alpha_pivot = pivot(z);
                pivots[p].erase(alpha_pivot); // Remove the pivot of Z[p][alpha] from pivots[p].
                if (I.size() != 0) 
                {
                    // Iterate over the rest of I and ensure that the pivots remain distinct:
                    int z_pivot = alpha_pivot;
                    for (auto a: I)
                    {
                        int a_pivot = pivot(Z[p][a]);
                        if (a_pivot > z_pivot) 
                        {
                            dynamic_xor(Z[p][a], z);
                            dynamic_xor(links[p][a], links[p][current_alpha]);
                        }   
                        else 
                        {
                            // We want to deep copy Z[p][a] to a temporary column and then copy z to Z[p][a].
                            column temp = make_shared<bitset>(*Z[p][a]);
                            int temp_pivot = a_pivot;
                            dynamic_xor(Z[p][a], z);
                            dynamic_xor(links[p][a], links[p][current_alpha]);
                            a_pivot = z_pivot;
                            // Update the information for z's.
                            z = temp;
                            current_alpha = a;
                            z_pivot = temp_pivot;
                            // Update the pivot of Z[p][a] in pivots[p].
                            pivots[p][a_pivot] = a;
                            pivots[p].erase(z_pivot);
                        }
                    }
                }
                /*
                OUTPUT the p-th interval [birth_timestamp_p[alpha], i] after gathering the relevant representative.
                */
                column wire_column = links[p][alpha];
                // Here, we need to gather all the wires pertaining to the cycle l and reconstruct the representatives at each index.
                vector<tuple<int, column>> wire_representatives;
                vector<tuple<int, vector<int>>> wire_representatives_ids;
                for (int col_idx = 0; col_idx < wire_column -> size(); ++col_idx) {
                    if ((*wire_column)[col_idx] == 1) {
                        wire_representatives.push_back(make_tuple(timestamp[p][col_idx], bundle[p][col_idx]));
                    }
                } 
                // Sort the representatives in the order of the timestamps.
                sort(wire_representatives.begin(), wire_representatives.end(), [&](tuple<int, column> &a, tuple<int, column> &b){ return (get<0>(a) < get<0>(b));});
                // Initialize a bitset column of size representative_max_size with all zeros.
                column current_representative = make_shared<bitset>(1, 0);
                uint interval_birth = birth_timestamp[p][alpha].first;
                // Add all the representatives with timestamps less than or equal to interval_birth.
                uint rep_idx;
                for (rep_idx = 0; rep_idx < wire_representatives.size(); ++rep_idx) {
                    if (get<0>(wire_representatives[rep_idx]) <= interval_birth) {
                        dynamic_xor(current_representative, get<1>(wire_representatives[rep_idx]));
                    }
                }
                vector<int> indices;
                for (int rep_idx = 0; rep_idx < current_representative -> size(); ++rep_idx) {
                    if ((*current_representative)[rep_idx] == 1) {
                        indices.push_back(unique_id[p][rep_idx]);
                    }
                }
                wire_representatives_ids.push_back(make_tuple(interval_birth, indices));
                for (; rep_idx < wire_representatives.size(); ++rep_idx) {
                    int time = get<0>(wire_representatives[rep_idx]);
                    vector<int> indices;
                    dynamic_xor(current_representative, get<1>(wire_representatives[rep_idx]));
                    for (int rep_idx = 0; rep_idx < current_representative -> size(); ++rep_idx) {
                        if ((*current_representative)[rep_idx] == 1) {
                            indices.push_back(unique_id[p][rep_idx]);
                        }
                    }
                    wire_representatives_ids.push_back(make_tuple(time, indices));
                }
                persistence -> push_back(make_tuple(interval_birth, i + 1, p, wire_representatives_ids));
                /*
                UPDATE:
                */ 
                // Remove the columns Z[p][alpha] from Z[p] and link[p][alpha] from links: assign these to null columns.
                Z[p][alpha] = nullptr;
                links[p][alpha] = nullptr;
                // Add this to the available columns for adding new cycles:
                available_columns_Z[p].push_back(alpha);
                used_columns_Z[p].erase(remove(used_columns_Z[p].begin(), used_columns_Z[p].end(), alpha), used_columns_Z[p].end());
                // We need to assign this to be invalid (redundant):
                birth_timestamp[p][alpha] = make_pair(-2, -1);
            }
        }
    }
    /*
     POST-PROCESSING:
    */
    // For each p and each column Z[p][a] of Z[p] with non-negative birth timestamp, output the p-th interval [birth_timestamp[p][a], n].
    for (int p = 0; p <= m; p++) {
        for (auto a: used_columns_Z[p]) {
            uint interval_birth = birth_timestamp[p][a].first;
            if (interval_birth >= 0) {
                column wire_column = links[p][a];
                // Here, we need to gather all the wires pertaining to the cycle l and reconstruct the representatives at each index.
                vector<tuple<int, column>> wire_representatives;
                vector<tuple<int, vector<int>>> wire_representatives_ids;
                for (int col_idx = 0; col_idx < wire_column -> size(); ++col_idx) {
                    if ((*wire_column)[col_idx] == 1) {
                        wire_representatives.push_back(make_tuple(timestamp[p][col_idx], bundle[p][col_idx]));
                    }
                } 
                // Sort the representatives in the order of the timestamps.
                sort(wire_representatives.begin(), wire_representatives.end(), [&](tuple<int, column> &a, tuple<int, column> &b){ return (get<0>(a) < get<0>(b));});
                // Initialize a bitset column of size representative_max_size with all zeros.
                column current_representative = make_shared<bitset>(1, 0);
                // Add all the representatives with timestamps less than or equal to interval_birth.
                uint rep_idx;
                for (rep_idx = 0; rep_idx < wire_representatives.size(); ++rep_idx) {
                    if (get<0>(wire_representatives[rep_idx]) <= interval_birth) {
                        dynamic_xor(current_representative, get<1>(wire_representatives[rep_idx]));
                    }
                }
                vector<int> indices;
                for (int rep_idx = 0; rep_idx < current_representative -> size(); ++rep_idx) {
                    if ((*current_representative)[rep_idx] == 1) {
                        indices.push_back(unique_id[p][rep_idx]);
                    }
                }
                wire_representatives_ids.push_back(make_tuple(interval_birth, indices));
                for (; rep_idx < wire_representatives.size(); ++rep_idx) {
                    int time = get<0>(wire_representatives[rep_idx]);
                    vector<int> indices;
                    dynamic_xor(current_representative, get<1>(wire_representatives[rep_idx]));
                    for (int rep_idx = 0; rep_idx < current_representative -> size(); ++rep_idx) {
                        if ((*current_representative)[rep_idx] == 1) {
                            indices.push_back(unique_id[p][rep_idx]);
                        }
                    }
                    wire_representatives_ids.push_back(make_tuple(time, indices));
                }
                persistence -> push_back(make_tuple(interval_birth, n, p, wire_representatives_ids));
            }
        }
    }
}

/*
IMPLEMENTATION OF HELPER FUNCTIONS:
*/
// Goes over a vector of ints and finds the position of the last non-zero element.
int pivot (column a)
{
    int pivot = (a -> find_first() != a -> npos) ? (a -> find_first()) : -1;
    // Find the next non-zero element.
    while (a -> find_next(pivot) != a -> npos) 
    {
        pivot = a -> find_next(pivot);
    }
    return pivot;
}

// Takes the XOR of bitsets of differing lengths. Instead of appending by 0s and then XORing, we XOR the bitsets directly.
/*
We assume that we compute a ^= b.
*/
void dynamic_xor(column a, column b)
{
    size_t a_size = a -> size();
    size_t b_size = b -> size();
    if (a_size < b_size) {
        a -> resize(b_size, 0);
    }
    for (int i = b -> find_first(); i != b -> npos; i = b -> find_next(i)) {
        (*a)[i] ^= (*b)[i];
    }
}
}// namespace ZZREP