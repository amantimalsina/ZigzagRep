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
    // If non_bd is true, then val_idx records the birth timestamp of the non-boundary cycle,
    // otherwise, val_idx records to the index of the corresponding chain
    bool non_bd;
    int val_idx;
    birth_timestamp(bool non_bd, int val_idx) : non_bd(non_bd), val_idx(val_idx) {} 
};

// Bitsets:
typedef boost::dynamic_bitset<> bitset;
typedef shared_ptr<bitset> pbitset;

// Maps and Pairs:
typedef std::map< vector<int>, int> SimplexIdMap;
typedef std::map<int, int> PivotMap;
typedef SimplexIdMap::value_type SimplexIdPair;
typedef PivotMap::value_type PivotPair;


/* DECLARATION OF HELPER FUNCTIONS: */
int pivot(pbitset a);

void dynamic_xor(pbitset a, pbitset b);

void boundary_as_cycles();

void add_new_cycle();

void add_new_chain();

void make_pivots_distinct(
    const std::vector<bool> *filt_op, 
    vector<pbitset> *Z_p,
    vector<pbitset> *links_p,
    vector<pbitset> *C_p,
    vector<pair <bool, int> > *birth_timestamp_p,
    PivotMap *pivots_p,
    const int prev_pivot,
    const int l
    );

void reduce_boundaries(); // TODO: Implement this as well.

void delete_cycle(); // TODO: Implement this as well.

void output_representatives(
    const int p, 
    const int birth, 
    const int death,
    pbitset link_interval,
    vector<int> *unique_id_p,
    vector<pbitset> *bundle,
    vector<int> *timestamp,
    std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > *persistence
    );

void ZigzagRep::compute(
        const std::vector<vector<int> > &filt_simp, 
        const std::vector<bool> &filt_op,
        std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > *persistence, 
        std::vector <std::map<int, int>> *i_to_id,
        const int m) {
    persistence -> clear();
    const int n = filt_simp.size();
    // DECLARATION & INITIALIZATION of Data Structures:
    vector<SimplexIdMap> id(m+1, SimplexIdMap()); //  The current integral id for each p-simplex.
    vector<vector<int>> unique_id(m+1, vector<int>()); // The unique integral id for each p-simplex.
    vector<vector<pbitset> > Z(m+1, vector<pbitset>()); // A cycle matrix Z[p] for each dimension p that is maintained such that each pbitset of Z[p] represents a p-cycle in K_i.
    vector<vector<int>> available_pbitsets_Z(m+1, vector<int>()); // A list of available pbitsets in Z[p] for each dimension p.
    vector<vector<int>> used_pbitsets_Z(m+1, vector<int>()); // A list of used pbitsets in Z[p] for each dimension p.
    vector<vector<pbitset> > C(m+1, vector<pbitset>()); // A chain matrix C[p] for each dimension p so that the chains in C[p] are associated with the the boundary cycles in Z[p-1].
    vector<vector<int>> available_pbitsets_C(m+1, vector<int>()); // A list of available pbitsets in Z[p] for each dimension p.
    vector<vector<int>> used_pbitsets_C(m+1, vector<int>()); // A list of used pbitsets in Z[p] for each dimension p.
    vector<vector<pair <bool, int> > > birth_timestamp(m+1, vector<pair <bool, int> >()); // For each pbitset Z[p][j], a birth timestamp is maintained, where a false value implies that the pbitset pertains to a boundary.
    vector<PivotMap> pivots(m+1, PivotMap()); // For each pbitset Z[p][j] we have a unique pivot entry; we then store pivots[p][i] = j where j = -1 implies that no pbitset exists in Z[p-1] with pivot i.
    // Data structures for storing the intervals and cycles in order to compute representatives: 
    vector<vector<pbitset> > bundle(m+1, vector<pbitset>()); // A collection of wires: bundle[p] stores p-dimensional wires.
    vector<vector<int> > timestamp(m+1, vector<int>()); // Map from the cycle matrix to the wires: links[p][j] is the collection of wires associated with Z[p][j].
    vector<vector<pbitset> > links(m+1, vector<pbitset>()); // Represents a map from the cycle matrix to the wires: links[p][j] is the collection of wires associated with Z[p][j]; we use the same used/available pbitsets as in Z.
    // COMPUTATION of zigzag persistence and the associated representatives:
    for (int i = 0; i < n; ++i) {
        const vector<int> &simp = filt_simp[i];
        const int p = simp.size() - 1; // p denotes the dimension of the simplex.
        if (filt_op[i]) { // INSERTION
            int unique_id_simp; // Find the unique id of the simplex.
            id[p].find(simp) != id[p].end() ? unique_id_simp = unique_id[p][id[p].at(simp)] : unique_id_simp = unique_id[p].size();
            id[p][simp] = unique_id[p].size(); // Change the id of the simplex to the new id.
            unique_id[p].push_back(unique_id_simp);
            (*i_to_id)[p][i] = unique_id_simp;
            // Represent the boundary of simp as a sum of pbitsets of Z_{p-1} by a reduction algorithm; I is such set of pbitsets.
            bool all_boundary = true;
            pbitset bd_simp;
            vector<int> I;
            if (p != 0) boundary_as_cycles();
            if (all_boundary) { // FORWARD BIRTH.
                add_new_cycle(); 
            }
            else { // The inserted simplex kills a (p-1)-cycle (an interval in dimension (p-1) dies).
                vector<int> J; // Let J consist of indices in I whose corresponding pbitsets in Z[p−1] have non-negative birth timestamps. 
                for (int a: I) if (birth_timestamp[p-1][a].first) J.push_back(a);
                sort(J.begin(), J.end());
                // Check if arrow at b^{p−1}[c]−1 points backward for all c in J
                int k;
                for (k = J.size()-1; k >= 0; k--) if (filt_op[birth_timestamp[p-1][J[k]].second - 1]) break;
                const int l = J[l]; // If the arrow at {b^{p−1}[c]−1} points backward for all c in J, let l be the smallest index in J. 
                /* OUTPUT the (p − 1)-th interval [b^{p−1}[l], i] after gathering the relevant representative. */
                output_representatives(p-1, birth_timestamp[p-1][l].second, i, links[p-1][l], &(unique_id[p-1]), &(bundle[p-1]), &(timestamp[p-1]), persistence);
                // Add new chain
                add_new_chain();
                /* Avoiding PIVOT CONFLICTS: (pbitset of bd_simp with that of another pbitset in Z[p-1]). We know that the pivots of Z[p-1] were unique before the insertion of bd_simp. Thus, we need to ensure that they remain unique after the insertion of bd_simp. */
                make_pivots_distinct(&filt_op, &Z[p-1], &links[p-1], &C[p-1], &birth_timestamp[p-1], &pivots[p-1], prev_pivot, l);
            }
        } 
        else { // DELETION:
            // Find the index of the simplex in id[p].
            int idx = id[p].at(simp);
            int used_idx;
            for (used_idx = 0; used_idx < used_pbitsets_Z[p].size(); --used_idx) if (idx < Z[p][used_idx] -> size() && (*Z[p][used_idx])[idx] == 1) break;
            if (used_idx == used_pbitsets_Z[p].size()) { // BACKWARD BIRTH.  
                int alpha, chain_alpha;
                reduce_boundaries();
                C[p][chain_alpha] = nullptr;
                available_pbitsets_C[p].push_back(chain_alpha); 
                used_pbitsets_C[p].erase(remove(used_pbitsets_C[p].begin(), used_pbitsets_C[p].end(), chain_alpha), used_pbitsets_C[p].end()); // TODO: Check that this is done correctly.
                birth_timestamp[p-1][alpha] = make_pair(true, i+1);
                // Update Bundles and Links:
                bundle[p-1].push_back(Z[p-1][alpha]);
                timestamp[p-1].push_back(i+1);
                pbitset new_link = make_shared<bitset>(bundle[p-1].size(), 0);
                new_link -> set(bundle[p-1].size()-1);
                links[p-1][alpha] = new_link;
            }
            else { // BACKWARD DEATH,
                int simp_col_idx = used_pbitsets_Z[p][used_idx];
                for (auto chain_idx : used_pbitsets_C[p]) if (idx < C[p][chain_idx] -> size() && (*C[p][chain_idx])[idx] == 1) dynamic_xor(C[p][chain_idx], Z[p][simp_col_idx]);  // Update C[p] so that no pbitsets contain the simplex.
                vector<int> I; // Remove the simplex from Z[p]
                for (auto col_idx: used_pbitsets_Z[p])  if (idx < (Z[p][col_idx] -> size()) && (*Z[p][col_idx])[idx] == 1) I.push_back(col_idx); // Gather indices of pbitsets that contain the simplex. 
                // sort I in the order of the birth timestamps where the order is the total order as above without using the sort function.
                sort(I.begin(), I.end(), [&](int &a, int &b){ int birth_a = birth_timestamp[p][a].second; int birth_b = birth_timestamp[p][b].second; return ((birth_a == birth_b) || ((birth_a < birth_b) && (filt_op[birth_b-1])) || ((birth_a > birth_b) && (!(filt_op[birth_a-1]))));});
                const int alpha = I[0]; // The pbitset to be deleted is the first pbitset in I.
                delete_cycle();
                /* OUTPUT the p-th interval [birth_timestamp_p[alpha], i] after gathering the relevant representative. */
                uint birth = birth_timestamp[p][alpha].second;
                output_representatives(p, birth, i, links[p][alpha], &(unique_id[p]), &(bundle[p]), &(timestamp[p]), persistence);
                /* UPDATE: Remove the pbitsets Z[p][alpha] from Z[p] and link[p][alpha] from links: assign these to null pbitsets. */ 
                Z[p][alpha] = nullptr;
                links[p][alpha] = nullptr;
                available_pbitsets_Z[p].push_back(alpha);       // Add this to the available pbitsets for adding new cycles:
                used_pbitsets_Z[p].erase(remove(used_pbitsets_Z[p].begin(), used_pbitsets_Z[p].end(), alpha), used_pbitsets_Z[p].end());
                birth_timestamp[p][alpha] = make_pair(false, -1); // We need to assign this to be invalid (redundant). TODO: Probably use assertions to test that boudaries are not using invalid pbitsets.
            }
        }
    }
    /*
     POST-PROCESSING:
    */
    // For each p and each pbitset Z[p][a] of Z[p] with non-negative birth timestamp, output the p-th interval [birth_timestamp[p][a], n].
    for (size_t p = 0; p <= m; p++) {
        for (auto a: used_pbitsets_Z[p]) {
            if (birth_timestamp[p][a].first) {
                uint birth = birth_timestamp[p][a].second;
                output_representatives(p, birth, n, links[p][a], &(unique_id[p]), &(bundle[p]), &(timestamp[p]), persistence);
            }
        }
    }
}

/*
IMPLEMENTATION OF HELPER FUNCTIONS:
*/
// Goes over a vector of ints and finds the position of the last non-zero element.
int pivot (pbitset a)
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
void dynamic_xor(pbitset a, pbitset b)
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

void boundary_as_cycles(){
    for (int i = 0; i <= p; i++) {
        // Remove the i-th vertex.
        vector<int> boundary_simplex = simp;
        boundary_simplex.erase(boundary_simplex.begin() + i);
        int idx = id[p-1].at(boundary_simplex);
        bd_simp -> set(idx);
    }
    /* 
    Find the pbitsets in Z[p-1] that sum to bd_simp:
    */
    // Copy bd_simp to a temporary pbitset.
    pbitset bd_simp_temp;
    bd_simp_temp = pbitset(new bitset ( *bd_simp)); 
    int pivot_bd = pivot(bd_simp_temp);
    bool zeroed = (pivot_bd == -1);
    while (!zeroed) {
        // Find the pbitset in M that has the same pivot as a.
        // assert(pivots[p-1].find(pivot_bd) != pivots[p-1].end());
        int conflict_idx = pivots[p-1].at(pivot_bd);
        // Add this pbitset to indices.
        dynamic_xor(bd_simp_temp, Z[p-1][conflict_idx]);
        I.push_back(conflict_idx);
        if (birth_timestamp[p-1][conflict_idx].first) {
            all_boundary = false; // Check the birth timestamps to check whether all of them are boundaries.
        }
        // Update the pivot of a.
        pivot_bd = pivot(bd_simp_temp);
        zeroed = (pivot_bd == -1);
    } 
            
}

void add_new_cycle(){
    pbitset new_pbitset = make_shared<bitset>(unique_id[p].size(), 0); // New pbitset with of size unique_id[p].size() with 1 appended at the end.
    new_pbitset -> set(id[p].at(simp));
    if (p != 0)  {
        for (auto a: I) {
            int chain_a = birth_timestamp[p-1][a].second; // Since a is a boundary, we know that the birth timestamp is non-negative and this is safe!
            dynamic_xor(new_pbitset, C[p][chain_a]); // Append a new pbitset simp + \sum_{a \in I} C^{p}[a] with birth timestamp i+1 to Z^p.
        }
    }
    // Add a wire to the bundle with timestamp i+1 and update the link from the cycle to the bundle.
    bundle[p].push_back(new_pbitset);
    timestamp[p].push_back(i+1);
    pbitset new_link = make_shared<bitset>(bundle[p].size(), 0);
    new_link -> set(bundle[p].size()-1);
    // Now, we add the new pbitset to the cycle matrix and the new link to the link matrix.
    uint pivot_new = pivot(new_pbitset);
    uint av_idx;
    if (available_pbitsets_Z[p].size() != 0) {
        av_idx = available_pbitsets_Z[p].back();
        available_pbitsets_Z[p].pop_back();
    }
    else {
        Z[p].push_back(nullptr);
        links[p].push_back(nullptr);
        birth_timestamp[p].push_back(make_pair(true, -1));
        av_idx = Z[p].size() - 1;
    }
    Z[p][av_idx] = new_pbitset;
    links[p][av_idx] = new_link;    
    birth_timestamp[p][av_idx] = make_pair(true, i+1);
    pivots[p][pivot_new] = av_idx;
    used_pbitsets_Z[p].push_back(av_idx);
}

void add_new_chain()
{
    /* Bundles and Links UPDATE: */
    bundle[p-1].push_back(bd_simp);
    timestamp[p-1].push_back(i+1);
    pbitset new_link = make_shared<bitset>(bundle[p-1].size(), 0);
    new_link -> set(bundle[p-1].size()-1);
    links[p-1][l] = new_link;
    // Populate the chain matrix and point to the boundary.
    pbitset new_pbitset = make_shared<bitset>(unique_id[p].size(), 0);
    uint av_idx;
    if (available_pbitsets_C[p].size() != 0) {
        av_idx = available_pbitsets_C[p].back();
        available_pbitsets_C[p].pop_back();
    }
    else {
        C[p].push_back(nullptr);
        av_idx = C[p].size() - 1;
    }
    C[p][av_idx] = new_pbitset;
    birth_timestamp[p-1][l] = make_pair(false, av_idx);
    used_pbitsets_C[p].push_back(av_idx);
    /* UPDATE: */
    const uint prev_pivot = pivot(Z[p-1][l]);
    Z[p-1][l] = bd_simp;
}

void reduce_boundaries(){
    vector<int> I; // Find all the pbitsets in C[p] that contain the simplex and get the indices of their corresponding boundaries.
    // alpha will be the only pbitset with negative birth timestamp such that C[p][chain_a] contains the simplex,  chain_alpha is needed to keep track of the column being added for link addition.
    int alpha, chain_alpha; 
    int smallest_pivot = unique_id[p-1].size();
    for (auto cyc_idx: used_pbitsets_Z[p-1]) {
        if (!birth_timestamp[p-1][cyc_idx].first) {// If the birth timestamp is false, then the pbitset is a boundary.
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
    for (auto cyc_idx: I) { // Add alpha to all the other pbitsets in I.
        int chain_idx = birth_timestamp[p-1][cyc_idx].second;
        if (cyc_idx != alpha) {
            dynamic_xor(Z[p-1][cyc_idx], Z[p-1][alpha]);
            dynamic_xor(links[p-1][cyc_idx], links[p-1][alpha]);
            dynamic_xor(C[p][chain_idx], C[p][chain_alpha]);
        }
    }
}

void delete_cycle(){
    pbitset z = Z[p][alpha];
    int z_pivot = pivot(z);
    pivots[p].erase(z_pivot); // Remove the pivot of Z[p][alpha] from pivots[p].
    I.erase(I.begin());
    int current_alpha = alpha; // Keeps track of the pbitset index being added in order to update links.
    if (I.size() != 0) {
        for (auto a: I) { // Iterate over the rest of I and ensure that the pivots remain distinct:
            int a_pivot = pivot(Z[p][a]);
            if (a_pivot > z_pivot) { // z does not change.
                dynamic_xor(Z[p][a], z);
                dynamic_xor(links[p][a], links[p][current_alpha]);
            }   
            else  { // We want to deep copy Z[p][a] to a temporary pbitset and then copy Z[p][a] to z.
                pbitset temp = make_shared<bitset>(*Z[p][a]);
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
}

void make_pivots_distinct(
    const std::vector<bool> *filt_op, 
    vector<pbitset> *Z_p,
    vector<pbitset> *links_p,
    vector<pbitset> *C_p,
    vector<pair <bool, int> > *birth_timestamp_p,
    PivotMap *pivots_p,
    const int prev_pivot,
    const int l
    )
{
    bool pivot_conflict;
    int a, b, pivot_a, pivot_b, current_idx, current_pivot;
    pivots_p -> erase(prev_pivot); // Remove the previous pivot
    current_pivot = pivot((*Z_p)[l]);
    if (pivots_p -> find(current_pivot) == pivots_p -> end())
    {   
        pivot_conflict = false;
        (*pivots_p)[current_pivot] = l;
    }
    else {
        pivot_conflict = true;
        // FIXME: Which one is a and which one is b?
        a = l;
        b = (*pivots_p)[current_pivot];
        pivot_a = current_pivot;
        pivot_b = current_pivot;
    }
    while (pivot_conflict)
    {
        // Check that a and b are valid indices in Z[p-1]:
        if ((*birth_timestamp_p)[a].first && !(*birth_timestamp_p)[b].first) { // Both a and b are boundaries.
            dynamic_xor((*Z_p)[a], (*Z_p)[b]);
            dynamic_xor((*links_p)[a], (*links_p)[b]);
            // Update chain matrices to reflect this addition (use the BdToChainMap):
            dynamic_xor((*C_p)[(*birth_timestamp_p)[a].second], (*C_p)[(*birth_timestamp_p)[b].second]);
            // Check for pivot conflicts again (pivot of b remains the same, but the pivot of a might have changed):
            (*pivots_p)[pivot_b] = b;
            pivot_a = pivot((*Z_p)[a]);
            if ((*pivots_p).find(pivot_a) == (*pivots_p).end())
            {
                pivot_conflict = false;
                (*pivots_p)[pivot_a] = a;
            }
            else 
            {
                if ((*pivots_p)[pivot_a] == a) {
                    pivot_conflict = false;
                }
                else 
                {
                    pivot_conflict = true;
                    b = (*pivots_p)[pivot_a];
                    pivot_b = pivot_a;
                }
            }

        }
        else if (!(*birth_timestamp_p)[a].first && (*birth_timestamp_p)[b].first) {
            dynamic_xor((*Z_p)[b], (*Z_p)[a]);
            dynamic_xor((*links_p)[b], (*links_p)[a]);
            // No need to update the chain matrices as a boundary is added to a cycle. Instead, check for pivot conflicts again:
            (*pivots_p)[pivot_a] = a;
            pivot_b = pivot((*Z_p)[b]);
            if ((*pivots_p).find(pivot_b) == (*pivots_p).end())
            {
                pivot_conflict = false;
                (*pivots_p)[pivot_b] = b;
            }
            else 
            {
                if ((*pivots_p)[pivot_b] == b) {
                    pivot_conflict = false;
                }
                else 
                {
                    pivot_conflict = true;
                    a = (*pivots_p)[pivot_b];
                    pivot_a = pivot_b;
                }
            }
        }
        else if ((*birth_timestamp_p)[a].first && !(*birth_timestamp_p)[b].first) {
            dynamic_xor((*Z_p)[a], (*Z_p)[b]);
            dynamic_xor((*links_p)[a], (*links_p)[b]);
            // Again, no need to update the chain matrices as a boundary is added to a cycle. Instead, check for pivot conflicts again:
            (*pivots_p)[pivot_b] = b;
            pivot_a = pivot((*Z_p)[a]);
            if ((*pivots_p).find(pivot_a) == (*pivots_p).end())
            {
                pivot_conflict = false;
                (*pivots_p)[pivot_a] = a;
            }
            else {
                if ((*pivots_p)[pivot_a] == a) {
                    pivot_conflict = false;
                }
                else 
                {
                    pivot_conflict = true;
                    b = (*pivots_p)[pivot_a];
                    pivot_b = pivot_a;
                }
            }
        }
        else {
            int birth_a = (*birth_timestamp_p)[a].second;
            int birth_b = (*birth_timestamp_p)[b].second;
            bool a_lessthan_b = ((birth_a == birth_b) || ((birth_a < birth_b) && ((*filt_op)[birth_b-1])) || ((birth_a > birth_b) && (!((*filt_op)[birth_a-1]))));   
            if (a_lessthan_b) {
                dynamic_xor((*Z_p)[b], (*Z_p)[a]);
                dynamic_xor((*links_p)[b], (*links_p)[a]);
                // Since both are cycles, no need to update the chain matrices. Instead, check for pivot conflicts again:
                (*pivots_p)[pivot_a] = a;
                pivot_b = pivot((*Z_p)[b]);
                if ((*pivots_p).find(pivot_b) == (*pivots_p).end())
                {
                    pivot_conflict = false;
                    (*pivots_p)[pivot_b] = b;
                }
                else 
                {
                    if ((*pivots_p)[pivot_b] == b) {
                        pivot_conflict = false;
                    }
                    else {
                        pivot_conflict = true;
                        a = (*pivots_p)[pivot_b];;
                        pivot_a = pivot_b;
                    }
                }
            }
            else {
                dynamic_xor((*Z_p)[a], (*Z_p)[b]);
                dynamic_xor((*links_p)[a], (*links_p)[b]);
                // Since both are cycles, no need to update the chain matrices. Instead, check for pivot conflicts again:
                (*pivots_p)[pivot_b] = b;
                pivot_a = pivot((*Z_p)[a]);
                if ((*pivots_p).find(pivot_a) == (*pivots_p).end())
                {
                    pivot_conflict = false;
                    (*pivots_p)[pivot_a] = a;
                }
                else 
                {
                    if ((*pivots_p)[pivot_a] == a) {
                        pivot_conflict = false;
                    }
                    else {
                        pivot_conflict = true;
                        b = (*pivots_p)[pivot_a];
                        pivot_b = pivot_a;
                    }
                }
            }
        }
    }
}

// TODO: Separate bundle and interval processing here.
void output_representatives(
    const int p, 
    const int birth, 
    const int death,
    pbitset link_interval,
    vector<int> *unique_id_p,
    vector<pbitset> *bundle,
    vector<int> *timestamp,
    std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > *persistence
    )
{
    // Here, we need to gather all the wires pertaining to the cycle l and reconstruct the representatives at each index.
    vector<tuple<int, pbitset>> wire_representatives;
    vector<tuple<int, vector<int>>> wire_representatives_ids;
    // Alternatively, just go over the non-zero values using the next operator
    size_t col_idx = link_interval -> find_first();
    while (col_idx != link_interval -> npos) {
        wire_representatives.push_back(make_tuple(timestamp[p][col_idx], bundle[p][col_idx]));
        col_idx = link_interval -> find_next(col_idx);
    }
    // Sort the representatives in the order of the timestamps.
    sort(wire_representatives.begin(), wire_representatives.end(), [&](tuple<int, pbitset> &a, tuple<int, pbitset> &b){ return (get<0>(a) < get<0>(b));});
    // Initialize a bitset pbitset of size representative_max_size with all zeros.
    pbitset current_representative = make_shared<bitset>(1, 0);
    // Add all the representatives with timestamps less than or equal to interval_birth.
    size_t wire_idx;
    for (wire_idx = 0; wire_idx < wire_representatives.size(); ++wire_idx) {
        if (get<0>(wire_representatives[wire_idx]) <= birth) {
            dynamic_xor(current_representative, get<1>(wire_representatives[wire_idx]));
        }
    }
    vector<int> indices;
    int rep_idx = current_representative -> find_first();
    while (rep_idx != current_representative -> npos) {
            indices.push_back((*unique_id_p)[rep_idx]);
            rep_idx = current_representative -> find_next(rep_idx);
    }
    wire_representatives_ids.push_back(make_tuple(birth, indices));
    for (; wire_idx < wire_representatives.size(); ++wire_idx) {
        int time = get<0>(wire_representatives[wire_idx]);
        dynamic_xor(current_representative, get<1>(wire_representatives[wire_idx]));
        vector<int> indices;
        int rep_idx = current_representative -> find_first();
        while (rep_idx != current_representative -> npos) {
                indices.push_back((*unique_id_p)[rep_idx]);
                rep_idx = current_representative -> find_next(rep_idx);
        }
        wire_representatives_ids.push_back(make_tuple(time, indices));
    }
    persistence -> push_back(make_tuple(birth, death, p, wire_representatives_ids));
}


// TODO: Extract the boundary and sum computatiion into a separate function.
}// namespace ZZREP
