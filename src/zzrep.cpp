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
#include <boost/dynamic_bitset.hpp>
#include <utility>
#include <vector>

using namespace std;

namespace ZZREP { 

// Template for hashing vectors.
template <class ElemType>
class VecHash { 
public:
    size_t operator()(const std::vector<ElemType>& v) const; 
};
template <class ElemType>
size_t VecHash<ElemType>
    ::operator()(const std::vector<ElemType>& v) const {
    std::size_t seed = 0;
    for (auto e : v) {boost::hash_combine(seed, e); }
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
// Bitsets:
typedef boost::dynamic_bitset<> bitset;
typedef shared_ptr<bitset> pbits;
// Maps and Pairs:
typedef std::map< vector<int>, int> SimplexIdMap;
typedef std::map<int, int> PivotMap;
typedef SimplexIdMap::value_type SimplexIdPair;
typedef PivotMap::value_type PivotPair;

// Data Structures:
class birth_timestamp {
public:
    // If non_bd is true, then val_idx records the birth timestamp of the non-boundary cycle,
    // otherwise, val_idx records to the index of the corresponding chain
    bool non_bd;
    int val_idx;
    birth_timestamp(bool non_bd, int val_idx) : non_bd(non_bd), val_idx(val_idx) {} 
};
class zigzag_matrices {
    public:
    vector<SimplexIdMap> id; // A collection of maps: id[p] stores the simplex ids of the p-dimensional cycles.
    vector<vector<int>> unique_id; // A collection of vectors: unique_id[p] stores the unique ids of the p-dimensional cycles.
    vector<vector<pbits> > Z; // A collection of matrices: Z[p][j] is the j-th cycle matrix of dimension p.
    vector<vector<int>> available_pbits_Z; // A collection of vectors: available_pbits_Z[p] stores the indices of the available pbitss for Z[p].
    vector<vector<int>> used_pbits_Z; // A collection of vectors: used_pbits_Z[p] stores the indices of the used pbitss for Z[p].
    vector<vector<pbits> > C;  // A collection of matrices: C[p][j] is the j-th coboundary matrix of dimension p.
    vector<vector<int>> available_pbits_C; // A collection of vectors: available_pbits_C[p] stores the indices of the available pbitss for C[p].
    vector<vector<int>> used_pbits_C; // A collection of vectors: used_pbits_C[p] stores the indices of the used pbitss for C[p].
    vector<vector<pair <bool, int> > > birth_timestamp; // A collection of vectors: birth_timestamp[p][j] is the birth timestamp of the j-th cycle of dimension p.
    vector<PivotMap> pivots;
    // Constructor:
    zigzag_matrices(const int m) 
    {
        id = vector<SimplexIdMap>(m+1, SimplexIdMap()); 
        unique_id = vector<vector<int>>(m+1, vector<int>()); 
        Z = vector<vector<pbits>>(m+1, vector<pbits>()); 
        available_pbits_Z = vector<vector<int>>(m+1, vector<int>()); 
        used_pbits_Z = vector<vector<int>>(m+1, vector<int>()); 
        C = vector<vector<pbits>>(m+1, vector<pbits>()); 
        available_pbits_C = vector<vector<int>>(m+1, vector<int>()); 
        used_pbits_C = vector<vector<int>>(m+1, vector<int>()); 
        birth_timestamp = vector<vector<pair<bool, int>>>(m+1, vector<pair<bool, int>>()); 
        pivots = vector<PivotMap>(m+1, PivotMap());
    };
};
class rep_matrices {
    public:
    vector<vector<pbits>> bundle; // A collection of wires: bundle[p] stores p-dimensional wires.
    vector<vector<int>> timestamp; // Map from the cycle matrix to the wires: links[p][j] is the collection of wires associated with Z[p][j].
    vector<vector<pbits> > links; // Represents a map from the cycle matrix to the wires: links[p][j] is the collection of wires associated with Z[p][j]; we use the same used/available pbitss as in Z.
    // Constructor:
    rep_matrices(const int m) 
    {
        bundle = vector<vector<pbits>>(m+1, vector<pbits>());
        timestamp = vector<vector<int>>(m+1, vector<int>());
        links = vector<vector<pbits>>(m+1, vector<pbits>());
    };
};

/* PRIMITIVE OPERATIONS */
int pivot(pbits a);
void dynamic_xor(pbits a, pbits b);

/* HELPER FUNCTIONS: */
void update_id(
    const int i,
    SimplexIdMap *id_p,
    vector<int> *unique_id_p,
    vector<int> *simp,
    std::map<int, int> *i_to_id_p
);
bool boundary_as_cycles(
    const int p,
    vector<int> *simp,
    pbits *bd_simp,
    vector<int> *I,
    zigzag_matrices *zz_mat
    );
void create_new_cycle(
    const int p,
    vector<int> *I,
    vector<int> *unique_id_p,
    SimplexIdMap *id_p,
    vector<int> *simp,
    vector<pair<bool, int>> *birth_timestamp_pm1,
    vector<pbits> *C_p,
    pbits *new_pbits
    );
void add_new_cycle(
    const int p,
    const int i,
    pbits *new_pbits,
    vector<pair<bool, int>> *birth_timestamp_p,
    vector<pbits> *Z,
    vector<int> *available_pbits_Z,
    vector<int> *used_pbits_Z,
    PivotMap *pivots,
    vector<pbits> *bundle,
    vector<pbits> *links,
    vector<int> *timestamp
    );
uint find_dying_cycle(
    const vector<bool> &filt_op,
    vector<pair<bool, int>> *birth_timestamp_pm1,
    vector<int> *I
);
uint add_chain_bd(
    const int i,
    const int l,
    pbits *bd_simp,
    vector<int> *unique_id,
    vector<pbits> *bundle,
    vector<int> *timestamp,
    vector<pbits> *links,
    vector<pbits> *Z_pm1,
    vector<pair<bool, int>> *birth_timestamp,
    vector<int> *available_pbits_C,
    vector<int> *used_pbits_C,
    vector<pbits> *C_p
    );
void make_pivots_distinct(
    const std::vector<bool> *filt_op, 
    vector<pbits> *Z_p,
    vector<pbits> *links_p,
    vector<pbits> *C_p,
    vector<pair <bool, int> > *birth_timestamp_p,
    PivotMap *pivots_p,
    const int l,
    const int prev_pivot
    );
void reduce_bd_update_mat(
    const int idx,
    const int i,
    vector<int> *unique_id_pm1,
    vector<pbits> *Z_pm1,
    vector<int> *used_pbits_Z_pm1,
    vector<pair<bool, int>> *birth_timestamp_pm1,
    vector<pbits> *C_p,
    vector<int> *available_pbits_C_p,
    vector<int> *used_pbits_C_p,
    vector<pbits> *bundle_pm1,
    vector<int> *timestamp_pm1,
    vector<pbits> *links_pm1
    );
uint delete_cycle(
    const int idx,
    const vector<bool> &filt_op,
    vector<pbits> *Z_p,
    vector<int> *used_pbits_Z_p,
    vector<int> *available_pbits_Z_p,
    vector<pair<bool, int>> *birth_timestamp_p,
    PivotMap *pivots_p,
    vector<pbits> *links_p
    );
void update_Z(
    const int alpha,
    vector<pbits> *Z_p,
    vector<pbits> *links_p,
    vector<int> *avaialable_pbits_Z_p,
    vector<int> *used_pbits_Z_p,
    vector<pair<bool, int>> *birth_timestamp_p
);
void output_representatives(
    const int p, 
    const int birth, 
    const int death,
    pbits link_interval,
    vector<int> *unique_id_p,
    vector<pbits> *bundle,
    vector<int> *timestamp,
    std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > *persistence
    );

void ZigzagRep::compute(
        const std::vector<vector<int> > &filt_simp, 
        const std::vector<bool> &filt_op,
        std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > *persistence, 
        std::vector <std::map<int, int>> *i_to_id,
        const int m) 
{

    persistence -> clear();
    const int n = filt_simp.size();

    // INITIALIZE DATA STRUCTURES:
    zigzag_matrices zz_mat(m);
    rep_matrices rp_mat(m);

    // COMPUTATION of zigzag persistence and the associated representatives:
    for (int i = 0; i < n; ++i) {
        vector<int> simp = filt_simp[i];
        const int p = simp.size() - 1; // p denotes the dimension of the simplex.

        if (filt_op[i]) { // INSERTION
            // Find and update the ID of the simplex.
            update_id(i, &zz_mat.id[p], &zz_mat.unique_id[p], &simp, &(*i_to_id)[p]);

            // Represent BOUNDARY AS A SUM of cycles.
            bool all_boundary = true;
            pbits bd_simp;
            vector<int> I;
            if (p != 0) all_boundary = boundary_as_cycles(p, &simp, &bd_simp, &I, &zz_mat);

            // FORWARD BIRTH.
            if (all_boundary) {
                // create NEW CYCLE
                pbits new_pbits = make_shared<bitset>(zz_mat.unique_id[p].size(), 0); // New pbits with of size unique_id[p].size() with 1 appended at the end.
                (*new_pbits).set(zz_mat.unique_id[p].size() - 1);
                if (p != 0)  {
                    for (auto a: I) {
                        int chain_a = (zz_mat.birth_timestamp[p-1])[a].second; // Since a is a boundary, we know that the birth timestamp is non-negative and this is safe!
                        dynamic_xor(new_pbits, zz_mat.C[p][chain_a]); // Append a new pbits simp + \sum_{a \in I} C^{p}[a] with birth timestamp i+1 to Z^p.
                    }
                }
                // ADD/UPDATE matrices.
                add_new_cycle(p, i, &new_pbits, &zz_mat.birth_timestamp[p], &zz_mat.Z[p], &zz_mat.available_pbits_Z[p], &zz_mat.used_pbits_Z[p], &zz_mat.pivots[p], &rp_mat.bundle[p], &rp_mat.links[p], &rp_mat.timestamp[p]);
            }   

            // FORWARD DEATH.
            else {
                // Find the DYING CYCLE:
                const uint l = find_dying_cycle(filt_op, &zz_mat.birth_timestamp[p-1], &I);

                // OUTPUT the (p − 1)-th interval [b^{p−1}[l], i] after gathering the relevant representative. */
                output_representatives(p-1, zz_mat.birth_timestamp[p-1][l].second, i, rp_mat.links[p-1][l], &(zz_mat.unique_id[p-1]), &(rp_mat.bundle[p-1]), &(rp_mat.timestamp[p-1]), persistence);
                
                // Add NEW CHAIN:
                uint prev_pivot = add_chain_bd(i, l, &bd_simp, &zz_mat.unique_id[p-1], &rp_mat.bundle[p-1], &rp_mat.timestamp[p-1], &rp_mat.links[p-1], &zz_mat.Z[p-1], &zz_mat.birth_timestamp[p-1], &zz_mat.available_pbits_C[p], &zz_mat.used_pbits_C[p], &zz_mat.C[p]);
                
                // Check for PIVOT CONFLICTS:
                make_pivots_distinct(&filt_op, &zz_mat.Z[p-1], &rp_mat.links[p-1], &zz_mat.C[p-1], &zz_mat.birth_timestamp[p-1], &zz_mat.pivots[p-1], l, prev_pivot);
            }
        } 
        else { // DELETION:
            // Find index of simplex.
            int idx = zz_mat.id[p].at(simp);
            int used_idx;
            for (used_idx = 0; used_idx < zz_mat.used_pbits_Z[p].size(); --used_idx) if (idx < zz_mat.Z[p][used_idx] -> size() && (*zz_mat.Z[p][used_idx])[idx] == 1) break;

            // BACKWARD BIRTH.  
            if (used_idx == zz_mat.used_pbits_Z[p].size()) { 
                // REDUCE Boundaries and update zz_mat and rp_mat:
                reduce_bd_update_mat(idx, i, &zz_mat.unique_id[p-1], &zz_mat.Z[p-1], &zz_mat.used_pbits_Z[p-1], &zz_mat.birth_timestamp[p-1], &zz_mat.C[p], &zz_mat.available_pbits_C[p], &zz_mat.used_pbits_C[p], &rp_mat.bundle[p-1], &rp_mat.timestamp[p-1], &rp_mat.links[p-1]);
            }
            // BACKWARD DEATH.
            else { 
                // REMOVE the simplex from C[p]:
                int simp_col_idx = zz_mat.used_pbits_Z[p][used_idx];
                for (auto chain_idx : zz_mat.used_pbits_C[p]) if (idx < zz_mat.C[p][chain_idx] -> size() && (*zz_mat.C[p][chain_idx])[idx] == 1) dynamic_xor(zz_mat.C[p][chain_idx], zz_mat.Z[p][simp_col_idx]);  // Update C[p] so that no pbitss contain the simplex.
                
                // REMOVE the simplex from Z[p]:
                uint alpha = delete_cycle(idx, filt_op, &zz_mat.Z[p], &zz_mat.used_pbits_Z[p], &zz_mat.available_pbits_Z[p], &zz_mat.birth_timestamp[p], &zz_mat.pivots[p], &rp_mat.links[p]);

                // OUTPUT the p-th interval [birth_timestamp_p[alpha], i] after gathering the relevant representative.
                uint birth = zz_mat.birth_timestamp[p][alpha].second;
                output_representatives(p, birth, i, rp_mat.links[p][alpha], &(zz_mat.unique_id[p]), &(rp_mat.bundle[p]), &(rp_mat.timestamp[p]), persistence);

                // UPDATE Z[p] and links[p]. 
                update_Z(alpha, &zz_mat.Z[p], &rp_mat.links[p], &zz_mat.available_pbits_Z[p], &zz_mat.used_pbits_Z[p], &zz_mat.birth_timestamp[p]);
            }
        }
    }
    // POST-PROCESSING:
    for (size_t p = 0; p <= m; p++) {
        for (auto a: zz_mat.used_pbits_Z[p]) {
            if (zz_mat.birth_timestamp[p][a].first) {
                uint birth = zz_mat.birth_timestamp[p][a].second;
                output_representatives(p, birth, n, rp_mat.links[p][a], &(zz_mat.unique_id[p]), &(rp_mat.bundle[p]), &(rp_mat.timestamp[p]), persistence);
            }
        }
    }
}

/* IMPLEMENTATION OF PRIMITIVES: */
// Goes over a vector of ints and finds the position of the last non-zero element.
int pivot (pbits a)
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
void dynamic_xor(pbits a, pbits b)
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

/* IMPLEMENTATION OF HELPER FUNCTIONS. */
void update_id(
    const int i,
    SimplexIdMap *id_p,
    vector<int> *unique_id_p,
    vector<int> *simp,
    std::map<int, int> *i_to_id_p
)
{
    int unique_id_simp; 
    (*id_p).find(*simp) != (*id_p).end() ? unique_id_simp = (*unique_id_p)[(*id_p).at(*simp)] : unique_id_simp = unique_id_p -> size();
    (*id_p)[*simp] = unique_id_p -> size(); // Change the id of the simplex to the new id.
    unique_id_p -> push_back(unique_id_simp);
    (*i_to_id_p)[i] = unique_id_simp;
}
bool boundary_as_cycles(
    const int p,
    vector<int> *simp,
    pbits *bd_simp,
    vector<int> *I,
    zigzag_matrices *zz_mat
){
    bool all_boundary = true;
    for (int i = 0; i <= p; i++) {
        // Remove the i-th vertex.
        vector<int> boundary_simplex = *simp;
        boundary_simplex.erase(boundary_simplex.begin() + i);
        int idx = (*zz_mat).id[p].at(boundary_simplex);
        (*bd_simp) -> set(idx);
    }
    // Find the pbits in Z[p-1] that sum to bd_simp:
    pbits bd_simp_temp;
    bd_simp_temp = pbits(new bitset (**bd_simp)); // Copy bd_simp to a temporary pbits.
    int pivot_bd = pivot(bd_simp_temp);
    bool zeroed = (pivot_bd == -1);
    while (!zeroed) {
        int conflict_idx = (*zz_mat).pivots[p-1].at(pivot_bd);
        dynamic_xor(bd_simp_temp, (*zz_mat).Z[p-1][conflict_idx]); // Add the conflicting pbits to bd_simp.
        (*I).push_back(conflict_idx);
        if ((*zz_mat).birth_timestamp[p-1][conflict_idx].first) all_boundary = false; // Check the birth timestamps to check whether all of them are boundaries.
        // Update the pivot of a.
        pivot_bd = pivot(bd_simp_temp);
        zeroed = (pivot_bd == -1);
    } 
    return all_boundary;
}

void add_new_cycle(
    const int p,
    const int i,
    pbits *new_pbits,
    zigzag_matrices *zz_mat,
    rep_matrices *rp_mat
){
    // Add a wire to the bundle with timestamp i+1 and update the link from the cycle to the bundle.
    (*rp_mat).bundle[p].push_back(*new_pbits);
    (*rp_mat).timestamp[p].push_back(i+1);
    pbits new_link = make_shared<bitset>((*rp_mat).bundle[p].size(), 0);
    new_link -> set((*rp_mat).bundle[p].size()-1);
    // Now, we add the new pbits to the cycle matrix and the new link to the link matrix.
    uint pivot_new = pivot(*new_pbits);
    uint av_idx;
    if ((*zz_mat).available_pbits_Z[p].size() != 0) {
        av_idx = (*zz_mat).available_pbits_Z[p].back();
        (*zz_mat).available_pbits_Z[p].pop_back();
    }
    else {
        (*zz_mat).Z[p].push_back(nullptr);
        (*rp_mat).links[p].push_back(nullptr);
        (*zz_mat).birth_timestamp[p].push_back(make_pair(true, -1));
        av_idx = (*zz_mat).Z[p].size() - 1;
    }
    (*zz_mat).Z[p][av_idx] = *new_pbits;
    (*rp_mat).links[p][av_idx] = new_link;    
    (*zz_mat).birth_timestamp[p][av_idx] = make_pair(true, i+1);
    (*zz_mat).pivots[p][pivot_new] = av_idx;
    (*zz_mat).used_pbits_Z[p].push_back(av_idx);
}
uint find_dying_cycle(
    const vector<bool> &filt_op,
    vector<pair<bool, int>> *birth_timestamp_pm1,
    vector<int> *I
)
{
    vector<int> J;
    for (int a: *I) if ((*birth_timestamp_pm1)[a].first) J.push_back(a);  // Gather all non-boundary cycles.
    sort(J.begin(), J.end());
    int k;
    for (k = J.size()-1; k >= 0; k--) if (filt_op[(*birth_timestamp_pm1)[J[k]].second - 1]) break; // conditional checks if arrow points forward.
    const int l = J[l];
    return l;
}
void make_pivots_distinct(
    const std::vector<bool> *filt_op, 
    vector<pbits> *Z_p,
    vector<pbits> *links_p,
    vector<pbits> *C_p,
    vector<pair <bool, int> > *birth_timestamp_p,
    PivotMap *pivots_p,
    const int l,
    const int prev_pivot
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
uint add_chain_bd(
    const int i,
    const int l,
    pbits *bd_simp,
    vector<int> *unique_id,
    vector<pbits> *bundle,
    vector<int> *timestamp,
    vector<pbits> *links,
    vector<pbits> *Z_pm1,
    vector<pair<bool, int>> *birth_timestamp,
    vector<int> *available_pbits_C,
    vector<int> *used_pbits_C,
    vector<pbits> *C_p
)
{
    /* Bundles and Links UPDATE: */
    bundle -> push_back(*bd_simp);
    timestamp -> push_back(i+1);
    pbits new_link = make_shared<bitset>(bundle -> size(), 0);
    new_link -> set(bundle ->size()-1);
    (*links)[l] = new_link;
    // Populate the chain matrix and point to the boundary.
    pbits new_chain = make_shared<bitset>(unique_id -> size(), 0);
    uint av_idx;
    if ((*available_pbits_C).size() != 0) {
        av_idx = (*available_pbits_C).back();
        (*available_pbits_C).pop_back();
    }
    else {
        C_p -> push_back(nullptr);
        av_idx = C_p -> size() - 1;
    }
    (*C_p)[av_idx] = new_chain;
    (*birth_timestamp)[l] = make_pair(false, av_idx);
    used_pbits_C -> push_back(av_idx);
    /* UPDATE: */
    const uint prev_pivot = pivot((*Z_pm1)[l]);
    (*Z_pm1)[l] = *bd_simp;
    return prev_pivot;
}
void reduce_bd_update_mat(
    const int idx,
    const int i,
    vector<int> *unique_id_pm1,
    vector<pbits> *Z_pm1,
    vector<int> *used_pbits_Z_pm1,
    vector<pair<bool, int>> *birth_timestamp_pm1,
    vector<pbits> *C_p,
    vector<int> *available_pbits_C_p,
    vector<int> *used_pbits_C_p,
    vector<pbits> *bundle_pm1,
    vector<int> *timestamp_pm1,
    vector<pbits> *links_pm1
){
    vector<int> I; // Find all the pbitss in C[p] that contain the simplex and get the indices of their corresponding boundaries.
    // alpha will be the only pbits with negative birth timestamp such that C[p][chain_a] contains the simplex,  chain_alpha is needed to keep track of the column being added for link addition.
    int alpha, chain_alpha; 
    int smallest_pivot = unique_id_pm1 -> size();
    for (auto cyc_idx: *used_pbits_Z_pm1) {
        if (!(*birth_timestamp_pm1)[cyc_idx].first) {// If the birth timestamp is false, then the pbits is a boundary.
            int chain_idx = (*birth_timestamp_pm1)[cyc_idx].second;
            if (idx < (*C_p)[chain_idx] -> size() && (*(*C_p)[chain_idx])[idx]==1) {
                I.push_back(cyc_idx);
                int pivot_cyc_idx = pivot((*Z_pm1)[cyc_idx]);
                if (pivot_cyc_idx < smallest_pivot) {
                    smallest_pivot = pivot_cyc_idx;
                    alpha = cyc_idx;
                    chain_alpha = chain_idx;
                }
            }
        }
    }
    for (auto cyc_idx: I) { // Add alpha to all the other pbitss in I.
        int chain_idx = (*birth_timestamp_pm1)[cyc_idx].second;
        if (cyc_idx != alpha) {
            dynamic_xor((*Z_pm1)[cyc_idx], (*Z_pm1)[alpha]);
            dynamic_xor((*links_pm1)[cyc_idx], (*links_pm1)[alpha]);
            dynamic_xor((*C_p)[chain_idx], (*C_p)[chain_alpha]);
        }
    }

    // UPDATE Chain matrices
    (*C_p)[chain_alpha] = nullptr;
    (*available_pbits_C_p).push_back(chain_alpha); 
    (*used_pbits_C_p).erase(remove((*used_pbits_C_p).begin(), (*used_pbits_C_p).end(), chain_alpha), (*used_pbits_C_p).end());
    (*birth_timestamp_pm1)[alpha] = make_pair(true, i+1);
    
    // UPDATE Bundles and Links:
    (*bundle_pm1).push_back((*Z_pm1)[alpha]);
    (*timestamp_pm1).push_back(i+1);
    pbits new_link = make_shared<bitset>(bundle_pm1 -> size(), 0);
    new_link -> set(bundle_pm1 -> size()-1);
    (*links_pm1)[alpha] = new_link;
}
uint delete_cycle(
    const int idx,
    const vector<bool> &filt_op,
    vector<pbits> *Z_p,
    vector<int> *used_pbits_Z_p,
    vector<int> *available_pbits_Z_p,
    vector<pair<bool, int>> *birth_timestamp_p,
    PivotMap *pivots_p,
    vector<pbits> *links_p
){

    vector<int> I; 
    for (auto col_idx: (*used_pbits_Z_p))
    {
        pbits column = (*Z_p)[col_idx];
        if (idx < (column -> size()) && (*column)[idx] == 1)
        {
            I.push_back(col_idx); // Gather indices of pbitss that contain the simplex. 
        }
    }  
    // sort I in the order of the birth timestamps where the order is the total order as above without using the sort function.
    sort(I.begin(), I.end(), [&](int &a, int &b){
        int birth_a = (*birth_timestamp_p)[a].second; int birth_b = (*birth_timestamp_p)[b].second;
        return ((birth_a == birth_b) || ((birth_a < birth_b) && (filt_op[birth_b-1])) || ((birth_a > birth_b) && (!(filt_op[birth_a-1]))));});
    const int alpha = I[0]; // The pbits to be deleted is the first pbits in I.
    pbits z = (*Z_p)[alpha];
    int z_pivot = pivot(z);
    pivots_p -> erase(z_pivot); // Remove the pivot of Z[p][alpha] from pivots[p].
    I.erase(I.begin());
    int current_alpha = alpha; // Keeps track of the pbits index being added in order to update links.
    if (I.size() != 0) {
        for (auto a: I) { // Iterate over the rest of I and ensure that the pivots remain distinct:
            int a_pivot = pivot((*Z_p)[a]);
            if (a_pivot > z_pivot) { // z does not change.
                dynamic_xor((*Z_p)[a], z);
                dynamic_xor((*links_p)[a], (*links_p)[current_alpha]);
            }   
            else  { // We want to deep copy Z[p][a] to a temporary pbits and then copy Z[p][a] to z.
                pbits temp = make_shared<bitset>(*(*Z_p)[a]);
                int temp_pivot = a_pivot;
                dynamic_xor((*Z_p)[a], z);
                dynamic_xor((*links_p)[a], (*links_p)[current_alpha]);
                a_pivot = z_pivot;
                // Update the information for z's.
                z = temp;
                current_alpha = a;
                z_pivot = temp_pivot;
                // Update the pivot of Z[p][a] in pivots[p].
                (*pivots_p)[a_pivot] = a;
                pivots_p -> erase(z_pivot);
            }
        }
    }
    return alpha;
}
void output_representatives(
    const int p, 
    const int birth, 
    const int death,
    pbits link_interval,
    vector<int> *unique_id_p,
    vector<pbits> *bundle,
    vector<int> *timestamp,
    std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > *persistence
    )
{
    // Here, we need to gather all the wires pertaining to the cycle l and reconstruct the representatives at each index.
    vector<tuple<int, pbits>> wire_representatives;
    vector<tuple<int, vector<int>>> wire_representatives_ids;
    // Alternatively, just go over the non-zero values using the next operator
    size_t col_idx = link_interval -> find_first();
    while (col_idx != link_interval -> npos) {
        wire_representatives.push_back(make_tuple(timestamp[p][col_idx], bundle[p][col_idx]));
        col_idx = link_interval -> find_next(col_idx);
    }
    // Sort the representatives in the order of the timestamps.
    sort(wire_representatives.begin(), wire_representatives.end(), [&](tuple<int, pbits> &a, tuple<int, pbits> &b){ return (get<0>(a) < get<0>(b));});
    // Initialize a bitset pbits of size representative_max_size with all zeros.
    pbits current_representative = make_shared<bitset>(1, 0);
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
void update_Z(
    const int alpha,
    vector<pbits> *Z_p,
    vector<pbits> *links_p,
    vector<int> *avaialable_pbits_Z_p,
    vector<int> *used_pbits_Z_p,
    vector<pair<bool, int>> *birth_timestamp_p
){
    (*Z_p)[alpha] = nullptr;
    (*links_p)[alpha] = nullptr;
    avaialable_pbits_Z_p -> push_back(alpha);       // Add this to the available pbitss for adding new cycles:
    used_pbits_Z_p -> erase(remove(used_pbits_Z_p -> begin(), used_pbits_Z_p -> end(), alpha), used_pbits_Z_p -> end());
    (*birth_timestamp_p)[alpha] = make_pair(false, -1); // We need to assign this to be invalid (redundant).
}

}// namespace ZZREP
