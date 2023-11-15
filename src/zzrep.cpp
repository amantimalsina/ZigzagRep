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
/* Template for hashing vectors: */
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
/* DATA STRUCTURES */
// Bitsets:
typedef boost::dynamic_bitset<> bitset;
typedef shared_ptr<bitset> chain;
// Maps and Pairs:
typedef std::map< vector<int>, int> SimplexIdMap;
typedef std::map<int, int> PivotMap;
typedef SimplexIdMap::value_type SimplexIdPair;
typedef PivotMap::value_type PivotPair;

// Cycle Record:
struct cycle_record {
    bool non_bd;
    int chain_idx; // The index of the chain in C[p] whose boundary is this cycle.
    int timestamp;
    cycle_record(bool non_bd, int chain_idx, int timestamp) : non_bd(non_bd), chain_idx(chain_idx), timestamp(timestamp) {}
};

// Collection of Matrices:
class zigzag_matrices {
    public:
    vector<SimplexIdMap> id; //  The current integral id for each p-simplex.
    vector<vector<int> > unique_id; // The unique integral id for each p-simplex.
    vector<vector<chain> > Z; // A cycle matrix Z[p] for each dimension p that is maintained such that each chain of Z[p] represents a p-cycle in K_i.
    vector<vector<int> > available_chain_Z; // A list of available chains in Z[p] for each dimension p.
    vector<vector<int> > used_chain_Z; // A list of used chains in Z[p] for each dimension p.
    vector<vector<chain> > C; // A chain matrix C[p] for each dimension p so that the chains in C[p] are associated with the the boundary cycles in Z[p-1].
    vector<vector<int> > available_chain_C; // A list of available chains in Z[p] for each dimension p.
    vector<vector<int> > used_chain_C; // A list of used chains in Z[p] for each dimension p.
    vector<vector<cycle_record > > Cycle_Record; // For each chain Z[p][j], a birth timestamp is maintained, where a false value implies that the chain pertains to a boundary.
    vector<PivotMap> pivots; // For each chain Z[p][j] we have a unique pivot entry; we then store pivots[p][i] = j where j = -1 implies that no chain exists in Z[p-1] with pivot i.
    /* Constructor */
    zigzag_matrices(const int m) {
        id = vector<SimplexIdMap>(m+1, SimplexIdMap());
        unique_id = vector<vector<int> >(m+1, vector<int>());
        Z = vector<vector<chain> >(m+1, vector<chain>());
        available_chain_Z = vector<vector<int> >(m+1, vector<int>());
        used_chain_Z = vector<vector<int> >(m+1, vector<int>());
        C = vector<vector<chain> >(m+1, vector<chain>());
        available_chain_C = vector<vector<int> >(m+1, vector<int>());
        used_chain_C = vector<vector<int> >(m+1, vector<int>());
        Cycle_Record = vector<vector<cycle_record > >(m+1, vector<cycle_record> ());
        pivots = vector<PivotMap>(m+1, PivotMap());
    }
};
class representative_matrices {
    public:
    vector<vector<chain> > bundle; // A collection of wires: bundle[p] stores p-dimensional wires.
    vector<vector<int> > timestamp; // Map from the cycle matrix to the wires: links[p][j] is the collection of wires associated with Z[p][j].
    vector<vector<chain> > links; // Represents a map from the cycle matrix to the wires: links[p][j] is the collection of wires associated with Z[p][j]; we use the same used/available chains as in Z.
    /* Constructor */
    representative_matrices(const int m) {
        bundle = vector<vector<chain> >(m+1, vector<chain>());
        timestamp = vector<vector<int> >(m+1, vector<int>());
        links = vector<vector<chain> >(m+1, vector<chain>());
    }
};

/* PRIMITIVE OPERATIONS: */
int pivot(chain a);
void dynamic_xor(chain a, chain b);

/* HELPER FUNCTIONS: */
void update_id(
    const int i,
    SimplexIdMap *id_p,
    vector<int> *unique_id_p,
    vector<int> *simp,
    std::map<int, int> *i_to_id_p
    );
bool compute_boundary(
    const int p,
    chain &bd_simp,
    vector<int> &simp,
    vector<int> &I,
    zigzag_matrices &zz_mat
    );
// FORWARD BIRTH:
void create_new_cycle(
    zigzag_matrices &zz_mat,
    chain &new_chain,
    vector<int> &I,
    vector<int> &simp,
    const int p,
    const int i
    );
void add_new_cycle(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    chain &new_chain,
    vector<int> &simp,
    const int p,
    const int i
    );
// FORWARD DEATH:
int find_dying_cycle(
    const vector<bool> &filt_op,
    zigzag_matrices &zz_mat,
    vector<int> &I,
    const int p
    );
void add_chain_bd(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    chain &bd_simp,
    vector<int> &simp,
    const int p,
    const int i,
    const int l
);
void make_pivots_distinct(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    const vector<bool> &filt_op,
    const int p,
    const int l
);
// BACKWARD BIRTH:
void reduce_bd_update(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    const int p,
    const int i,
    const int idx
);
// BACKWARD DEATH:
uint delete_cycle(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    const vector<bool> &filt_op,
    const int p,
    const int idx
);
void delete_update(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    const int p,
    const int alpha,
    const int idx
);
// OUTPUT:
void output_representatives(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    vector<tuple<int, int, int, vector<tuple<int, vector<int>>>>> *persistence,
    const int p,
    const int i,
    const int a,
    const int birth,
    const int death
);    


void ZigzagRep::compute(
        const std::vector<vector<int> > &filt_simp, 
        const std::vector<bool> &filt_op,
        std::vector <std::tuple <int, int, int, std::vector<std::tuple<int, std::vector<int>>> > > *persistence, 
        std::vector <std::map<int, int>> *i_to_id,
        int m) {
    
    persistence -> clear();
    int n = filt_simp.size();

    /* DECLARATION & INITIALIZATION of Data Structures: */
    zigzag_matrices zz_mat(m);
    representative_matrices rep_mat(m);

    /* COMPUTATION of zigzag persistence and the associated representatives: */
    for (int i = 0; i < n; ++i) {
        vector<int> simp = filt_simp[i];
        int p = simp.size() - 1; // p denotes the dimension of the simplex.

        if (filt_op[i]) { // INSERTION
            // Find and update ID of the simplex.
            update_id(i, &zz_mat.id[p], &zz_mat.unique_id[p], &simp, &(*i_to_id)[p]);

            // Represent the BOUNDARY as a SUM OF CYLES.
            bool all_boundary = true;
            chain bd_simp;
            vector<int> I;
            if (p != 0) all_boundary = compute_boundary(p, bd_simp, simp, I, zz_mat);

            // FORWARD BIRTH.
            if (all_boundary)
            {
                // Create NEW CYCLE:
                chain new_chain; 
                create_new_cycle(zz_mat, new_chain, I, simp, p, i);

                // ADD & UPDATE:
                add_new_cycle(zz_mat, rep_mat, new_chain, simp, p, i);
            }
            // FORWARD DEATH.
            else {
                // Find the DYING CYCLE:
                int l = find_dying_cycle(filt_op, zz_mat, I, p);

                // OUTPUT the (p − 1)-th interval [b^{p−1}[l], i] after gathering the relevant representative:
                uint birth = zz_mat.Cycle_Record[p-1][l].timestamp;
                output_representatives(zz_mat, rep_mat, persistence, p-1, i, l, birth, i);

                // Add NEW CHAIN:
                add_chain_bd(zz_mat, rep_mat, bd_simp, simp, p, i, l);

                // Avoid PIVOT CONFLICTS:
                make_pivots_distinct(zz_mat, rep_mat, filt_op, p, l);
            }
        } 
        else { // DELETION:
            // Find simplex in Z.
            int idx = zz_mat.id[p].at(simp);
            uint simp_col_idx;
            bool cycle_found = false;
            for (auto used_col_idx: zz_mat.used_chain_Z[p]) {
                if (idx < zz_mat.Z[p][used_col_idx] -> size() && (*zz_mat.Z[p][used_col_idx])[idx] == 1) {
                    simp_col_idx = used_col_idx;
                    cycle_found = true;
                    break;
                }
            }

            // BACKWARD BIRTH:
            if (!cycle_found)
            {   
                // REDUCE boundaries to get solitary boundary and update data structures.
                reduce_bd_update(zz_mat, rep_mat, p, i, idx);
            }
            // BACKWARD DEATH:
            else
            {
                // REMOVE simplex from C[p].
                for (auto chain_idx : zz_mat.used_chain_C[p]) if (idx < zz_mat.C[p][chain_idx] -> size() && (*zz_mat.C[p][chain_idx])[idx] == 1) dynamic_xor(zz_mat.C[p][chain_idx], zz_mat.Z[p][simp_col_idx]);

                // REMOVE simplex from Z[p]
                uint alpha = delete_cycle(zz_mat, rep_mat, filt_op, p, idx);

                // OUTPUT the p-th interval [Cycle_Record_p[alpha], i] after gathering the relevant representative.
                uint birth = zz_mat.Cycle_Record[p][alpha].timestamp;
                output_representatives(zz_mat, rep_mat, persistence, p, i, alpha, birth, i);

                // UPDATE zz_mat and rep_mat.
                delete_update(zz_mat, rep_mat, p, alpha, idx);
            }
        }
        // Assert that for each column in the cycle matrix, summing up the links should give the corresponding column back:
        for (size_t p = 0; p <= m; p++) {
            for (auto col_idx: zz_mat.used_chain_Z[p]) {
                chain link = rep_mat.links[p][col_idx];
                chain current_representative = make_shared<bitset>(1, 0);
                size_t bundle_idx = link -> find_first();
                while (bundle_idx != link -> npos) {
                    dynamic_xor(current_representative, rep_mat.bundle[p][bundle_idx]);
                    bundle_idx = link -> find_next(bundle_idx);
                }
                // Now check that we get the same column as the column in Z[p]:
                size_t nonzero_idx = current_representative -> find_first();
                while (nonzero_idx != current_representative -> npos) {
                    if ((*zz_mat.Z[p][col_idx])[nonzero_idx] == 1) {
                        // Print the column:
                        cout << "Column: ";
                        for (size_t i = 0; i < zz_mat.unique_id[p].size(); i++) {
                            if ((*zz_mat.Z[p][col_idx])[i] == 1) {
                                cout << zz_mat.unique_id[p][i] << " ";
                            }
                        }
                        cout << endl;
                        // Print the representative:
                        cout << "Representative: ";
                        for (size_t i = 0; i < zz_mat.unique_id[p].size(); i++) {
                            if ((*current_representative)[i] == 1) {
                                cout << zz_mat.unique_id[p][i] << " ";
                            }
                        }
                        cout << endl;
                        assert(false);
                    }
                    nonzero_idx = current_representative -> find_next(nonzero_idx);
                }
            }
        }
    }
    // POST-PROCESSING: 
    for (size_t p = 0; p <= m; p++) {
        for (auto a: zz_mat.used_chain_Z[p]) {
            if (zz_mat.Cycle_Record[p][a].non_bd) {
                uint birth = zz_mat.Cycle_Record[p][a].timestamp;
                output_representatives(zz_mat, rep_mat, persistence, p, n, a, birth, n);            
            }
        }
    }
}

/* IMPLEMENTATION OF HELPER FUNCTIONS: */
void update_id(
    const int i,
    SimplexIdMap *id_p,
    vector<int> *unique_id_p,
    vector<int> *simp,
    std::map<int, int> *i_to_id_p
)
{
    int unique_id_simp;
    if ((*id_p).find(*simp) != (*id_p).end()) { // If the simplex is being re-inserted.
        unique_id_simp = (*unique_id_p)[(*id_p).at(*simp)]; 
    }
    else {
        unique_id_simp = unique_id_p -> size();
    }
    (*id_p)[*simp] = unique_id_p -> size(); // Change the id of the simplex to the new id.
    unique_id_p -> push_back(unique_id_simp);
    (*i_to_id_p)[i] = unique_id_simp;   
}

bool compute_boundary(
    const int p,
    chain &bd_simp,
    vector<int> &simp,
    vector<int> &I,
    zigzag_matrices &zz_mat
    )
{
    bool all_boundary = true;
    bd_simp = make_shared<bitset>(zz_mat.unique_id[p-1].size(), 0);
    for (int i = 0; i <= p; i++) {
        // Remove the i-th vertex.
        vector<int> boundary_simplex = simp;
        boundary_simplex.erase(boundary_simplex.begin() + i);
        int idx = zz_mat.id[p-1].at(boundary_simplex);
        bd_simp -> set(idx);
    }
    /* Find the columns in Z[p-1] that sum to bd_simp: */
    chain bd_simp_temp;
    bd_simp_temp = chain(new bitset ( *bd_simp)); 
    int pivot_bd = pivot(bd_simp_temp);
    bool zeroed = (pivot_bd == -1);
    while (!zeroed) {
        // Find the chain in M that has the same pivot as a.
        int conflict_idx = zz_mat.pivots[p-1].at(pivot_bd);
        // Add this chain to indices.
        dynamic_xor(bd_simp_temp, zz_mat.Z[p-1][conflict_idx]);
        I.push_back(conflict_idx);
        if (zz_mat.Cycle_Record[p-1][conflict_idx].non_bd) {
            all_boundary = false; // Check the birth timestamps to check whether all of them are boundaries.
        }
        // Update the pivot of a.
        pivot_bd = pivot(bd_simp_temp);
        zeroed = (pivot_bd == -1);
    }
    return all_boundary; 
}

void create_new_cycle(
    zigzag_matrices &zz_mat,
    chain &new_chain,
    vector<int> &I,
    vector<int> &simp,
    const int p,
    const int i
)
{
    new_chain = make_shared<bitset>(zz_mat.unique_id[p].size(), 0);
    new_chain -> set(zz_mat.id[p].at(simp));
    if (p != 0)  {
        for (auto a: I) {
            int chain_a = zz_mat.Cycle_Record[p-1][a].chain_idx; // Since a is a boundary, we know that the birth timestamp is non-negative and this is safe!
            dynamic_xor(new_chain, zz_mat.C[p][chain_a]);
        }
    }
}

void add_new_cycle(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    chain &new_chain,
    vector<int> &simp,
    const int p,
    const int i
)
{
    uint pivot_new = pivot(new_chain);
    // Add a wire to the bundle with timestamp i+1 and update the link from the cycle to the bundle.
    rep_mat.bundle[p].push_back(new_chain);
    rep_mat.timestamp[p].push_back(i+1);
    chain new_link;
    new_link = make_shared<bitset>(rep_mat.bundle[p].size(), 0);
    new_link -> set(rep_mat.bundle[p].size()-1);
    // Now, we add the new chain to the cycle matrix and the new link to the link matrix.
    if (zz_mat.available_chain_Z[p].size() != 0) 
    {
        uint av_idx = zz_mat.available_chain_Z[p].back();
        zz_mat.Z[p][av_idx] = new_chain;
        rep_mat.links[p][av_idx] = new_link;    
        zz_mat.Cycle_Record[p][av_idx] = cycle_record(true, -1, i + 1);
        zz_mat.pivots[p][pivot_new] = av_idx;
        // Update the used and available chains:
        zz_mat.used_chain_Z[p].push_back(av_idx);
        zz_mat.available_chain_Z[p].pop_back();
    }
    else 
    {
        zz_mat.Z[p].push_back(new_chain);
        rep_mat.links[p].push_back(new_link);
        zz_mat.Cycle_Record[p].push_back(cycle_record(true, -1, i+1));
        zz_mat.pivots[p][pivot_new] = zz_mat.Z[p].size() - 1;
        // Update the used and available chains:
        zz_mat.used_chain_Z[p].push_back(zz_mat.Z[p].size() - 1);
    }
}
int find_dying_cycle(
    const vector<bool> &filt_op,
    zigzag_matrices &zz_mat,
    vector<int> &I,
    const int p
)
{
    /* Let J consist of indices in I whose corresponding chains in Z[p−1] have non-negative birth timestamps. */ 
    vector<int> J;
    for (int a: I) {
        if (zz_mat.Cycle_Record[p-1][a].non_bd) { // Gather all the non-boundary cycles.
            J.push_back(a);
        }
    }
    sort(J.begin(), J.end());
    // Check if arrow at b^{p−1}[c]−1 points backward for all c in J
    bool arrow_backward = true;        
    int l;
    for (int j_idx = J.size()-1; j_idx >= 0; j_idx--) {
        int birth_j = zz_mat.Cycle_Record[p-1][J[j_idx]].timestamp;
        if (filt_op[birth_j-1]) {
                l = J[j_idx]; // l will be the largest c in J if the arrow at {b^{p−1}[c]−1} points forward.
                arrow_backward = false;
                /*
                // Assert that l is the largest index in J such that the arrow at {b^{p−1}[c]−1} points forward.
                for (int k_idx =  J.size()-1; k_idx >= 0; k_idx--) {
                    int birth_k = zz_mat.Cycle_Record[p-1][J[k_idx]].timestamp;
                    if (filt_op[birth_k-1]) {
                        assert(k_idx <= l);
                    }
                }*/
                break;
        }
    }    
    if (arrow_backward) { // If the arrow at {b^{p−1}[c]−1} points backward for all c in J, let l be the smallest index in J.  
        l = *J.begin();
        /*
        // assert that the arrows in J all point backward.
        for (auto j: J) {
            int birth_j = zz_mat.Cycle_Record[p-1][j].timestamp;
            assert(!filt_op[birth_j-1]);
        }
        */
    }
    return l;
}
void add_chain_bd(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    chain &bd_simp,
    vector<int> &simp,
    const int p,
    const int i,
    const int l
)
{
    /* Bundles and Links UPDATE: */
    rep_mat.bundle[p-1].push_back(bd_simp);
    rep_mat.timestamp[p-1].push_back(i+1);
    chain new_link = make_shared<bitset>(rep_mat.bundle[p-1].size(), 0);
    new_link -> set(rep_mat.bundle[p-1].size()-1);
    rep_mat.links[p-1][l] = new_link;
    /* UPDATE: */
    uint prev_pivot = pivot(zz_mat.Z[p-1][l]);
    zz_mat.Z[p-1][l] = bd_simp;
    // Set C[p][l] = simp and update the boundary-to-chain map.
    chain current_chain;
    current_chain = make_shared<bitset>(zz_mat.unique_id[p].size(), 0);
    current_chain -> set(zz_mat.id[p].at(simp)); // Set last entry to 1.
    if (zz_mat.available_chain_C[p].size() != 0) 
    {
        uint av_idx = zz_mat.available_chain_C[p].back();
        zz_mat.C[p][av_idx] = current_chain;
        zz_mat.Cycle_Record[p-1][l] = cycle_record(false, av_idx, i+1);
        // Update the used and available chains:
        zz_mat.used_chain_C[p].push_back(av_idx);
        zz_mat.available_chain_C[p].pop_back();
    }
    else 
    {
        zz_mat.C[p].push_back(current_chain);
        uint av_idx = zz_mat.C[p].size() - 1;
        zz_mat.Cycle_Record[p-1][l] = cycle_record(false, av_idx, i+1);
        // Update the used and available chains:
        zz_mat.used_chain_C[p].push_back(av_idx);
    }
    zz_mat.pivots[p-1].erase(prev_pivot); // Remove the previous pivot
}
void make_pivots_distinct(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    const vector<bool> &filt_op,
    const int p,
    const int l
)
{
    bool pivot_conflict;
    int a, b, pivot_a, pivot_b, current_idx, current_pivot;
    current_pivot = pivot(zz_mat.Z[p-1][l]);
    if (zz_mat.pivots[p-1].find(current_pivot) == zz_mat.pivots[p-1].end()) {   
        pivot_conflict = false;
        zz_mat.pivots[p-1][current_pivot] = l;
    }
    else {
        pivot_conflict = true;
        a = l;
        b = zz_mat.pivots[p-1][current_pivot];;
        pivot_a = current_pivot;
        pivot_b = current_pivot;
    }
    while (pivot_conflict)
    {
        bool a_lessthan_b; // true implies we add a to b and false implies we add b to a.
        bool both_boundary = false;
        if (!zz_mat.Cycle_Record[p-1][a].non_bd && zz_mat.Cycle_Record[p-1][b].non_bd) {
            a_lessthan_b = true;
        }
        else if (zz_mat.Cycle_Record[p-1][a].non_bd && !zz_mat.Cycle_Record[p-1][b].non_bd) {
            a_lessthan_b = false;
        }
        else 
        {
            int birth_a = zz_mat.Cycle_Record[p-1][a].timestamp;
            int birth_b = zz_mat.Cycle_Record[p-1][b].timestamp;
            a_lessthan_b = ((birth_a == birth_b) || ((birth_a < birth_b) && (filt_op[birth_b-1])) || ((birth_a > birth_b) && (!(filt_op[birth_a-1]))));   
            both_boundary = (!zz_mat.Cycle_Record[p-1][a].non_bd && !zz_mat.Cycle_Record[p-1][b].non_bd);
            if (both_boundary) {
                a_lessthan_b = false;
            }
        }
        if (a_lessthan_b) {
            int temp = b;
            int pivot_temp = pivot_b;
            b = a;
            pivot_b = pivot_a;
            a = temp;
            pivot_a = pivot_temp;
        }
        dynamic_xor(zz_mat.Z[p-1][a], zz_mat.Z[p-1][b]);
        dynamic_xor(rep_mat.links[p-1][a], rep_mat.links[p-1][b]);
        if (both_boundary) {
            int chain_a = zz_mat.Cycle_Record[p-1][a].chain_idx;
            int chain_b = zz_mat.Cycle_Record[p-1][b].chain_idx;
            dynamic_xor(zz_mat.C[p][chain_a], zz_mat.C[p][chain_b]);
        }
        // Check for pivot conflicts again:
        zz_mat.pivots[p-1][pivot_b] = b;
        pivot_a = pivot(zz_mat.Z[p-1][a]);
        if (zz_mat.pivots[p-1].find(pivot_a) == zz_mat.pivots[p-1].end())
        {
            pivot_conflict = false;
            zz_mat.pivots[p-1][pivot_a] = a;
        }
        else 
        {
            if (zz_mat.pivots[p-1][pivot_a] == a) {
                pivot_conflict = false;
            }
            else {
                pivot_conflict = true;
                b = zz_mat.pivots[p-1][pivot_a];
                pivot_b = pivot_a;
            }
        }
    }
}
void reduce_bd_update(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    const int p,
    const int i,
    const int idx
)
{
    // Find all the chains in C[p] that contain the simplex and get the indices of their corresponding boundaries.
    vector<int> I;
    int alpha, chain_alpha; // alpha will be the only chain with negative birth timestamp such that C[p][chain_a] contains the simplex.
    int smallest_pivot = zz_mat.unique_id[p-1].size();
    for (auto cyc_idx: zz_mat.used_chain_Z[p-1])
    {
        if (!zz_mat.Cycle_Record[p-1][cyc_idx].non_bd) // If the birth timestamp is false, then the chain is a boundary.
        {
            int chain_idx = zz_mat.Cycle_Record[p-1][cyc_idx].chain_idx;
            if (idx < zz_mat.C[p][chain_idx] -> size() && (*zz_mat.C[p][chain_idx])[idx]==1) {
                I.push_back(cyc_idx);
                int pivot_cyc_idx = pivot(zz_mat.Z[p-1][cyc_idx]);
                if (pivot_cyc_idx < smallest_pivot) {
                    smallest_pivot = pivot_cyc_idx;
                    alpha = cyc_idx;
                    chain_alpha = chain_idx;
                }
            }
        }
    }
    // Add alpha to all the other chains in I.
    for (auto cyc_idx: I) {
        int chain_idx = zz_mat.Cycle_Record[p-1][cyc_idx].chain_idx;
        if (cyc_idx != alpha) {
            dynamic_xor(zz_mat.Z[p-1][cyc_idx], zz_mat.Z[p-1][alpha]);
            dynamic_xor(rep_mat.links[p-1][cyc_idx], rep_mat.links[p-1][alpha]);
            dynamic_xor(zz_mat.C[p][chain_idx], zz_mat.C[p][chain_alpha]);
        }
    }
    zz_mat.C[p][chain_alpha] = nullptr; // Zero out the chain C[p][chain_alpha] from C[p].
    // Update the used and available indices for C[p]:
    zz_mat.available_chain_C[p].push_back(chain_alpha);
    zz_mat.used_chain_C[p].erase(remove(zz_mat.used_chain_C[p].begin(), zz_mat.used_chain_C[p].end(), chain_alpha), zz_mat.used_chain_C[p].end()); // TODO: Check that this is done correctly.
    zz_mat.Cycle_Record[p-1][alpha] = cycle_record(true, -1, i+1);
    // Add a new wire to the (p-1)-th bundle:
    rep_mat.bundle[p-1].push_back(zz_mat.Z[p-1][alpha]);
    rep_mat.timestamp[p-1].push_back(i+1);
    chain new_link = make_shared<bitset>(rep_mat.bundle[p-1].size(), 0);
    new_link -> set(rep_mat.bundle[p-1].size()-1);
    rep_mat.links[p-1][alpha] = new_link;
}
uint delete_cycle(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    const vector<bool> &filt_op,
    const int p,
    const int idx
)
{
    vector<int> I; // Gather indices of chains that contain the simplex. 
    for (auto col_idx: zz_mat.used_chain_Z[p]) 
    {
        if (idx < (zz_mat.Z[p][col_idx] -> size()) && (*zz_mat.Z[p][col_idx])[idx] == 1)
        {
            I.push_back(col_idx);
        }
    }
    // sort I in the order of the birth timestamps where the order is the total order as above without using the sort function.
    sort(I.begin(), I.end(), [&](int &a, int &b){
        int birth_a = zz_mat.Cycle_Record[p][a].timestamp;
        int birth_b = zz_mat.Cycle_Record[p][b].timestamp;
            return ((birth_a == birth_b) || ((birth_a < birth_b) && (filt_op[birth_b-1])) || ((birth_a > birth_b) && (!(filt_op[birth_a-1]))));
            });
    // The chain to be deleted is the first chain in I.
    int alpha = I[0];
    chain z = zz_mat.Z[p][alpha];
    int z_pivot = pivot(z);
    zz_mat.pivots[p].erase(z_pivot); // Remove the pivot of Z[p][alpha] from pivots[p].
    I.erase(I.begin());
    // DELETE: chain alpha from Z[p] and links[p].
    int current_alpha = alpha; // Keeps track of the chain index being added in order to update links.
    if (I.size() != 0) 
    {
        // Iterate over the rest of I and ensure that the pivots remain distinct:
        for (auto a: I)
        {
            int a_pivot = pivot(zz_mat.Z[p][a]);
            if (a_pivot > z_pivot) // z does not change.
            {
                dynamic_xor(zz_mat.Z[p][a], z);
                dynamic_xor(rep_mat.links[p][a], rep_mat.links[p][current_alpha]);
            }   
            else 
            {
                // We want to deep copy Z[p][a] to a temporary chain and then copyZ[p][a] to z.
                chain temp = make_shared<bitset>(*zz_mat.Z[p][a]);
                int temp_pivot = a_pivot;
                dynamic_xor(zz_mat.Z[p][a], z);
                dynamic_xor(rep_mat.links[p][a], rep_mat.links[p][current_alpha]);
                a_pivot = z_pivot;
                // Update the information for z's.
                z = temp;
                current_alpha = a;
                z_pivot = temp_pivot;
                // Update the pivot of Z[p][a] in pivots[p].
                zz_mat.pivots[p][a_pivot] = a;
                zz_mat.pivots[p].erase(z_pivot);
            }
        }
    }
    return alpha;
}
void delete_update(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    const int p,
    const int alpha,
    const int idx
)
{
    // Remove the chains Z[p][alpha] from Z[p] and link[p][alpha] from links: assign these to null chains.
    zz_mat.Z[p][alpha] = nullptr;
    rep_mat.links[p][alpha] = nullptr;
    // Add this to the available chains for adding new cycles:
    zz_mat.available_chain_Z[p].push_back(alpha);
    zz_mat.used_chain_Z[p].erase(remove(zz_mat.used_chain_Z[p].begin(), zz_mat.used_chain_Z[p].end(), alpha), zz_mat.used_chain_Z[p].end());
    // We need to assign this to be invalid (redundant):
    zz_mat.Cycle_Record[p][alpha] = cycle_record(false, -1, -1); // TODO: Probably use assertions to test that boudaries are not using invalid chains.
}
void output_representatives(
    zigzag_matrices &zz_mat,
    representative_matrices &rep_mat,
    vector<tuple<int, int, int, vector<tuple<int, vector<int>>>>> *persistence,
    const int p,
    const int i,
    const int a,
    const int birth,
    const int death
)
{
    // Here, we need to gather all the wires pertaining to the cycle l and reconstruct the representatives at each index.
    vector<tuple<int, chain>> wire_representatives;
    vector<tuple<int, vector<int>>> wire_representatives_ids;
    // Alternatively, just go over the non-zero values using the next operator
    size_t col_idx = rep_mat.links[p][a] -> find_first();
    while (col_idx != rep_mat.links[p][a] -> npos) {
        wire_representatives.push_back(make_tuple(rep_mat.timestamp[p][col_idx], rep_mat.bundle[p][col_idx]));
        col_idx = rep_mat.links[p][a] -> find_next(col_idx);
    }
    // Sort the representatives in the order of the timestamps.
    sort(wire_representatives.begin(), wire_representatives.end(), [&](tuple<int, chain> &a, tuple<int, chain> &b){ return (get<0>(a) < get<0>(b));});
    // Initialize a bitset chain of size representative_max_size with all zeros.
    chain current_representative = make_shared<bitset>(1, 0);
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
            indices.push_back(zz_mat.unique_id[p][rep_idx]);
            rep_idx = current_representative -> find_next(rep_idx);
    }
    wire_representatives_ids.push_back(make_tuple(birth, indices));
    for (; wire_idx < wire_representatives.size(); ++wire_idx) {
        int time = get<0>(wire_representatives[wire_idx]);
        dynamic_xor(current_representative, get<1>(wire_representatives[wire_idx]));
        vector<int> indices;
        int rep_idx = current_representative -> find_first();
        while (rep_idx != current_representative -> npos) {
                indices.push_back(zz_mat.unique_id[p][rep_idx]);
                rep_idx = current_representative -> find_next(rep_idx);
        }
        wire_representatives_ids.push_back(make_tuple(time, indices));
    }
    persistence -> push_back(make_tuple(birth, death, p, wire_representatives_ids));
}





/* IMPLEMENTATION OF PRIMITIVE OPERATIONS: */
// Goes over a vector of ints and finds the position of the last non-zero element.
int pivot (chain a)
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
void dynamic_xor(chain a, chain b)
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
