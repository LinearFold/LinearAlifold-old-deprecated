
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <set>

#include "Linearalifold.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"
#include "Utils/ribo.h"

// #define SPECIAL_HP


using namespace std;


#ifdef lv
    bool comparefunc(std::pair<int,State> a, std::pair<int,State> b) {
        return a.first > b.first;
    }

    void BeamCKYParser::sort_keys(std::unordered_map<int, State> &map, std::vector<std::pair<int,State>> &sorted_keys) {
        sorted_keys.clear();
        for(auto &kv : map) {
            sorted_keys.push_back(kv);
        }
        sort(sorted_keys.begin(), sorted_keys.end(), comparefunc);    
    }
#endif

void BeamCKYParser::get_parentheses(char* result, string& seq) {
    memset(result, '.', seq_length);
    result[seq_length] = 0;

    stack<tuple<int, int, State>> stk;
    stk.push(make_tuple(0, seq_length-1, bestC[seq_length-1]));

    if(is_verbose){
        printf(">verbose\n");
    }
    // verbose stuff
    vector<pair<int,int>> multi_todo;
    unordered_map<int,int> mbp; // multi bp
    double total_energy = .0;
    double external_energy = .0;

    while ( !stk.empty() ) {
        tuple<int, int, State> top = stk.top();
        int i = get<0>(top), j = get<1>(top);
        State& state = get<2>(top);
        stk.pop();

        int k, p, q;

        switch (state.manner) {
            case MANNER_H:
                // this state should not be traced
                break;
            case MANNER_HAIRPIN:
                {
                    result[i] = '(';
                    result[j] = ')';
                }
                break;
            case MANNER_SINGLE:
                {
                    result[i] = '(';
                    result[j] = ')';
                    p = i + state.trace.paddings.l1;
                    q = j - state.trace.paddings.l2;
                    stk.push(make_tuple(p, q, bestP[q][p]));
                }
                break;
            case MANNER_HELIX:
                {
                    result[i] = '(';
                    result[j] = ')';
                    stk.push(make_tuple(i+1, j-1, bestP[j-1][i+1]));
                }
                break;
            case MANNER_MULTI: 
                p = i + state.trace.paddings.l1;
                q = j - state.trace.paddings.l2;
                stk.push(make_tuple(p, q, bestM2[q][p]));
                break;
            case MANNER_MULTI_eq_MULTI_plus_U:
                p = i + state.trace.paddings.l1;
                q = j - state.trace.paddings.l2;
                stk.push(make_tuple(p, q, bestM2[q][p]));
                break;
            case MANNER_P_eq_MULTI:
                result[i] = '(';
                result[j] = ')';
                stk.push(make_tuple(i, j, bestMulti[j][i]));
                break;
            case MANNER_M2_eq_M_plus_P:
                k = state.trace.split;
                stk.push(make_tuple(i, k, bestM[k][i]));
                stk.push(make_tuple(k+1, j, bestP[j][k+1]));
                break;
            case MANNER_M_eq_M2:
                stk.push(make_tuple(i, j, bestM2[j][i]));
                break;
            case MANNER_M_eq_M_plus_U:
                stk.push(make_tuple(i, j-1, bestM[j-1][i]));
                break;
            case MANNER_M_eq_P:
                stk.push(make_tuple(i, j, bestP[j][i]));
                break;
            case MANNER_C_eq_C_plus_U:
                k = j - 1;
                if (k != -1)
                    stk.push(make_tuple(0, k, bestC[k]));
                break;
            case MANNER_C_eq_C_plus_P:
                {
                    k = state.trace.split;
                    if (k != -1) {
                        stk.push(make_tuple(0, k, bestC[k]));
                        stk.push(make_tuple(k+1, j, bestP[j][k+1]));
                    }
                    else {
                        stk.push(make_tuple(i, j, bestP[j][i]));
                    }
                }
                break;
            default:  // MANNER_NONE or other cases
                printf("wrong manner at %d, %d: manner %d\n", i, j, state.manner); fflush(stdout);
                assert(false);
                
        }
    }

    return;
}

unsigned long quickselect_partition(vector<pair<value_type, int>>& scores, unsigned long lower, unsigned long upper) {
    value_type pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);

    }
    return upper;
}

// in-place quick-select
value_type quickselect(vector<pair<value_type, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}


value_type BeamCKYParser::beam_prune(std::unordered_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        value_type newscore;
        // lisiz: for _V, avoid -inf-int=+inf
        if ((k >= 0) && (bestC[k].score == VALUE_MIN)) newscore = VALUE_MIN;
        else newscore = (k >= 0 ? bestC[k].score : 0) + cand.score;
        scores.push_back(make_pair(newscore, i));
    }
    if (scores.size() <= beam) return VALUE_MIN;
    value_type threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}

void BeamCKYParser::sortM(value_type threshold,
                          std::unordered_map<int, State> &beamstep,
                          std::vector<std::pair<value_type, int>> &sorted_stepM) {
    sorted_stepM.clear();
    if (threshold == VALUE_MIN) {
        // no beam pruning before, so scores vector not usable
        for (auto &item : beamstep) {
            int i = item.first;
            State &cand = item.second;
            int k = i - 1;
            value_type newscore;
            // lisiz: constraints may cause all VALUE_MIN, sorting has no use
            if ((use_constraints) && (k >= 0) && (bestC[k].score == VALUE_MIN)) newscore = cand.score;
            else newscore = (k >= 0 ? bestC[k].score : 0) + cand.score;
            sorted_stepM.push_back(make_pair(newscore, i));
        }
    } else {
        for (auto &p : scores) {
            if (p.first >= threshold) sorted_stepM.push_back(p);
        }
    }

    sort(sorted_stepM.begin(), sorted_stepM.end(), std::greater<pair<value_type, int>>());
}

void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;

    bestH.clear();
    bestH.resize(seq_length);
    bestP.clear();
    bestP.resize(seq_length);
    bestM2.clear();
    bestM2.resize(seq_length);
    bestM.clear();
    bestM.resize(seq_length);
    bestC.clear();
    bestC.resize(seq_length);
    bestMulti.clear();
    bestMulti.resize(seq_length);

#ifdef is_cube_pruning
        sorted_bestM.clear();
        sorted_bestM.resize(seq_length);
#endif

    nucs.clear();
    nucs.resize(seq_length);

    scores.reserve(seq_length);
}

bool check_pairable_ij(vector<int> & SS_fast_i, vector<int> & SS_fast_j, float ** ribo, vector<vector<ribo_state>> & pscore, int i, int j){ //it is hc_decompose  = fc->hc->mx[n * i + j]; in mfe.c

    // in hard.c, RNAalifold use this to check if two pairs can be paired, 0 no, 63(VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS = 63) yes
    // if ((sn[i] != sn[j]) ||
    //       (((j - i + 1) <= md->max_bp_span) && ((j - i - 1) >= md->min_loop_size))) {
    //     printf("in hard.c line 772 (%d %d) (%d %d) (%d %d) (%d %d)\n", sn[i] , sn[j], i,j,j - i - 1,md->max_bp_span,j - i - 1,md->min_loop_size);
    //     int min_score = md->cv_fact * MINPSCORE;
    //     int act_score = (fc->hc->type == VRNA_HC_WINDOW) ?
    //                     fc->pscore_local[i][j - i] :
    //                     fc->pscore[fc->jindx[j] + i];


    //     if (act_score >= min_score)
    //       constraint = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;


    // md->cv_fact = 1.0
    // MINPSCORE = -2 * UNIT = -200

    int min_score = -200; //
    int act_score = pscore[i][j].ribo_score;
    if (act_score == VALUE_MIN) {
        act_score = make_pscores_ij(SS_fast_i, SS_fast_j, ribo);
        pscore[i][j].ribo_score = act_score;
    }

    if (act_score >= min_score) return true;

    return false;
}


vector<pair<int,int>> get_pairs(string& structure){

    vector<int> stack;
    vector<pair<int,int>> results;
    for (int i = 0 ; i < structure.size() ; i++){

        if (structure[i] == '('){
            stack.push_back(i);
        }

        else if (structure[i] == ')'){
            auto left = stack.back();
            stack.pop_back();
            results.push_back(make_pair(left,i));
        }

    }

    return results;
}

BeamCKYParser::DecoderResult BeamCKYParser::parse_alifold(std::vector<std::string> & MSA, float ** ribo, vector<vector<ribo_state>> & pscore, vector<vector<int>> & a2s_fast, vector<vector<int>> & s5_fast, vector<vector<int>> & s3_fast, vector<vector<int>> & SS_fast, vector<float> & smart_gap) {

    struct timeval parse_starttime, parse_endtime;

    int n_seq = MSA.size();

    // number of states
    unsigned long nos_H = 0, nos_P = 0, nos_M2 = 0,
            nos_M = 0, nos_C = 0, nos_Multi = 0;

    gettimeofday(&parse_starttime, NULL);

    prepare(static_cast<unsigned>(MSA[0].length()));

    vector<int> next_pair[NOTON];
    {
 
        for (int nuci = 0; nuci < NOTON; ++nuci) {
            next_pair[nuci].resize(MSA[0].size(), -1);
            int next = -1;
            for (int j = MSA[0].size()-1; j >=0; --j) {
                next_pair[nuci][j] = next;

                int temp_next = 1<<30;
                for (int s = 0; s < MSA.size(); s++){
                    // auto sequence_level_j = a2s[s][j];
                    if (_allowed_pairs[nuci][SS_fast[j][s]]){

                        temp_next = min(temp_next, j);
                    }

                }
                if (temp_next != 1<<30){
                    next = temp_next;
                }
            }
        }
    }

    vector<vector<int>> next_pair_MSA[MSA.size()];
    {
        for (int s = 0 ; s < MSA.size() ; s++){
            next_pair_MSA[s].resize(NOTON+1);
            for (int nuci = 0; nuci < NOTON+1; ++nuci) {
                next_pair_MSA[s][nuci].resize(seq_length, -1);
                int next = -1;
                for (int j = MSA[0].size()-1; j >=0; --j) {
                    next_pair_MSA[s][nuci][j] = next;
                    if (_allowed_pairs[nuci][SS_fast[j][s]]) next = j;
                }
            }
        }
    }

    auto next_position = init_next_position_only(seq_length);

    std::vector<std::string> seq_MSA_no_gap;
    seq_MSA_no_gap.resize(MSA.size());
    std::vector<std::vector<int>> if_tetraloops_MSA;
    std::vector<std::vector<int>> if_hexaloops_MSA;
    std::vector<std::vector<int>> if_triloops_MSA;
    if_tetraloops_MSA.resize(MSA.size());
    if_hexaloops_MSA.resize(MSA.size());
    if_triloops_MSA.resize(MSA.size());

#ifdef SPECIAL_HP
    for (int s = 0 ; s < MSA.size() ; s++){

        seq_MSA_no_gap[s] = "";

        for (auto nuc : MSA[s]){
            if (nuc != '-'){
                seq_MSA_no_gap[s] += nuc;
            }
        }
        seq_MSA_no_gap[s] += '\0';
    }

    v_init_tetra_hex_tri(seq_MSA_no_gap, if_tetraloops_MSA, if_hexaloops_MSA, if_triloops_MSA);
#endif

    // start CKY decoding

    if(seq_length > 0) bestC[0].set(- v_score_external_unpaired(0, 0), MANNER_C_eq_C_plus_U);
    if(seq_length > 1) bestC[1].set(- v_score_external_unpaired(0, 1), MANNER_C_eq_C_plus_U);

    ++nos_C;

    struct timeval starttime, endtime;

    gettimeofday(&starttime, NULL);

    float smart_gap_threshold = 0.5;
    // from left to right
    for(int j = 0; j < seq_length; ++j) {

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        auto & SS_j = SS_fast[j];
        auto & s3_j = s3_fast[j];
        auto & s5_j = s5_fast[j];
        auto & a2s_j = a2s_fast[j];
        auto & a2s_seq_length_1 = a2s_fast[seq_length-1];

        // beam of H
        {
            if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);
            
            
            if (smart_gap[j] - smart_gap[j-1] > smart_gap_threshold){

                int jnext = -1;

                if (next_position[j][j] != 0){
                    jnext = next_position[j][j];
                }

                else{

                    jnext = jnext_org;

                    for (int s = 0; s< MSA.size(); s++){
                        auto jnext_temp = next_pair_MSA[s][SS_j[s]][j];
                        while (jnext_temp -  j < 4 && jnext_temp != -1) jnext_temp = next_pair_MSA[s][SS_j[s]][jnext_temp];

                        if (jnext_temp != -1) jnext = min(jnext, jnext_temp);
                    }

                    if (jnext == jnext_org) jnext = -1;

                    next_position[j][j] = jnext;

                }


                while (1){
                    if (jnext == -1) break;
                    if (check_pairable_ij(SS_j,SS_fast[jnext], ribo, pscore, j, jnext)) break;

   
                    if (next_position[j][jnext] != 0){
                        jnext = next_position[j][jnext];
                    }


                    else{
                        int jnext_min = jnext_org;

                        auto prev_jnext = jnext;


                        for (int s = 0; s< MSA.size(); s++){
                            auto jnext_temp = next_pair_MSA[s][SS_j[s]][jnext];
                            if (jnext_temp != -1) jnext_min = min(jnext_min, jnext_temp);
                        }

                        if (jnext_min != jnext_org) jnext = jnext_min;
                        else jnext = -1;

                        next_position[j][prev_jnext] = jnext;

                    }

                }

                if (jnext != -1) {

                        value_type newscore = 0;

                        int tetra_hex_tri = -1;
                        auto & s5_jnext = s5_fast[jnext];
                        auto & SS_jnext = SS_fast[jnext];
                        auto & a2s_jnext_1 = a2s_fast[jnext - 1];
                        int new_nucj, new_nucj1, new_nucjnext_1, new_nucjnext, u;

                        for (int s = 0; s < n_seq; s++){
                            new_nucj = SS_j[s];
                            new_nucj1 = s3_j[s];
                            new_nucjnext_1 = s5_jnext[s];
                            new_nucjnext =  SS_jnext[s];
                            u = a2s_jnext_1[s] - a2s_j[s];

#ifdef SPECIAL_HP
                            if (a2s_j[s] >= 0 and a2s_j[s] < seq_MSA_no_gap[s].size()){
                                if (u == 4){ // 6:tetra
                                    tetra_hex_tri = if_tetraloops_MSA[s][a2s_j[s]];
                                }
                                else if (u == 6){ // 8:hexa
                                    tetra_hex_tri = if_hexaloops_MSA[s][a2s_j[s]];
                                }
                                else if (u == 3){ // 5:tri
                                    tetra_hex_tri = if_triloops_MSA[s][a2s_j[s]];
                                }
                            }
#endif


                            if (u < 3) newscore += -600;                            
                            else newscore += - v_score_hairpin(0, u + 1, new_nucj, new_nucj1, new_nucjnext_1, new_nucjnext, tetra_hex_tri);


                        }

                        if (pscore[j][jnext].ribo_score == VALUE_MIN){
                            pscore[j][jnext].ribo_score = make_pscores_ij(SS_j, SS_jnext, ribo);
                        }



                        update_if_better(bestH[jnext][j], newscore + pscore[j][jnext].ribo_score, MANNER_H);


                }
            }
            {
                // for every state h in H[j]
                //   1. extend h(i, j) to h(i, jnext)
                //   2. generate p(i, j)

                sort_keys(beamstepH, keys);
                for (auto &item : keys) {
                    int i = item.first;
                    auto & SS_i = SS_fast[i];
                    auto & s3_i = s3_fast[i];
                    auto & a2s_i = a2s_fast[i];

                    State &state = item.second;

                    // 2. generate p(i, j)
                    // lisiz, change the order because of the constriants
                    {
                        // add pscore when a P is generated.

                        // from H to P, must pairable
                        update_if_better(beamstepP[i], state.score , MANNER_HAIRPIN);

                        ++ nos_P;  
                    }

                    int jnext;

                    if (next_position[i][j] != 0){
                        jnext = next_position[i][j];
                    }

                    jnext = jnext_org;
                    for (int s = 0; s< n_seq; s++){
                        auto jnext_temp = next_pair_MSA[s][SS_i[s]][j];
                        if (jnext_temp != -1){
                            jnext = min(jnext, jnext_temp);
                        }
                    }

                    if (jnext == jnext_org){
                        jnext = -1;
                    }

                    next_position[i][j] = jnext;

                    while (1){

                        if (jnext == -1) break;
                        if (check_pairable_ij(SS_i,SS_fast[jnext], ribo, pscore, i, jnext)) break;

                        if (next_position[i][jnext] != 0){
                            jnext = next_position[i][jnext];
                        }

                        else{
                            auto prev_jnext = jnext;
                            int jnext_min = jnext_org;
                            for (int s = 0; s< MSA.size(); s++){
                                auto jnext_temp = next_pair_MSA[s][SS_i[s]][jnext];
                                if (jnext_temp != -1) jnext_min = min(jnext_min, jnext_temp);
                            }

                            if (jnext_min != jnext_org) jnext = jnext_min;
                            else jnext = -1;

                            next_position[i][prev_jnext] = jnext;

                        }
                    }

                    if (jnext != -1) {

                        // 1. extend h(i, j) to h(i, jnext)
                        value_type newscore = 0;

                        int tetra_hex_tri = -1;
                        auto & s5_jnext = s5_fast[jnext];
                        auto & SS_jnext = SS_fast[jnext];
                        auto & a2s_jnext_1 = a2s_fast[jnext - 1];

                        int new_nuci, new_nuci1, new_nucjnext_1, new_nucjnext, u;

                        for (int s = 0; s < n_seq; s++){
                            new_nuci = SS_i[s];
                            new_nuci1 = s3_i[s];
                            new_nucjnext_1 = s5_jnext[s];
                            new_nucjnext = SS_jnext[s];
                            u = a2s_jnext_1[s] - a2s_i[s];

#ifdef SPECIAL_HP
                            if (a2s_i[s] >= 0 and a2s_i[s] < seq_MSA_no_gap[s].size()){
                                if (u == 4){ // 6:tetra
                                    tetra_hex_tri = if_tetraloops_MSA[s][a2s_i[s]];
                                }
                                else if (u == 6){ // 8:hexa
                                    tetra_hex_tri = if_hexaloops_MSA[s][a2s_i[s]];
                                }
                                else if (u == 3){ // 5:tri
                                    tetra_hex_tri = if_triloops_MSA[s][a2s_i[s]];
                                }
                            }
#endif
                            if (u < 3) newscore += -600;
                            else newscore += - v_score_hairpin(0, u + 1, new_nuci, new_nuci1, new_nucjnext_1, new_nucjnext, tetra_hex_tri);
                        }

                        if (pscore[i][jnext].ribo_score == VALUE_MIN){
                            pscore[i][jnext].ribo_score = make_pscores_ij(SS_i, SS_jnext, ribo);
                        }
                        update_if_better(bestH[jnext][i], newscore + pscore[i][jnext].ribo_score, MANNER_H);

                    }
                }
            }
        }
        if (j == 0) continue;
        // beam of Multi
        {
            if (beam > 0 && beamstepMulti.size() > beam) beam_prune(beamstepMulti);
            // for every state in Multi[j]
            //   1. extend (i, j) to (i, jnext)
            //   2. generate P (i, j)
            sort_keys(beamstepMulti, keys);
            for (auto &item : keys) {
                int i = item.first;
                State& state = item.second;

                auto & SS_i = SS_fast[i];
                auto & s3_i = s3_fast[i];

                {
                    value_type newscore = 0;
                    newscore = state.score;
                    int new_nuci, new_nuci1, new_nucj_1, new_nucj;

                    for (int s = 0; s < n_seq; s++){
                        new_nuci = SS_i[s];
                        new_nuci1 = s3_i[s];
                        new_nucj_1 = s5_j[s];
                        new_nucj = SS_j[s];
                        newscore += - v_score_multi(-1, -1, new_nuci, new_nuci1, new_nucj_1, new_nucj, -1);
                    }


                    if (pscore[i][j].ribo_score == VALUE_MIN){
                        pscore[i][j].ribo_score = make_pscores_ij(SS_fast[i], SS_fast[j], ribo);
                    }

                    newscore += pscore[i][j].ribo_score;

                    update_if_better(beamstepP[i], newscore, MANNER_P_eq_MULTI);

                }

                {
                    int jnext;

                    if (next_position[i][j] != 0){
                        jnext = next_position[i][j];
                    }

                    else{
                        jnext = jnext_org;
                        for (int s = 0; s< MSA.size();s++){
                            auto jnext_temp = next_pair_MSA[s][SS_i[s]][j];
                            if (jnext_temp != -1){
                                jnext = min(jnext, jnext_temp);
                            }
                        }

                        if (jnext == jnext_org){
                            jnext = -1;
                        }
                        next_position[i][j] = jnext;
                    }

                    while (1){

                        if (jnext == -1) break;
                        if (check_pairable_ij(SS_i,SS_fast[jnext], ribo, pscore, i, jnext)) break;

                        if (next_position[i][jnext] != 0){
                            jnext = next_position[i][jnext];
                        }

                        else{
                            auto prev_jnext = jnext;
                            int jnext_min = jnext_org;

                            for (int s = 0; s< MSA.size(); s++){
                                auto jnext_temp = next_pair_MSA[s][SS_i[s]][jnext];
                                if (jnext_temp != -1) jnext_min = min(jnext_min, jnext_temp);
                            }

                            if (jnext_min != jnext_org) jnext = jnext_min;
                            else jnext = -1;

                             next_position[i][prev_jnext] = jnext;
                        }
                    }

                    if (jnext != -1) {
                        int new_l1 = state.trace.paddings.l1;
                        int new_l2 = state.trace.paddings.l2 + jnext - j;
                        // 1. extend (i, j) to (i, jnext)
                        value_type newscore = 0;
                        newscore = state.score;
                        update_if_better(bestMulti[jnext][i], newscore, MANNER_MULTI_eq_MULTI_plus_U, new_l1, new_l2);
                    }
                }
            }
        }
        // beam of P
        {
            if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

            // for every state in P[j]
            //   1. generate new helix/bulge
            //   2. M = P
            //   3. M2 = M + P
            //   4. C = C + P
#ifdef is_cube_pruning
            bool use_cube_pruning = beam > MIN_CUBE_PRUNING_SIZE
                                    && beamstepP.size() > MIN_CUBE_PRUNING_SIZE;
#else
            bool use_cube_pruning = false;
#endif               

            sort_keys(beamstepP, keys);
            for (auto &item : keys) {
                int i = item.first;
                State& state = item.second;

                auto & s5_i = s5_fast[i];
                auto & SS_i = SS_fast[i];
                auto & a2s_i = a2s_fast[i];
                auto & a2s_j = a2s_fast[j];

                // 2. M = P
                if(i > 0 && j < seq_length-1){
                    value_type newscore;
                    newscore = state.score;
                    int new_nuci_1, new_nuci, new_nucj, new_nucj1;

                    for (int s = 0; s < n_seq; s++){

                        new_nuci_1 = ((i - 1) > -1)? s5_i[s] : -1;
                        new_nuci = SS_i[s];
                        new_nucj = SS_j[s];
                        new_nucj1 = (j + 1) < seq_length? s3_j[s] : -1;
                        newscore += - v_score_M1(-1, -1, -1, new_nuci_1, new_nuci, new_nucj, new_nucj1, -1); // no position information needed
                    }
                    update_if_better(beamstepM[i], newscore, MANNER_M_eq_P);
                }
                // 3. M2 = M + P
                if(!use_cube_pruning) {
                    int k = i - 1;
                    if ( k > 0 && !bestM[k].empty()) {
                        value_type M1_score;

                        M1_score = state.score;

                        int new_nuci_1, new_nuci, new_nucj, new_nucj1;

                        for (int s = 0; s < n_seq; s++){

                            new_nuci_1 = s5_i[s]; //TODO, need to check boundary?
                            new_nuci = SS_i[s];
                            new_nucj = SS_j[s];
                            new_nucj1 = (j + 1) < seq_length ? s3_j[s] : -1; //TODO, need to check boundary? it may diff. from RNAalifold
                            M1_score += - v_score_M1(-1, -1, -1, new_nuci_1, new_nuci, new_nucj, new_nucj1, -1);
                        }
                        // candidate list
                        auto bestM2_iter = beamstepM2.find(i);
#ifndef is_candidate_list
                        for (auto &m : bestM[k]) {
                            int newi = m.first;
                            // eq. to first convert P to M1, then M2/M = M + M1
                            value_type newscore = M1_score + m.second.score;
                            update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                        }
#else
                        if (bestM2_iter==beamstepM2.end() || M1_score > bestM2_iter->second.score) {
                            for (auto &m : bestM[k]) {
                                int newi = m.first;
                                // eq. to first convert P to M1, then M2/M = M + M1
                                value_type newscore = M1_score + m.second.score;
                                update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                            }
                        }
#endif
                    }
                }
                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                      State& prefix_C = bestC[k];
                      if (prefix_C.manner != MANNER_NONE) {
                        value_type newscore;
                        newscore = prefix_C.score + state.score;
                        int new_nuck, new_nuck1, new_nucj, new_nucj1;

                        for (int s = 0; s < n_seq; s++){

                            new_nuck = (a2s_i[s] > 0) ? s5_i[s] : -1; //external.c line 1165, weird
                            new_nuck1 = SS_i[s];
                            new_nucj = SS_j[s];
                            new_nucj1 = (a2s_j[s] < a2s_seq_length_1[s]) ? s3_j[s] : -1; //external.c line 1165, weird
                            newscore += - v_score_external_paired(-1, -1, new_nuck, new_nuck1, new_nucj, new_nucj1, -1);
                        }

                        update_if_better(beamstepC, newscore, MANNER_C_eq_C_plus_P, k);
                      }
                    } else {
                        value_type newscore;
                        newscore = state.score;

                        int new_nuck1, new_nucj, new_nucj1;

                        for (int s = 0; s < MSA.size(); s++){

                            new_nuck1 = SS_i[s];
                            new_nucj = SS_j[s];
                            new_nucj1 = (a2s_j[s] < a2s_seq_length_1[s]) ? s3_j[s] : -1; //external.c line 1165, weird
                            newscore += - v_score_external_paired(0, j, -1, new_nuck1, new_nucj, new_nucj1, -1);
                        }
                        update_if_better(beamstepC, newscore, MANNER_C_eq_C_plus_P, -1);
                    }
                }
                // 1. P2P generate new helix / single_branch
                // new state is of shape p..i..j..q
                if (i >0 && j<seq_length-1) {

                    auto & a2s_i_1 = a2s_fast[i-1];

                    int *tt2;
                    tt2 = (int *)vrna_alloc(sizeof(int) * n_seq);
                    auto & SS_i = SS_fast[i];
                    auto & SS_j = SS_fast[j];
                    for (int s = 0; s < n_seq; s++){
                        
                        tt2[s] = NUM_TO_PAIR(SS_j[s], SS_i[s]);
                    }

                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {

                        auto & SS_p = SS_fast[p];
                        auto & s5_p = s5_fast[p];
                        auto & s3_p = s3_fast[p];   
                        auto & a2s_p = a2s_fast[p];                    

                        int q;

                        if (next_position[p][j] != 0){
                            q = next_position[p][j];
                        }

                        else {

                            q = jnext_org;
                            
                            for (int s = 0; s < MSA.size(); s++){
                                auto jnext_temp = next_pair_MSA[s][SS_p[s]][j];
                                if (jnext_temp != -1){
                                    q = min(q, jnext_temp);
                                }
                            }

                            if (q == jnext_org){
                                q = -1;
                            }

                            next_position[p][j] = q;

                        }

                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {

                            auto & SS_q = SS_fast[q];
                            auto & s5_q = s5_fast[q];
                            auto & s3_q = s3_fast[q];

                            if (check_pairable_ij(SS_p,SS_q, ribo, pscore,p,q)){

                                if (p == i - 1 && q == j + 1) {
                                    int nucp, nucq, nucp1, nucq_1, nuci_1, nucj1, u1_local, u2_local, type;
                                    value_type newscore = state.score;

                                    for (int s = 0; s < n_seq; s++){

                                        nucp = SS_p[s];
                                        nucq = SS_q[s];
                                        nucp1 = s3_p[s];
                                        nucq_1 = s5_q[s];
                                        nuci_1 = s5_i[s];
                                        nucj1 = s3_j[s];

                                        type = NUM_TO_PAIR(nucp, nucq);
                                        newscore += -v_score_single_alifold(0, 0, type, tt2[s], nucp1, nucq_1, nuci_1, nucj1); //internal.c line 476, left, right gaps are all 0
                                    }

                                    if (pscore[p][q].ribo_score == VALUE_MIN){
                                        pscore[p][q].ribo_score = make_pscores_ij(SS_fast[p], SS_fast[q], ribo);
                                    }

                                    newscore += pscore[p][q].ribo_score;
                                    update_if_better(bestP[q][p], newscore, MANNER_HELIX);

                                } else {
                                    auto & a2s_q_1 = a2s_fast[q-1];
                                    int nucp, nucq, nucp1, nucq_1, nuci_1, nucj1, u1_local, u2_local, type;
                                    value_type newscore = state.score;
                                    for (int s = 0; s < n_seq; s++){
                                        nucp = SS_p[s];
                                        nucq = SS_q[s];
                                        nucp1 = s3_p[s];
                                        nucq_1 = s5_q[s];
                                        nuci_1 = s5_i[s];
                                        nucj1 = s3_j[s];

                                        // in internal.c
                                        // i    k     l   j RNAalifold
                                        // p    i     j   q ours

                                        u1_local  = a2s_i_1[s] - a2s_p[s];
                                        u2_local  = a2s_q_1[s] - a2s_j[s];

                                        type = NUM_TO_PAIR(nucp, nucq);
                                        newscore += -v_score_single_alifold(u1_local, u2_local, type, tt2[s], nucp1, nucq_1, nuci_1, nucj1); 



                                    }

                                    if (pscore[p][q].ribo_score == VALUE_MIN){
                                        pscore[p][q].ribo_score = make_pscores_ij(SS_fast[p], SS_fast[q], ribo);
                                    }
                                    newscore += pscore[p][q].ribo_score;
                                    update_if_better(bestP[q][p], newscore, MANNER_SINGLE, (i - p), q - j);

                                }

                            }

                            if (next_position[p][q] != 0){
                                q = next_position[p][q];
                            }

                            else{

                                auto prev_q = q;

                                int temp_q = jnext_org;
                                
                                for (int s = 0; s < MSA.size(); s++){
                                    auto jnext_temp = next_pair_MSA[s][SS_p[s]][q];
                                    if (jnext_temp != -1){
                                        temp_q = min(temp_q, jnext_temp);
                                    }
                                }

                                if (temp_q == jnext_org){
                                    q = -1;
                                }
                                else q = temp_q;

                                next_position[p][prev_q] = q;
                            }

                        }
                    }
                    free(tt2);
                }

            }
        }
        // beam of M2
        {
            if (beam > 0 && beamstepM2.size() > beam) beam_prune(beamstepM2);

            // for every state in M2[j]
            //   1. multi-loop  (by extending M2 on the left)
            //   2. M = M2
            sort_keys(beamstepM2, keys);
            for (auto &item : keys) {

                int i = item.first;
                State& state = item.second;

                // 2. M = M2
                {
                    update_if_better(beamstepM[i], state.score, MANNER_M_eq_M2);
                    ++ nos_M;
                }

                // 1. multi-loop
                {
                    auto smart_i = smart_gap[i-1];
                    for (int p = i-1; ((smart_i -  smart_gap[p]) <= 2*SINGLE_MAX_LEN) and p >= 0; --p) {

                        if (smart_gap[p] - smart_gap[p-1] < smart_gap_threshold){

                            continue;
                        }

                        auto & SS_p = SS_fast[p];


                        int q;

                        if (next_position[p][j] != 0){
                            q = next_position[p][j];
                        }

                        else{
                            q = jnext_org;
                            for (int s = 0; s< MSA.size();s++){
                                auto jnext_temp = next_pair_MSA[s][SS_p[s]][j];

                                if (jnext_temp != -1){
                                    q = min(q, jnext_temp);
                                }
                            }

                            if (q == jnext_org){
                                q = -1;
                            }
                            next_position[p][j] = q;
                        }

                        while (1){

                            if (q == -1) break;
                            if (check_pairable_ij(SS_p, SS_fast[q], ribo, pscore,p,q)) break;


                            if (next_position[p][q] != 0){
                                q = next_position[p][q];
                            }

                            else{
                                auto prev_q = q;
                                int jnext_min = jnext_org;

                                for (int s = 0; s< MSA.size(); s++){
                                    auto jnext_temp = next_pair_MSA[s][SS_p[s]][q];
                                    if (jnext_temp != -1) jnext_min = min(jnext_min, jnext_temp);
                                }

                                if (jnext_min != jnext_org) q = jnext_min;
                                else q = -1;

                                next_position[p][prev_q] = q;

                            }
                        }

                        if (q != -1) {
                            // the current shape is p..i M2 j ..q
                            value_type newscore;

                            newscore = state.score;

                            update_if_better(bestMulti[q][p], newscore, MANNER_MULTI, (i - p), q - j);

                        }
                    }
                }
            }
        }
        // beam of M
        {
            value_type threshold = VALUE_MIN;
            if (beam > 0 && beamstepM.size() > beam) threshold = beam_prune(beamstepM);

#ifdef is_cube_pruning
                sortM(threshold, beamstepM, sorted_bestM[j]);
#endif

            // for every state in M[j]
            //   1. M = M + unpaired
            sort_keys(beamstepM, keys);
            for (auto &item : keys) {

                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
                    value_type newscore;
                    newscore = state.score;
                    update_if_better(bestM[j+1][i], newscore, MANNER_M_eq_M_plus_U);
                }
            }
        }

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {

                value_type newscore;
                newscore = beamstepC.score;
                update_if_better(bestC[j+1], newscore, MANNER_C_eq_C_plus_U);
            }
        }

    }  // end of for-loo j

    State& viterbi = bestC[seq_length-1];
    char result[seq_length+1];
    get_parentheses(result, MSA[0]);
    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
    unsigned long nos_tot = nos_H + nos_P + nos_M2 + nos_Multi + nos_M + nos_C;
    fflush(stdout);
    return {string(result), viterbi.score, nos_tot, parse_elapsed_time};
}


BeamCKYParser::BeamCKYParser(int beam_size,
                             bool nosharpturn,
                             bool verbose)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose){
      initialize();
}


int main(int argc, char** argv){


    struct timeval parse_alifold_starttime, parse_alifold_endtime;

    gettimeofday(&parse_alifold_starttime, NULL);

    int beamsize = 100;
    bool sharpturn = false;
    bool is_verbose = false;


    if (argc >= 1) {
        beamsize = atoi(argv[1]);
        is_verbose = atoi(argv[2]) == 1;
    }

    std::vector<std::string> MSA;


    for (string seq; getline(cin, seq);) {
        // printf("Input: %s\n", seq.c_str());
        if (seq.length() == 0)
            continue;

        if (seq[0] == ';' || seq[0] == '>') {
            // printf("%s\n", seq.c_str());
            continue;
        }
        
        // convert to uppercase
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

        // convert T to U
        replace(seq.begin(), seq.end(), 'T', 'U');

        MSA.push_back(seq);

        // printf("%s\n", seq.c_str());
    }

    auto n_seq = MSA.size();
    auto MSA_seq_length = MSA[0].size();


    auto ribo = get_ribosum(MSA, n_seq, MSA_seq_length);
    auto pscore = init_pscores_only(MSA_seq_length);
    vector<float> smart_gap;
    vector<vector<int>> a2s_fast, s5_fast, s3_fast, SS_fast;
    a2s_prepare_is(MSA, n_seq, MSA_seq_length, a2s_fast, s5_fast, s3_fast, SS_fast, smart_gap);
    BeamCKYParser parser(beamsize, !sharpturn, is_verbose);
    BeamCKYParser::DecoderResult result_alifold = parser.parse_alifold(MSA, ribo, pscore, a2s_fast, s5_fast, s3_fast, SS_fast, smart_gap);
    gettimeofday(&parse_alifold_endtime, NULL);
    double parse_elapsed_time = parse_alifold_endtime.tv_sec - parse_alifold_starttime.tv_sec + (parse_alifold_endtime.tv_usec-parse_alifold_starttime.tv_usec)/1000000.0;

    double printscore = (result_alifold.score / -100.0);
    auto result_pairs = get_pairs(result_alifold.structure);
    float pscore_f = 0.;
    for (auto & pair : result_pairs){
        pscore_f += pscore[pair.first][pair.second].ribo_score;
    }
    pscore_f = -pscore_f/n_seq/100.;

    printf("%s (%.2f = %.2f + %.2f)\n", result_alifold.structure.c_str(), printscore/n_seq, printscore/n_seq - pscore_f, pscore_f);
    if (is_verbose){
        printf("beam size %d\n", beamsize);
        printf("runtime %.2f seconds\n", parse_elapsed_time);
    }
    return 0;
}
