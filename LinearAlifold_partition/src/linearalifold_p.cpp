
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
#include <stdio.h> 

#include "linearalifold_p.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"
#include "bpp.cpp"
// #include "Utils/ribo.h"

#define SPECIAL_HP

using namespace std;

unsigned long quickselect_partition(vector<pair<pf_type, int>>& scores, unsigned long lower, unsigned long upper) {
    pf_type pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);
    }
    return upper;
}

// in-place quick-select
pf_type quickselect(vector<pair<pf_type, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}


pf_type BeamCKYParser::beam_prune(std::unordered_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        pf_type newalpha = (k >= 0 ? bestC[k].alpha : pf_type(0.0)) + cand.alpha;
        scores.push_back(make_pair(newalpha, i));
    }
    if (scores.size() <= beam) return VALUE_MIN;
    pf_type threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}


void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;

    nucs = new int[seq_length];
    bestC = new State[seq_length];
    bestH = new unordered_map<int, State>[seq_length];
    bestP = new unordered_map<int, State>[seq_length];
    bestM = new unordered_map<int, State>[seq_length];
    bestM2 = new unordered_map<int, State>[seq_length];
    bestMulti = new unordered_map<int, State>[seq_length];
    
    scores.reserve(seq_length);
}

void BeamCKYParser::postprocess() {

    delete[] bestC;  
    delete[] bestH;  
    delete[] bestP;  
    delete[] bestM;  
    delete[] bestM2;  
    delete[] bestMulti;  

    delete[] nucs;  
}


bool BeamCKYParser::check_pairable_ij(vector<int> & SS_fast_i, vector<int> & SS_fast_j, float ** ribo, vector<vector<ribo_state>> & pscore, int i, int j){ //it is hc_decompose  = fc->hc->mx[n * i + j]; in mfe.c

    // in hard.c, RNAalifold use this to check if column i and column j can be paired, 0 no, 63(VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS = 63) yes
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


    if (act_score == std::numeric_limits<int>::lowest()) {
        act_score = make_pscores_ij(SS_fast_i, SS_fast_j, ribo);
        pscore[i][j].ribo_score = act_score;
    }

    if (act_score >= min_score) return true;

    return false;
    // return true;

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


void BeamCKYParser::parse_alifold(std::vector<std::string> & MSA_, vector<vector<int>> & a2s_, vector<vector<ribo_state>> & pscore, vector<vector<int>> & s5_, vector<vector<int>> & s3_, vector<vector<int>> & SS_, float ** ribo_, vector<float> & smart_gap_) {
    
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    auto seq = MSA_[0];

    prepare(static_cast<unsigned>(seq.length()));

    MSA = MSA_;
    a2s_fast = a2s_;
    s5_fast = s5_;
    s3_fast = s3_;
    SS_fast = SS_;
    ribo = ribo_;
    smart_gap = smart_gap_;

    int n_seq = MSA.size();

    double kTn = double(kT) * n_seq;


    nucs_MSA.clear();
    nucs_MSA.resize(MSA.size());


    for (int s = 0; s < MSA.size(); s++){
        nucs_MSA[s].resize(seq_length);
        for (int i = 0; i < seq_length; ++i){
            nucs_MSA[s][i] = GET_ACGU_NUM(MSA[s][i]);
        }
    }


    vector<vector<vector<int>>> next_pair_MSA_;

    next_pair_MSA_.resize(MSA.size());


    for (int s = 0 ; s < MSA.size() ; s++){
        next_pair_MSA_[s].resize(NOTON+1);
        for (int nuci = 0; nuci < NOTON+1; ++nuci) {
            next_pair_MSA_[s][nuci].resize(seq_length, -1);
            int next = -1;
            for (int j = MSA[0].size()-1; j >=0; --j) {
                next_pair_MSA_[s][nuci][j] = next;
                if (_allowed_pairs[nuci][nucs_MSA[s][j]]) next = j;
            }
        }
    }

    next_pair_MSA = next_pair_MSA_;

    vector<vector<int>>().swap(next_pair_ij);


    auto next_position = init_next_position_only(seq_length);



    for (int i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    vector<int> next_pair[NOTON];
    {
        for (int nuci = 0; nuci < NOTON; ++nuci) {
            // next_pair
            next_pair[nuci].resize(seq_length, -1);
            int next = -1;
            for (int j = seq_length-1; j >=0; --j) {
                next_pair[nuci][j] = next;
                if (_allowed_pairs[nuci][nucs[j]]) next = j;
            }
        }
    }



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

    v_init_tetra_hex_tri(MSA, if_tetraloops_MSA, if_hexaloops_MSA, if_triloops_MSA);
#endif



        if(seq_length > 0) bestC[0].alpha = 0.0;
        if(seq_length > 1) bestC[1].alpha = 0.0;

    float smart_gap_threshold = 0.5;

    value_type newscore;


    struct timeval starttime, endtime;

    gettimeofday(&starttime, NULL);


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

                    int tetra_hex_tri = -1;
                    auto & s5_jnext = s5_fast[jnext];
                    auto & SS_jnext = SS_fast[jnext];
                    auto & a2s_jnext_1 = a2s_fast[jnext - 1];
                    int new_nucj, new_nucj1, new_nucjnext_1, new_nucjnext, u;

                    newscore = 0;
                    for (int s = 0; s < MSA.size(); s++){

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

                        if (u < 3) {
                            newscore += -600;
                        }
                        else newscore += - v_score_hairpin(0, u + 1, new_nucj, new_nucj1, new_nucjnext_1, new_nucjnext, tetra_hex_tri);

                    }

                    Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore/kTn);
                }
            }

            {
                // for every state h in H[j]
                //   1. extend h(i, j) to h(i, jnext)
                //   2. generate p(i, j)
                for (auto &item : beamstepH) {
                    int i = item.first;
                    State &state = item.second;
                    auto & SS_i = SS_fast[i];
                    auto & s3_i = s3_fast[i];
                    auto & a2s_i = a2s_fast[i];

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
                        if (check_pairable_ij(SS_i, SS_fast[jnext], ribo, pscore, i, jnext)) break;

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

                        int tetra_hex_tri = -1;
                        newscore = 0;
                        auto & s5_jnext = s5_fast[jnext];
                        auto & SS_jnext = SS_fast[jnext];
                        auto & a2s_jnext_1 = a2s_fast[jnext - 1];
                        int new_nuci, new_nuci1, new_nucjnext_1, new_nucjnext, u;

                        for (int s = 0; s < MSA.size(); s++){

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

                        Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore/kTn);

                    }

                    // 2. generate p(i, j)
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha);

                }
            }
        }
        if (j == 0) continue;

        // beam of Multi
        {
            if (beam > 0 && beamstepMulti.size() > beam) beam_prune(beamstepMulti);

            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;

                auto & SS_i = SS_fast[i];
                auto & s3_i = s3_fast[i];


                // 1. extend (i, j) to (i, jnext)
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
                        Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha);
                    }
                }

                // 2. generate P (i, j)
                {
                    int new_nuci, new_nuci1, new_nucj_1, new_nucj;

                    newscore = 0;
                    for (int s = 0; s < MSA.size(); s++){

                        new_nuci = SS_i[s];
                        new_nuci1 = s3_i[s];
                        new_nucj_1 = s5_j[s];
                        new_nucj = SS_j[s];

                        newscore += (- v_score_multi(-1, -1, new_nuci, new_nuci1, new_nucj_1, new_nucj, -1));

                    }

                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha + newscore / kTn);

                }
            }
        }

        // beam of P
        {   

            for(auto& item : beamstepP) { // can we merge it into an existing DP?
                int i = item.first;
                State& state = item.second;

                if (pscore[i][j].ribo_score == std::numeric_limits<int>::lowest()){
                    pscore[i][j].ribo_score = make_pscores_ij(SS_fast[i], SS_fast[j], ribo);
                }

                state.alpha = state.alpha + pscore[i][j].ribo_score / (kTn);

            }

            if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

            // for every state in P[j]
            //   1. generate new helix/bulge
            //   2. M = P
            //   3. M2 = M + P
            //   4. C = C + P
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;

                auto & s5_i = s5_fast[i];
                auto & SS_i = SS_fast[i];

                auto & a2s_i = a2s_fast[i];
                auto & a2s_j = a2s_fast[j];

                // 1. generate new helix / single_branch
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

                                    newscore = 0;

                                    for (int s = 0; s < MSA.size(); s++){

                                        nucp = SS_p[s];
                                        nucq = SS_q[s];

                                        nucp1 = s3_p[s];
                                        nucq_1 = s5_q[s];
                                        nuci_1 = s5_i[s];
                                        nucj1 = s3_j[s];

                                        type = NUM_TO_PAIR(nucp, nucq);


                                        newscore += -v_score_single_alifold(0, 0, type, tt2[s], nucp1, nucq_1, nuci_1, nucj1); //internal.c line 476, left, right gaps are all 0

                                    }

                                    Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore / kTn);

                                } else {
                                    // single branch

                                    auto & a2s_q_1 = a2s_fast[q-1];
                                    int nucp, nucq, nucp1, nucq_1, nuci_1, nucj1, u1_local, u2_local, type;

                                    newscore = 0;
                                    for (int s = 0; s < MSA.size(); s++){

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


                                    Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore / kTn);

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

                // 2. M = P
                if(i > 0 && j < seq_length-1){

                    newscore = 0;
                    int new_nuci_1, new_nuci, new_nucj, new_nucj1;
                    auto & s5_i = s5_fast[i];
                    auto & SS_i = SS_fast[i];
                    auto & SS_j = SS_fast[j];
                    auto & s3_j = s3_fast[j];

                    for (int s = 0; s < MSA.size(); s++){

                        new_nuci_1 = ((i - 1) > -1)? s5_i[s] : -1;
                        new_nuci = SS_i[s];
                        new_nucj = SS_j[s];
                        new_nucj1 = (j + 1) < seq_length? s3_j[s] : -1;
                        newscore += (- v_score_M1(-1, -1, -1, new_nuci_1, new_nuci, new_nucj, new_nucj1, -1)); // no position information needed
                    }
                        Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha + newscore/kTn);

                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {

                    newscore = 0;

                    int new_nuci_1, new_nuci, new_nucj, new_nucj1;

                    for (int s = 0; s < MSA.size(); s++){

                        new_nuci_1 = s5_i[s]; //TODO, need to check boundary?
                        new_nuci = SS_i[s];
                        new_nucj = SS_j[s];
                        new_nucj1 = (j + 1) < seq_length ? s3_j[s] : -1; //TODO, need to check boundary? it may diff. from RNAalifold

                        newscore += (- v_score_M1(-1, -1, -1, new_nuci_1, new_nuci, new_nucj, new_nucj1, -1));
                    }

                    pf_type m1_alpha = state.alpha + newscore / kTn;
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(beamstepM2[newi].alpha, m_state.alpha + m1_alpha);
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        State& prefix_C = bestC[k];

                        newscore = 0;

                        int new_nuck, new_nuck1, new_nucj, new_nucj1;

                        for (int s = 0; s < MSA.size(); s++){

                            new_nuck = (a2s_i[s] > 0) ? s5_i[s] : -1; //external.c line 1165, weird
                            new_nuck1 = SS_i[s];
                            new_nucj = SS_j[s];
                            new_nucj1 = (a2s_j[s] < a2s_seq_length_1[s]) ? s3_j[s] : -1; //external.c line 1165, weird

                            newscore += (- v_score_external_paired(-1, -1, new_nuck, new_nuck1, new_nucj, new_nucj1, -1));
                        }

                        Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + newscore/kTn);

                    } else {

                        newscore = 0;
                        int new_nuck1, new_nucj, new_nucj1;
                        for (int s = 0; s < MSA.size(); s++){

                            new_nuck1 = SS_i[s];
                            new_nucj = SS_j[s];
                            new_nucj1 = (a2s_j[s] < a2s_seq_length_1[s]) ? s3_j[s] : -1; //external.c line 1165, weird
                            newscore += - v_score_external_paired(0, j, -1, new_nuck1, new_nucj, new_nucj1, -1);
                        }
                        Fast_LogPlusEquals(beamstepC.alpha, state.alpha + newscore/kTn);

                    }
                }
            }
        }

        // beam of M2
        {
            if (beam > 0 && beamstepM2.size() > beam) beam_prune(beamstepM2);

            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

                // 1. multi-loop

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
                        Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha);      
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha);  
            }
        }

        // beam of M
        {
            if (beam > 0 && beamstepM.size() > beam) beam_prune(beamstepM);

            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
                    Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha); 
                }
            }
        }

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
                Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha);    
            }
        }
    }  // end of for-loo j


    State& viterbi = bestC[seq_length-1];

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    fprintf(stdout,"Free Energy of Ensemble: %.2f kcal/mol\n", -kTn * viterbi.alpha / 100.0 / MSA.size());
    if(is_verbose) fprintf(stdout,"Partition Function Calculation Time: %.2f seconds.\n", parse_elapsed_time);
    fflush(stdout);

    // lhuang
    if(pf_only && !forest_file.empty()) dump_forest(seq, true); // inside-only forest

    if(!pf_only){

        outside_alifold(next_pair, pscore, smart_gap, smart_gap_threshold, next_position);

        if (!forest_file.empty())
          dump_forest(seq, false); // inside-outside forest
            cal_PairProb(viterbi, pscore);

        if (mea_) {
            PairProb_MEA(seq);
        }

        if (threshknot_){
            ThreshKnot(seq);
        }
    }
    postprocess();
    return;
}


void BeamCKYParser::print_states(FILE *fptr, unordered_map<int, State>& states, int j, string label, bool inside_only, double threshold) {    
    for (auto & item : states) {
        int i = item.first;
        State & state = item.second;
        if (inside_only) fprintf(fptr, "%s %d %d %.5lf\n", label.c_str(), i+1, j+1, state.alpha);
        else if (state.alpha + state.beta > threshold) // lhuang : alpha + beta - totalZ < ...
            fprintf(fptr, "%s %d %d %.5lf %.5lf\n", label.c_str(), i+1, j+1, state.alpha, state.beta);
    }
}

void BeamCKYParser::dump_forest(string seq, bool inside_only) {  
    printf("Dumping (%s) Forest to %s...\n", (inside_only ? "Inside-Only" : "Inside-Outside"), forest_file.c_str());
    FILE *fptr = fopen(forest_file.c_str(), "w");  // lhuang: should be fout >>
    fprintf(fptr, "%s\n", seq.c_str());
    int n = seq.length(), j;
    for (j = 0; j < n; j++) {
        if (inside_only) fprintf(fptr, "E %d %.5lf\n", j+1, bestC[j].alpha);
        else fprintf(fptr, "E %d %.5lf %.5lf\n", j+1, bestC[j].alpha, bestC[j].beta);
    }
    double threshold = bestC[n-1].alpha - 9.91152; // lhuang -9.xxx or ?
    for (j = 0; j < n; j++) 
        print_states(fptr, bestP[j], j, "P", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM[j], j, "M", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM2[j], j, "M2", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestMulti[j], j, "Multi", inside_only, threshold);
}

BeamCKYParser::BeamCKYParser(int beam_size,
                             bool nosharpturn,
                             bool verbose,
                             string bppfile,
                             string bppfileindex,
                             bool pfonly,
                             float bppcutoff,
			                 string forestfile,
                             bool mea,
                             float MEA_gamma,
                             string MEA_file_index,
                             bool MEA_bpseq,
                             bool ThreshKnot,
                             float ThreshKnot_threshold,
                             string ThreshKnot_file_index)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose),
      bpp_file(bppfile),
      bpp_file_index(bppfileindex),
      pf_only(pfonly),
      bpp_cutoff(bppcutoff),
      forest_file(forestfile), 
      mea_(mea),
      gamma(MEA_gamma),
      mea_file_index(MEA_file_index),
      bpseq(MEA_bpseq),
      threshknot_(ThreshKnot),
      threshknot_threshold(ThreshKnot_threshold),
      threshknot_file_index(ThreshKnot_file_index){
#ifdef lpv
        initialize();
#else
        initialize();
        initialize_cachesingle();
#endif

}

int main(int argc, char** argv){

    struct timeval total_starttime, total_endtime;
    gettimeofday(&total_starttime, NULL);

    int beamsize = 100;
    bool sharpturn = false;
    bool is_verbose = false;
    string bpp_file;
    string bpp_prefix;
    bool pf_only = false;
    float bpp_cutoff = 0.0;
    string forest_file;

    float MEA_gamma = 3.0;
    bool mea = false;
    bool MEA_bpseq = false;
    string MEA_prefix;
    float ThreshKnot_threshold = 0.3;
    bool ThreshKnot = false;
    string ThresKnot_prefix;


    if (argc > 1) {
        beamsize = atoi(argv[1]);
        sharpturn = atoi(argv[2]) == 1;
        is_verbose = atoi(argv[3]) == 1;
        bpp_file = argv[4];
        bpp_prefix = argv[5];
        pf_only = atoi(argv[6]) == 1;
        bpp_cutoff = atof(argv[7]);
    	forest_file = argv[8];
        mea = atoi(argv[9]) == 1;
        MEA_gamma = atof(argv[10]);
        ThreshKnot = atoi(argv[11]) == 1;
        ThreshKnot_threshold = atof(argv[12]);
        ThresKnot_prefix = argv[13];
        MEA_prefix = argv[14];
        MEA_bpseq = atoi(argv[15]) == 1;
    }


    if (is_verbose) printf("beam size: %d\n", beamsize);

    // variables for decoding
    int num=0, total_len = 0;
    unsigned long long total_states = 0;
    double total_score = .0;
    double total_time = .0;

    int seq_index = 0;
    string bpp_file_index = "";
    string ThreshKnot_file_index = "";
    string MEA_file_index = "";


    std::vector<std::string> MSA_;

    for (string seq; getline(cin, seq);) {
        if (seq.length() == 0)
            continue;

        if (seq[0] == ';' || seq[0] == '>') {
            // printf("%s\n", seq.c_str());
            if (!bpp_file.empty()) {
                FILE *fptr = fopen(bpp_file.c_str(), "a"); 
                if (fptr == NULL) { 
                    printf("Could not open file!\n"); 
                    return 0; 
                }
                fprintf(fptr, "%s\n", seq.c_str());
                fclose(fptr); 
            }
            continue;
        }

        seq_index ++;
        if (!bpp_prefix.empty()) bpp_file_index = bpp_prefix + to_string(seq_index);

        if (!ThresKnot_prefix.empty()) ThreshKnot_file_index = ThresKnot_prefix + to_string(seq_index);

        if (!MEA_prefix.empty()) MEA_file_index = MEA_prefix + to_string(seq_index);
        
        // convert to uppercase
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

        // convert T to U
        replace(seq.begin(), seq.end(), 'T', 'U');

        MSA_.push_back(seq);

        // printf("%s\n", seq.c_str());

    }

    auto n_seq = MSA_.size();
    auto MSA_seq_length = MSA_[0].size();
    auto ribo_ = get_ribosum(MSA_, n_seq, MSA_seq_length);
    auto pscore = init_pscores_only(MSA_seq_length);
    vector<float> smart_gap;
    vector<vector<int>> a2s_fast, s5_fast, s3_fast, SS_fast;
    a2s_prepare_is(MSA_, n_seq, MSA_seq_length, a2s_fast, s5_fast, s3_fast, SS_fast, smart_gap);
    BeamCKYParser parser(beamsize, !sharpturn, is_verbose, bpp_file, bpp_file_index, pf_only, bpp_cutoff, forest_file, mea, MEA_gamma, MEA_file_index, MEA_bpseq, ThreshKnot, ThreshKnot_threshold, ThreshKnot_file_index);
    parser.parse_alifold(MSA_, a2s_fast, pscore, s5_fast, s3_fast, SS_fast, ribo_, smart_gap);

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;

    return 0;
}
