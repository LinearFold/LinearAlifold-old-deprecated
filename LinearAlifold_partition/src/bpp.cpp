
#include <stdio.h> 
#include <set>
#include <algorithm>
#include "linearalifold_p.h"
#include "Utils/ribo.h"

using namespace std;

void BeamCKYParser::output_to_file(string file_name, const char * type) {
    if(!file_name.empty()) {
        printf("Outputing base pairing probability matrix to %s...\n", file_name.c_str()); 
        FILE *fptr = fopen(file_name.c_str(), type); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }

        // int turn = no_sharp_turn?3:0;
        for (int i = 1; i <= seq_length; i++) {
            for (int j = i + turn + 1; j <= seq_length; j++) {
                pair<int, int> key = make_pair(i,j);
                auto got = Pij.find(key);
                if (got != Pij.end()){
                    fprintf(fptr, "%d %d %.4e\n", i, j, got->second);
                }
            }
        }
        fprintf(fptr, "\n");
        fclose(fptr); 
        printf("Done!\n"); 
    }

    return;
}

void BeamCKYParser::output_to_file_MEA_threshknot_bpseq(string file_name, const char * type, map<int,int>& pairs, string & seq) {

    int i,j;
    char nuc;
    if(!file_name.empty()) {
        printf("Outputing base pairs in bpseq format to %s...\n", file_name.c_str()); 
        FILE *fptr = fopen(file_name.c_str(), type); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }

        for (int i = 1; i <= seq_length; i++) {
            if (pairs.find(i) != pairs.end()){
                j = pairs[i];
            }
            else{
                j = 0;
            }
            nuc = seq[i-1];
            fprintf(fptr, "%d %c %d\n", i, nuc, j);
        }

        fprintf(fptr, "\n");
        fclose(fptr); 
        printf("Done!\n"); 
    }
    else{
        for (int i = 1; i <= seq_length; i++) {
            if (pairs.find(i) != pairs.end()){
                j = pairs[i];
            }
            else{
                j = 0;
            }
            nuc = seq[i-1];
            printf("%d %c %d\n", i, nuc, j);
        }
        printf("\n");
    }

}

void BeamCKYParser::cal_PairProb(State& viterbi, vector<vector<ribo_state>> & pscore) {
    
    double kTn = double(kT) * MSA.size();
    
    for(int j=0; j<seq_length; j++){
        for(auto &item : bestP[j]){
            int i = item.first;
            State state = item.second;
            pf_type temp_prob_inside = state.alpha + state.beta - viterbi.alpha - pscore[i][j].ribo_score/kTn;
            if (temp_prob_inside > pf_type(-9.91152)) {
                pf_type prob = Fast_Exp(temp_prob_inside);
                if(prob > pf_type(1.0)) prob = pf_type(1.0);
                if(prob < pf_type(bpp_cutoff)) continue;
                Pij[make_pair(i+1, j+1)] = prob;
            }
        }
    }

    // -o mode: output to a single file with user specified name;
    // bpp matrices for different sequences are separated with empty lines
    if (!bpp_file.empty()){
        output_to_file(bpp_file, "a");
    } 

    // -prefix mode: output to multiple files with user specified prefix;
    else if (!bpp_file_index.empty()) {
        output_to_file(bpp_file_index, "w");
    }
    return;
}


string BeamCKYParser::back_trace(const int i, const int j, const vector<vector<int> >& back_pointer){

    if (i>j) return "";
    if (back_pointer[i][j] == -1){
        if (i == j) return ".";
        else return "." + back_trace(i+1,j, back_pointer);
    }else if (back_pointer[i][j] != 0){
        int k = back_pointer[i][j];
        assert(k + 1 > 0 && k + 1 <= seq_length);
        string temp;
        if (k == j) temp = "";
        else temp = back_trace(k+1,j, back_pointer);
        return "(" + back_trace(i+1,k-1, back_pointer) + ")" + temp;
    }
    assert(false);
    return "";
}

map<int, int> BeamCKYParser::get_pairs(string & structure){
    map<int, int> pairs;
    stack<int> s;
    int index = 1;
    int pre_index = 0;
    for (auto & elem : structure){
        if (elem == '(') s.push(index);
        else if(elem == ')'){
            pre_index = s.top();
            pairs[pre_index] = index;
            pairs[index] = pre_index;
            s.pop();
        }
        index++;
    }
    return pairs;
}

void BeamCKYParser::ThreshKnot(string & seq){
    
    map<int, pf_type> rowprob;
    vector<tuple<int, int, pf_type> > prob_list;

    map<int, int> pairs;
    set<int> visited;

    for(auto& pij : Pij){
        auto i = pij.first.first; //index starts from 1
        auto j = pij.first.second; 
        auto score = pij.second;

        if (score < threshknot_threshold) continue;

        prob_list.push_back(make_tuple(i,j,score));

        rowprob[i] = max(rowprob[i], score);
        rowprob[j] = max(rowprob[j], score);

    }

    for(auto& elem : prob_list){

        auto i = std::get<0>(elem);
        auto j = std::get<1>(elem);
        auto score =  std::get<2>(elem);

        if (score == rowprob[i] && score == rowprob[j]){

            if ((visited.find(i) != visited.end()) || (visited.find(j) != visited.end())) continue;
            visited.insert(i);
            visited.insert(j);

            pairs[i] = j;
            pairs[j] = i;
        }
    }

    // fprintf(stdout, "%s\n", seq.c_str());
    output_to_file_MEA_threshknot_bpseq(threshknot_file_index, "w", pairs, seq);
}


void BeamCKYParser::PairProb_MEA(string & seq) {
    
    vector<vector<pf_type> > OPT;
    OPT.resize(seq_length);

    for (int i = 0; i < seq_length; ++i) OPT[i].resize(seq_length);

    vector<vector<pf_type>> P;
    P.resize(seq_length);

    for (int i = 0; i < seq_length; ++i) P[i].resize(seq_length);

    vector<vector<int> > back_pointer;
    back_pointer.resize(seq_length);

    for (int i = 0; i < seq_length; ++i) back_pointer[i].resize(seq_length);

    vector<vector<int>> paired;
    paired.resize(seq_length);

    vector<pf_type> Q;
    for (int i = 0; i < seq_length; ++i) Q.push_back(pf_type(1.0));

    for(auto& pij : Pij){
        auto i = pij.first.first-1;
        auto j = pij.first.second-1;
        auto score = pij.second;

        P[i][j] = score;

        paired[i].push_back(j);
        Q[i] -= score;
        Q[j] -= score;
    }

    for (int i = 0; i < seq_length; ++i) std::sort (paired[i].begin(), paired[i].end());
    for (int l = 0; l< seq_length; l++){
        for (int i = 0; i<seq_length - l; i++){
            int j = i + l;
            if (i == j){
                OPT[i][j] = Q[i];
                back_pointer[i][j] = -1;
                continue;
            }
            OPT[i][j] = OPT[i][i] + OPT[i+1][j];
            back_pointer[i][j] = -1;
            for (int k : paired[i]){
                if (k>j) break;
                pf_type temp_OPT_k1_j;
                if (k<j) temp_OPT_k1_j = OPT[k+1][j];
                else temp_OPT_k1_j = pf_type(0.);
                auto temp_score = 2 * gamma * P[i][k] + OPT[i+1][k-1] + temp_OPT_k1_j;
                if (OPT[i][j] < temp_score){
                    OPT[i][j] = temp_score;
                    back_pointer[i][j] = k;
                }
            }
        }
    }

    auto structure = back_trace(0,seq_length-1, back_pointer);

    if (!bpseq){
        if(!mea_file_index.empty()) {
            FILE *fptr = fopen(mea_file_index.c_str(), "w"); 
            if (fptr == NULL) { 
                printf("Could not open file!\n"); 
                return; 
            }
            // fprintf(fptr, "%s\n", seq.c_str());
            fprintf(fptr, "%s\n\n", structure.c_str());
        }

        else{
            // printf("%s\n", seq.c_str());
            printf("%s\n\n", structure.c_str());
        }
    }

    else{
        auto pairs = get_pairs(structure);
        output_to_file_MEA_threshknot_bpseq(mea_file_index, "w", pairs, seq);
    }
}


void BeamCKYParser::outside_alifold(vector<int> next_pair[], vector<vector<ribo_state>> & pscore, vector<float> & smart_gap, float smart_gap_threshold, vector<vector<int>> & next_position){
      
    struct timeval bpp_starttime, bpp_endtime;
    gettimeofday(&bpp_starttime, NULL);

    bestC[seq_length-1].beta = 0.0;

    int n_seq = MSA.size();

    double kTn = double(kT) * n_seq;



    // from right to left
    value_type newscore;
    for(int j = seq_length-1; j > 0; --j) {

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


        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
            Fast_LogPlusEquals(beamstepC.beta, (bestC[j+1].beta));
            }
        }
    
        // beam of M
        {
            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
                    Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta);
                }
            }
        }

        // beam of M2
        {
            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

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

                            Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta);

                        }
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(state.beta, beamstepM[i].beta);
            }
        }

        // beam of P
        {  

            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;

                auto & s5_i = s5_fast[i];
                auto & SS_i = SS_fast[i];

                auto & a2s_i = a2s_fast[i];
                auto & a2s_j = a2s_fast[j];

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
                                    // helix

                                    newscore = 0;
                                    int nucp, nucq, nucp1, nucq_1, nuci_1, nucj1, u1_local, u2_local, type;

                                    for (int s = 0; s < MSA.size(); s++){

                                        nucp = SS_p[s];
                                        nucq = SS_q[s];

                                        nucp1 = s3_p[s];
                                        nucq_1 = s5_q[s];
                                        nuci_1 = s5_i[s];
                                        nucj1 = s3_j[s];

                                        type = NUM_TO_PAIR(nucp, nucq);

                                        newscore += -v_score_single_alifold(0, 0, type, tt2[s], nucp1, nucq_1, nuci_1, nucj1);
                                    }

                                    Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore/kTn);

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
                                    Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore/kTn);



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

                    int new_nuci_1, new_nuci, new_nucj, new_nucj1;


                    auto & s5_i = s5_fast[i];
                    auto & SS_i = SS_fast[i];
                    auto & SS_j = SS_fast[j];
                    auto & s3_j = s3_fast[j];                    

                    newscore = 0;

                    for (int s = 0; s < MSA.size(); s++){
                        new_nuci_1 = ((i - 1) > -1)? s5_i[s] : -1;
                        new_nuci = SS_i[s];
                        new_nucj = SS_j[s];
                        new_nucj1 = (j + 1) < seq_length? s3_j[s] : -1;
                        newscore += (- v_score_M1(-1, -1, -1, new_nuci_1, new_nuci, new_nucj, new_nucj1, -1)); // no position information needed
                    }


                    Fast_LogPlusEquals(state.beta, beamstepM[i].beta + newscore/kTn);


                }
                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {

                    newscore = 0;

                    int new_nuci_1, new_nuci, new_nucj, new_nucj1;

                    auto & s5_i = s5_fast[i];
                    auto & SS_i = SS_fast[i];
                    auto & SS_j = SS_fast[j];
                    auto & s3_j = s3_fast[j];


                    for (int s = 0; s < MSA.size(); s++){

                         new_nuci_1 = s5_i[s]; //TODO, need to check boundary?
                         new_nuci = SS_i[s];
                         new_nucj = SS_j[s];
                         new_nucj1 = (j + 1) < seq_length ? s3_j[s] : -1; //TODO, need to check boundary? it may diff. from RNAalifold

                        newscore += (- v_score_M1(-1, -1, -1, new_nuci_1, new_nuci, new_nucj, new_nucj1, -1));
                    }                    

                    pf_type m1_alpha = newscore/kTn;
                    pf_type m1_plus_P_alpha = state.alpha + m1_alpha;
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(state.beta, (beamstepM2[newi].beta + m_state.alpha + m1_alpha));
                        Fast_LogPlusEquals(m_state.beta, (beamstepM2[newi].beta + m1_plus_P_alpha));
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {

                        int new_nuck, new_nuck1, new_nucj, new_nucj1;

                        auto & s5_i = s5_fast[i];
                        auto & SS_i = SS_fast[i];
                        auto & SS_j = SS_fast[j];
                        auto & s3_j = s3_fast[j];
                        auto & a2s_i = a2s_fast[i];
                        auto & a2s_j = a2s_fast[j];
                        auto & a2s_seq_length_1 = a2s_fast[seq_length-1];


                        newscore = 0;

                        for (int s = 0; s < MSA.size(); s++){
                            new_nuck = (a2s_i[s] > 0) ? s5_i[s] : -1; //external.c line 1165, weird
                            new_nuck1 = SS_i[s];
                            new_nucj = SS_j[s];
                            new_nucj1 = (a2s_j[s] < a2s_seq_length_1[s]) ? s3_j[s] : -1; //external.c line 1165, weird

                            newscore += (- v_score_external_paired(-1, -1, new_nuck, new_nuck1, new_nucj, new_nucj1, -1));
                        }      

                        pf_type external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + newscore/kTn;

                        Fast_LogPlusEquals(bestC[k].beta, state.alpha + external_paired_alpha_plus_beamstepC_beta);
                        Fast_LogPlusEquals(state.beta, bestC[k].alpha + external_paired_alpha_plus_beamstepC_beta);
                    } else {

                        int new_nuck1, new_nucj, new_nucj1;

                        auto & SS_i = SS_fast[i];
                        auto & SS_j = SS_fast[j];
                        auto & s3_j = s3_fast[j];
                        auto & a2s_j = a2s_fast[j];
                        auto & a2s_seq_length_1 = a2s_fast[seq_length-1];

                        newscore = 0;

                        for (int s = 0; s < MSA.size(); s++){
                            new_nuck1 = SS_i[s];
                            new_nucj = SS_j[s];
                            new_nucj1 = (a2s_j[s] < a2s_seq_length_1[s]) ? s3_j[s] : -1; //external.c line 1165, weird

                            newscore += - v_score_external_paired(0, j, -1, new_nuck1, new_nucj, new_nucj1, -1);

                        }

                        Fast_LogPlusEquals(state.beta, (beamstepC.beta + newscore/kTn));

                    }
                }
                
                if (pscore[i][j].ribo_score == std::numeric_limits<int>::lowest()){
                    pscore[i][j].ribo_score = make_pscores_ij(SS_fast[i], SS_fast[j], ribo);
                }

                state.beta = state.beta + pscore[i][j].ribo_score / kTn;
            }

        }

        // beam of Multi
        {
            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;

                auto & SS_i = SS_fast[i];
                auto & s3_i = s3_fast[i];
                auto & s5_j = s5_fast[j];
                auto & SS_j = SS_fast[j];

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
                        Fast_LogPlusEquals(state.beta, (bestMulti[jnext][i].beta));
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

                    Fast_LogPlusEquals(state.beta, beamstepP[i].beta + newscore/kTn);
                }
            }
        }

    }  // end of for-loo j


    gettimeofday(&bpp_endtime, NULL);
    double bpp_elapsed_time = bpp_endtime.tv_sec - bpp_starttime.tv_sec + (bpp_endtime.tv_usec-bpp_starttime.tv_usec)/1000000.0;
    if(is_verbose) fprintf(stdout,"Base Pairing Probabilities Calculation Time: %.2f seconds.\n", bpp_elapsed_time);
    fflush(stdout);

    return;
}

